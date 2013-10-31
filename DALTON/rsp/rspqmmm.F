!
!...   Copyright (c) 2013 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2013 (2013), see http://daltonprogram.org"
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
C**************************************************************************
C  /* Deck qmmmlno */
      SUBROUTINE QMMMLNO(NOSIM,BOVECS,CREF,CMO,XINDX,UDV,DV,UDVTR,
     &                  DVTR,EVECS,WRK,LWRK)
C
C  JK, Nov.08
C
C  Purpose:  Calculate SCF/DFT E^2 contribution from a surrounding
C            polarizable MM medium to an orbital trial vector. New QMMM
C            code.
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "mmtimes.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "infdim.h"

      DIMENSION BOVECS(*), CMO(*), XINDX(*), UDV(NASHDI,NASHDI)
      DIMENSION UDVTR(N2ASHX), DVTR(*), EVECS(KZYVAR,*)
      DIMENSION WRK(*), DV(*), CREF(*)

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LOCDEB
      LOGICAL FNDLAB
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)

      LOCDEB = .FALSE.

      CALL QENTER('QMMMLNO')

      IF (MMTIME) THEN 
        DTIME = SECOND()
        BTIME = SECOND()
      ENDIF

      IF (NNZAL .EQ. 0) THEN ! no pol. and nothing to do here
        CALL QEXIT('QMMMLNO')
        RETURN
      ENDIF

      IF (LGSPOL) THEN
        CALL QEXIT('QMMMLNO')
        RETURN
      ENDIF

      IF  ( ((NASHT.EQ.0) .AND. (TRPLET))  ) THEN
        CALL QEXIT('QMMMLNO')
        RETURN
      ENDIF


      KUCMO   = 1
      KUBO    = KUCMO   + NORBT*NBAST
      KRXYO   = KUBO    + NOSIM*N2ORBX
      IF (TRPLET) THEN
        KRXYOT  = KRXYO   + NOSIM*N2ORBX
        KEFIELD = KRXYOT  + NOSIM*N2ORBX
      ELSE
        KEFIELD = KRXYO   + NOSIM*N2ORBX
      ENDIF
      KINDMOM = KEFIELD + 3*NNZAL
      IF (MMMAT) THEN
        KINVMAT = KINDMOM + 3*NNZAL
        KLAST    = KINVMAT + 3*NNZAL*(3*NNZAL+1)/2
      ELSE IF (MMITER) THEN
        KLAST    = KINDMOM + 3*NNZAL
      ENDIF
      LWRK1 = LWRK - KLAST + 1      

      IF (LWRK1 .LT. 0) CALL ERRWRK('QMMMLNO',-KLAST,LWRK)

      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KUBO),NOSIM*N2ORBX)
      CALL DZERO(WRK(KRXYO),NOSIM*N2ORBX)
      IF (TRPLET) CALL DZERO(WRK(KRXYOT),NOSIM*N2ORBX)
      CALL DZERO(WRK(KEFIELD),3*NNZAL)
      CALL DZERO(WRK(KINDMOM),3*NNZAL)
      IF (MMMAT) CALL DZERO(WRK(KINVMAT),3*NNZAL*(3*NNZAL+1)/2)

C     Read the relay matrix from file

      N = 3*NNZAL 

      IF (MMMAT) THEN
        LUQMMM = -1
        IF (LUQMMM .LT. 0) THEN
          CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
        ENDIF
        REWIND(LUQMMM)

        IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
          CALL READT(LUQMMM,N*(N+1)/2,WRK(KINVMAT))
        ELSE
          CALL QUIT('Problem reading the matrix from the QMMMIM file.')
        ENDIF

        CALL GPCLOSE(LUQMMM,'KEEP')

        IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
          WRITE(LUPRI,*) ' The classical response matrix is'//
     &                   ' read from file'
          DO I = 1, N*(N+1)/2
            WRITE(LUPRI,*) WRK(KINVMAT+I-1)
          END DO
        ENDIF
      ENDIF

      CALL UPKCMO(CMO,WRK(KUCMO))

      IF (NOSIM.GT.0) THEN
         CALL RSPZYM(NOSIM,BOVECS,WRK(KUBO))
         CALL DSCAL(NOSIM*N2ORBX,-1.0D0,WRK(KUBO),1)
      END IF

      IF (MMTIME) THEN 
        DTIME = SECOND() - DTIME
        TMMLNO0 = TMMLNO0 + DTIME
      ENDIF
C     In contrast to the old code, the zero-order contribution, i.e.
C     - sum_a mu_ind_a t_a_pq has been included since we in 
C     sirius do not subtract this from the fock/ks-operator. We here only 
C     need to calculate the explicit electronic contribution, i.e.
C     - sum_a mu_ind_a[delta rho] t_q_pa, where delta rho is the first order
C     response/change in the density due to the applied perturbation. 
C
      DO 200 IOSIM = 1,NOSIM

         IF (MMTIME) DTIME = SECOND()
         JRXYO = KRXYO + (IOSIM-1)*N2ORBX
         IF (TRPLET) JRXYOT = KRXYOT + (IOSIM-1)*N2ORBX
         JUBO  = KUBO  + (IOSIM-1)*N2ORBX

C-------------------------------------------------------
C 1. Construct vector containing expectation value of electric 
C    field due to first order change in the density. Length of vector 
C    is NNZAL = no. of pol. sites.
C-------------------------------------------------------

         CALL DZERO(WRK(KEFIELD),3*NNZAL)

#if defined(VAR_MPI)
         IF (NODTOT .GE. 1) THEN
            CALL QMMMLNO_M1(UDV,UDVTR,WRK(KEFIELD),WRK(KUCMO),
     *                      WRK(JUBO),WRK(KLAST),LWRK1)
         ELSE
#endif
           LRI = 1             ! group-index in vector
           DO 201 J=1,MMCENT

             IF (ZEROAL(J) .EQ. -1) GOTO 201

             CALL MMLNO_ITER1(J,WRK(KEFIELD+LRI-1+0),
     *                        WRK(KEFIELD+LRI-1+1),WRK(KEFIELD+LRI-1+2),
     *                        UDV,UDVTR,WRK(KUCMO),WRK(JUBO),
     *                        WRK(KLAST),LWRK1,IPRRSP)
             LRI = LRI + 3
             GOTO 201
  201      CONTINUE        
#if defined(VAR_MPI)
         ENDIF
#endif
         NDIM = 3*NNZAL
         NTOTI = MAX(NDIM,1)

         IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
           WRITE(LUPRI,*) 'Transformed electric fields'
           CALL OUTPUT(WRK(KEFIELD),1,NDIM,1,1,NDIM,1,1,LUPRI)
           WRITE(LUPRI,*)
         ENDIF
         IF (MMTIME) THEN
           DTIME = SECOND() - DTIME
           TMMLNO1 = TMMLNO1 + DTIME
         ENDIF
C        Calculate the induced dipoles corresponding to this field

        IF (MMTIME) DTIME = SECOND()
         IF (MMMAT) THEN
           CALL DSPMV('L', NDIM, ONE, WRK(KINVMAT), WRK(KEFIELD), 1,
     &                ZERO, WRK(KINDMOM), 1)
         ELSE IF (MMITER) THEN
           IOPT = 2 ! Do not read from file any previuos induced moments
           CALL F2QMMM(WRK(KEFIELD),NNZAL,WRK(KINDMOM),WRK(KLAST),LWRK1,
     *                 IOPT,IPQMMM)
         ENDIF

         IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
           WRITE(LUPRI,*) 'Transformed induced dipole moments' 
           CALL OUTPUT(WRK(KINDMOM),1,N,1,1,N,1,1,LUPRI)
           WRITE(LUPRI,*)
         ENDIF

C        Now we need the integrals again transformed to the MO basis. 
C        These could have been stored
C        before, but we choose to do all of this integral-direct.
         IF (MMTIME) THEN
           DTIME = SECOND() - DTIME
           TMMLNO2 = TMMLNO2 + DTIME
         ENDIF

         IF (MMTIME) DTIME = SECOND()
         
         IF (TRPLET) THEN
            LRXYO = JRXYOT
         ELSE
            LRXYO = JRXYO
         ENDIF

#if defined(VAR_MPI)
         IF (NODTOT .GE. 1) THEN
            CALL QMMMLNO_M2(WRK(LRXYO),WRK(KUCMO),WRK(KINDMOM),
     *                      WRK(KLAST),LWRK1)
         ELSE
#endif
            LRI = 1
            DO 202 J=1,MMCENT
               IF (ZEROAL(J) .EQ. -1) GOTO 202
               CALL MMLNO_ITER2(J,WRK(KINDMOM+LRI-1),
     *                        WRK(KINDMOM+LRI),WRK(KINDMOM+LRI+1),
     *                        WRK(KUCMO),WRK(LRXYO),
     *                        WRK(KLAST),LWRK1,IPRRSP)
               LRI = LRI + 3
 202        CONTINUE 
#if defined(VAR_MPI)
         ENDIF
#endif
         IF (MMTIME) THEN
           DTIME = SECOND() - DTIME
           TMMLNO3 = TMMLNO3 + DTIME
         ENDIF
 200  CONTINUE

C--------------------------------------------
C 3.  Add effective operators to response
C--------------------------------------------

      IF (MMTIME) DTIME = SECOND()

      IF (TRPLET) THEN
         CALL SLVSOR(.TRUE.,.FALSE.,NOSIM,UDVTR,EVECS(1,1),WRK(KRXYO))
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYOT))
      ELSE
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYO))
      ENDIF

      IF (MMTIME) THEN
        DTIME   = SECOND() - DTIME
        TMMLNO3 = TMMLNO3  + DTIME
        BTIME   = SECOND() - BTIME
        TMMLNO  = TMMLNO   + BTIME
        TMMRSP  = TMMRSP   + BTIME
      ENDIF

      CALL QEXIT('QMMMLNO')
      RETURN
      END 
C-----------------------------------------------------------------------
C
C  /* Deck qmmmqro */
      SUBROUTINE QMMMQRO(VEC1,VEC2,ETRS,XINDX,ZYM1,ZYM2,
     &                 UDV,WRK,LWRK,KZYVR,KZYV1,KZYV2,
     &                 IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     &                 ISPIN0,ISPIN1,ISPIN2)
C
C  JK, Dec. 08
C
C  Purpose:  Calculate SCF/DFT E^3 contribution from a surrounding
C            polarizable MM medium to an orbital trial vector. New 
C            QMMM code
C
#include "implicit.h"
#include "dummy.h"
#include "maxorb.h"
#include "inforb.h"
#include "infdim.h"
#include "infinp.h"
#include "infvar.h"
#include "infrsp.h"
#include "infpri.h"
#include "rspprp.h"
#include "infcr.h"
#include "inftap.h"
#include "qrinf.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "mmtimes.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)

      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(LWRK),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2)
      DIMENSION MJWOP(2,MAXWOP,8)
      LOGICAL LCON, LORB, LREF
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LOCDEB, FNDLAB
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMQRO')

      IF (MMTIME) DTIME = SECOND()
      IF (MMTIME) BTIME = SECOND()

      LOCDEB = .FALSE.

      IF (NNZAL .EQ. 0) THEN ! no pol. and nothing to do here
        CALL QEXIT('QMMMQRO')
        RETURN
      ENDIF

      IF (LGSPOL) THEN
        CALL QEXIT('QMMMQRO')
        RETURN
      ENDIF

      IF (NASHT.GE.1) CALL QUIT('ONLY CLOSED SHELL FOR QR-QM/MM')

      IF (TRPLET) CALL QUIT('NO TRPLET IN QR-QM/MM YET')

      KCREF   = 1
      KTRES   = KCREF   + NCREF
      KUCMO   = KTRES   + N2ORBX
      KTLMA   = KUCMO   + NORBT*NBAST
      KTLMB   = KTLMA   + N2ORBX
      KTRMO   = KTLMB   + N2ORBX
      KUTR    = KTRMO   + NNORBX
      IF (MMMAT) THEN
        KINVMAT = KUTR    + N2ORBX
        KEF1    = KINVMAT + 3*NNZAL*(3*NNZAL+1)/2
      ELSE
        KEF1    = KUTR    + N2ORBX
      ENDIF
      KEF2    = KEF1    + 3*NNZAL
      KIND1   = KEF2    + 3*NNZAL
      KIND2   = KIND1   + 3*NNZAL
      KWRK    = KIND2   + 3*NNZAL
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QMMMQRO',-KWRK,LWRK1)

      CALL DZERO(WRK(KCREF),NCREF)
      CALL DZERO(WRK(KTRES),N2ORBX)
      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KTLMA),N2ORBX)
      CALL DZERO(WRK(KTLMB),N2ORBX)
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      IF (MMMAT) CALL DZERO(WRK(KINVMAT),3*NNZAL*(3*NNZAL+1)/2)
      CALL DZERO(WRK(KEF1),3*NNZAL)
      CALL DZERO(WRK(KEF2),3*NNZAL)
      CALL DZERO(WRK(KIND1),3*NNZAL)
      CALL DZERO(WRK(KIND2),3*NNZAL)

      NSIM = 1

C     We assume no symmetry in the DFT/MM calculations although we have kept 
C     some symmetry options below ...

      ISYMT = 1

C     Get the reference state

      CALL GETREF(WRK(KCREF),MZCONF(1))

C     Unpack the response vectors

      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)

C     Unpack symmetry blocked CMO

      CALL UPKCMO(CMO,WRK(KUCMO))

      N = 3*NNZAL

      IF (MMMAT) THEN

        LUQMMM = -1
        IF (LUQMMM .LT. 0) THEN
          CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
        ENDIF
        REWIND(LUQMMM)

        IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
          CALL READT(LUQMMM,N*(N+1)/2,WRK(KINVMAT))
        ELSE
          CALL QUIT('Problem reading the matrix from the QMMMIM file.')
        ENDIF

        CALL GPCLOSE(LUQMMM,'KEEP')

        IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
          WRITE(LUPRI,*) ' The classical response matrix is'//
     &                   ' read from file'
          DO I = 1, N*(N+1)/2
            WRITE(LUPRI,*) WRK(KINVMAT+I-1)
          END DO
        ENDIF

      ENDIF

      KMAT = KWRK
      KLAST = KMAT + 3*NNBASX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMQRO',-KLAST,LWRK2)

      IF (MMTIME) THEN 
        DTIME = SECOND() - DTIME
        TMMQRO0 = TMMQRO0 + DTIME
      ENDIF

      IF (MMTIME) DTIME = SECOND()

#if defined(VAR_MPI)
      IF (NODTOT .GE. 1) THEN
         CALL QMMMQRO_M1(WRK(KEF1),WRK(KEF2),UDV,WRK(KUCMO),
     &                   ISYMT,ISYMV1,ISYMV2,ZYM1,ZYM2,WRK(KLAST),LWRK2)
      ELSE
#endif
        LRI = 1                ! group-index in vector
        DO 100 J = 1,MMCENT 
          IF (ZEROAL(J) .EQ. -1) GOTO 100
          CALL MMQRO_ITER1(J,WRK(KEF1+LRI-1),WRK(KEF2+LRI-1),UDV,
     *                     WRK(KUCMO),ISYMT,ISYMV1,ISYMV2,ZYM1,ZYM2,
     *                     WRK(KLAST),LWRK2,IPRRSP)
          LRI = LRI + 3
  100   CONTINUE
#if defined(VAR_MPI)
      ENDIF
#endif
      IF (MMTIME) THEN 
        DTIME = SECOND() - DTIME
        TMMQRO1 = TMMQRO1 + DTIME
      ENDIF
C     Calculate the induced dipoles corresponding to the 
C     (transformed) fields 1 and 2
      IF (MMTIME) DTIME = SECOND()

      NDIM = 3*NNZAL
      NTOTI = MAX(NDIM,1)

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
        WRITE(LUPRI,*) 'Induced field 1'
        CALL OUTPUT(WRK(KEF1),1,NDIM,1,1,NDIM,1,1,LUPRI)
        WRITE(LUPRI,*)
      ENDIF

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
        WRITE(LUPRI,*) 'Induced field 2'
        CALL OUTPUT(WRK(KEF2),1,NDIM,1,1,NDIM,1,1,LUPRI)
        WRITE(LUPRI,*)
      ENDIF

      IF (MMMAT) THEN
        CALL DSPMV('L', NDIM, ONE, WRK(KINVMAT), WRK(KEF1), 1, ZERO,
     &             WRK(KIND1), 1)
      ELSE IF (MMITER) THEN
        IOPT = 2 ! Do not read from file any previuos induced moments
        CALL F2QMMM(WRK(KEF1),NNZAL,WRK(KIND1),WRK(KLAST),LWRK2,
     *              IOPT,IPQMMM)
      ENDIF

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
        WRITE(LUPRI,*) 'Transformed induced dipole moments field 1'
        CALL OUTPUT(WRK(KIND1),1,NDIM,1,1,NDIM,1,1,LUPRI)
        WRITE(LUPRI,*)
      ENDIF

      IF (MMMAT) THEN
        CALL DSPMV('L', NDIM, ONE, WRK(KINVMAT), WRK(KEF2), 1, ZERO, 
     &             WRK(KIND2), 1)
      ELSE IF (MMITER) THEN
        IOPT = 2 ! Do not read from file any previuos induced moments
        CALL F2QMMM(WRK(KEF2),NNZAL,WRK(KIND2),WRK(KLAST),LWRK2,
     *              IOPT,IPQMMM)
      ENDIF

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
        WRITE(LUPRI,*) 'Transformed induced dipole moments field 2'
        CALL OUTPUT(WRK(KIND2),1,NDIM,1,1,NDIM,1,1,LUPRI)
        WRITE(LUPRI,*)
      ENDIF

      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMQRO2 = TMMQRO2 + DTIME
      ENDIF
C     Now we need the integrals again transformed to the MO basis 
C     and also one-index transformed integrals. 
C     These could have been stored
C     before, but we choose to do all of this integral-direct.
      IF (MMTIME) DTIME = SECOND()

#if defined(VAR_MPI)
      IF (NODTOT .GE. 1) THEN
         CALL QMMMQRO_M2(WRK(KIND1),WRK(KIND2),WRK(KUCMO),
     *                   WRK(KTRES),ISYMT,ISYMV2,ZYM2,
     *                   WRK(KLAST),LWRK2)
      ELSE
#endif
        LRI = 1                ! group-index in vector
        DO 101 J = 1,MMCENT 
           IF (ZEROAL(J) .EQ. -1) GOTO 101
           CALL MMQRO_ITER2(J,WRK(KIND1+LRI-1),WRK(KIND2+LRI-1),
     *                   WRK(KTRES),WRK(KUCMO),ISYMT,ISYMV2,ZYM2,
     *                   WRK(KLAST),LWRK2,IPRRSP)
           LRI = LRI + 3
 101    CONTINUE
#if defined(VAR_MPI)
      ENDIF
#endif
      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMQRO3 = TMMQRO3 + DTIME
      ENDIF

      IF (MMTIME) DTIME = SECOND()
C       Make the gradient
C
C     / <0| [qj ,TRES] |0> \
C     |          0         |
C     | <0| [qj+,TRES] |0> |
C      \         0         /
C
      ISYMDN = 1
      OVLAP  = ONE
      JSPIN = 0
      ISYMV  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      NZYVEC = NCREF
      NZCVEC = NCREF

      CALL RSP1GR(NSIM,KZYVR,IDUM,JSPIN,IGRSYM,JSPIN,ISYMV,ETRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,WRK(KTRES),
     *            XINDX,MJWOP,WRK(KWRK),LWRK1,LORB,LCON,LREF)
C
      IF (MMTIME) THEN
        DTIME   = SECOND() - DTIME
        TMMQRO4 = TMMQRO4  + DTIME
        BTIME   = SECOND() - BTIME
        TMMQRO  = TMMQRO   + BTIME
        TMMRSP  = TMMRSP   + BTIME
      ENDIF
      CALL QEXIT('QMMMQRO')
      RETURN
      END
C
C-------------------------------------------------------------------------------
C  /* Deck qmmmcro */
      SUBROUTINE QMMMCRO(VEC1,VEC2,VEC3,ETRS,XINDX,
     &                 ZYM1,ZYM2,ZYM3,UDV,WRK,LWRK,
     &                 KZYVR,KZYV1,KZYV2,KZYV3,
     &                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,MJWOP)
C
C  JMO, Apr. 09
C
C  Purpose:  Calculate SCF/DFT E^4 contribution from a surrounding
C            polarizable MM medium to an orbital trial vector. New 
C            QMMM code
C
#include "implicit.h"
#include "maxorb.h"
#include "inforb.h"
#include "infdim.h"
#include "infinp.h"
#include "infvar.h"
#include "infrsp.h"
#include "infpri.h"
#include "rspprp.h"
#include "infcr.h"
#include "inftap.h"
#include "qrinf.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "mmtimes.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "dummy.h"
#include "infpar.h"

      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*),WRK(LWRK),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3)
      DIMENSION MJWOP(2,MAXWOP,8)
      LOGICAL LCON, LORB, LREF
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LOCDEB, FNDLAB
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMCRO')

      IF (MMTIME) DTIME = SECOND()
      IF (MMTIME) BTIME = SECOND()

      LOCDEB = .FALSE.

      IF (NNZAL .EQ. 0) THEN ! no pol. and nothing to do here
        CALL QEXIT('QMMMCRO')
        RETURN
      ENDIF

      IF (NASHT.GE.1) CALL QUIT('ONLY CLOSED SHELL FOR CR-QM/MM')

      IF (TRPLET) CALL QUIT('NO TRPLET IN CR-QM/MM YET')

      KCREF   = 1
      KTRES   = KCREF   + NCREF
      KUCMO   = KTRES   + N2ORBX
      KTLMA   = KUCMO   + NORBT*NBAST
      KTLMB   = KTLMA   + N2ORBX
      KTLMC   = KTLMB   + N2ORBX
C Change here KTLMC
      KTRMO   = KTLMC   + N2ORBX
      KUTR    = KTRMO   + NNORBX
      IF (MMMAT) THEN
        KINVMAT = KUTR    + N2ORBX
        KEF1    = KINVMAT + 3*NNZAL*(3*NNZAL+1)/2
      ELSE
        KEF1    = KUTR    + N2ORBX
      ENDIF
      KEF2    = KEF1    + 3*NNZAL
      KEF3    = KEF2    + 3*NNZAL
      KIND1   = KEF3    + 3*NNZAL
      KIND2   = KIND1   + 3*NNZAL
      KIND3   = KIND2   + 3*NNZAL
      KWRK    = KIND3   + 3*NNZAL
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QMMMCRO',-KWRK,LWRK1)

      CALL DZERO(WRK(KCREF),NCREF)
      CALL DZERO(WRK(KTRES),N2ORBX)
      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KTLMA),N2ORBX)
      CALL DZERO(WRK(KTLMB),N2ORBX)
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      IF (MMMAT) CALL DZERO(WRK(KINVMAT),3*NNZAL*(3*NNZAL+1)/2)
      CALL DZERO(WRK(KEF1),3*NNZAL)
      CALL DZERO(WRK(KEF2),3*NNZAL)
      CALL DZERO(WRK(KEF3),3*NNZAL)
      CALL DZERO(WRK(KIND1),3*NNZAL)
      CALL DZERO(WRK(KIND2),3*NNZAL)
      CALL DZERO(WRK(KIND3),3*NNZAL)

      NSIM = 1

C     We assume no symmetry in the DFT/MM calculations although we have kept 
C     some symmetry options below ...

      ISYMT = 1

C     Get the reference state

      CALL GETREF(WRK(KCREF),MZCONF(1))

C     Unpack the response vectors

      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)
      CALL GTZYMT(NSIM,VEC3,KZYV3,ISYMV3,ZYM3,MJWOP)

C     Unpack symmetry blocked CMO

      CALL UPKCMO(CMO,WRK(KUCMO))

      N = 3*NNZAL

      IF (MMMAT) THEN

        LUQMMM = -1
        IF (LUQMMM .LT. 0) THEN
          CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
        ENDIF
        REWIND(LUQMMM)

        IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
          CALL READT(LUQMMM,N*(N+1)/2,WRK(KINVMAT))
        ELSE
          CALL QUIT('Problem reading the matrix from the QMMMIM file.')
        ENDIF

        CALL GPCLOSE(LUQMMM,'KEEP')

        IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
          WRITE(LUPRI,*) ' The classical response matrix is'//
     &                   ' read from file'
          DO I = 1, N*(N+1)/2
            WRITE(LUPRI,*) WRK(KINVMAT+I-1)
          END DO
        ENDIF

      ENDIF 

      KMAT = KWRK
      KLAST = KMAT + 3*NNBASX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMCRO',-KLAST,LWRK2)

      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMCRO0 = TMMCRO0 + DTIME
      ENDIF

      IF (MMTIME) DTIME = SECOND()

#if defined(VAR_MPI)
      IF (NODTOT .GE. 1) THEN
         CALL QMMMCRO_M1(WRK(KEF1),WRK(KEF2),WRK(KEF3),UDV,WRK(KUCMO),
     &                   ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM1,ZYM2,ZYM3,
     &                   WRK(KLAST),LWRK2)
      ELSE
#endif
        LRI = 1                ! group-index in vector
        DO 100 J = 1,MMCENT 
           IF (ZEROAL(J) .EQ. -1) GOTO 100
           CALL MMCRO_ITER1(J,WRK(KEF1+LRI-1),WRK(KEF2+LRI-1),
     *                   WRK(KEF3+LRI-1),UDV,WRK(KUCMO),
     *                   ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM1,ZYM2,ZYM3,
     *                   WRK(KLAST),LWRK2,IPRRSP)
           LRI = LRI + 3
 100    CONTINUE
#if defined(VAR_MPI)
      ENDIF
#endif
      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMCRO1 = TMMCRO1 + DTIME
      ENDIF

      IF (MMTIME) DTIME = SECOND()
C     Calculate the induced dipoles corresponding to the 
C     (transformed) fields 1, 2 and 3 either by matrix inversion or 
C     iteratively.

      NDIM = 3*NNZAL
      NTOTI = MAX(NDIM,1)

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
        WRITE(LUPRI,*) 'Induced field 1'
        CALL OUTPUT(WRK(KEF1),1,NDIM,1,1,NDIM,1,1,LUPRI)
        WRITE(LUPRI,*)
      ENDIF

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
        WRITE(LUPRI,*) 'Induced field 2'
        CALL OUTPUT(WRK(KEF2),1,NDIM,1,1,NDIM,1,1,LUPRI)
        WRITE(LUPRI,*)
      ENDIF

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
        WRITE(LUPRI,*) 'Induced field 3'
        CALL OUTPUT(WRK(KEF3),1,NDIM,1,1,NDIM,1,1,LUPRI)
        WRITE(LUPRI,*)
      ENDIF

      IF (MMMAT) THEN
        CALL DSPMV('L', NDIM, ONE, WRK(KINVMAT), WRK(KEF1), 1, ZERO, 
     &             WRK(KIND1), 1)
      ELSE IF (MMITER) THEN
        IOPT = 2 ! Do not read from file any previuos induced moments
        CALL F2QMMM(WRK(KEF1),NNZAL,WRK(KIND1),WRK(KLAST),LWRK2,
     *              IOPT,IPQMMM)
      ENDIF

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
        WRITE(LUPRI,*) 'Transformed induced dipole moments field 1'
        CALL OUTPUT(WRK(KIND1),1,NDIM,1,1,NDIM,1,1,LUPRI)
        WRITE(LUPRI,*)
      ENDIF

      IF (MMMAT) THEN
        CALL DSPMV('L', NDIM, ONE, WRK(KINVMAT), WRK(KEF2), 1, ZERO, 
     &             WRK(KIND2), 1)
      ELSE IF (MMITER) THEN
        IOPT = 2 ! Do not read from file any previuos induced moments
        CALL F2QMMM(WRK(KEF2),NNZAL,WRK(KIND2),WRK(KLAST),LWRK2,
     *              IOPT,IPQMMM)
      ENDIF

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
        WRITE(LUPRI,*) 'Transformed induced dipole moments field 2'
        CALL OUTPUT(WRK(KIND2),1,NDIM,1,1,NDIM,1,1,LUPRI)
        WRITE(LUPRI,*)
      ENDIF

      IF (MMMAT) THEN
        CALL DSPMV('L', NDIM, ONE, WRK(KINVMAT), WRK(KEF3), 1, ZERO, 
     &             WRK(KIND3), 1)
      ELSE IF (MMITER) THEN
        IOPT = 2 ! Do not read from file any previuos induced moments
        CALL F2QMMM(WRK(KEF3),NNZAL,WRK(KIND3),WRK(KLAST),LWRK2,
     *              IOPT,IPQMMM)
      ENDIF

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
        WRITE(LUPRI,*) 'Transformed induced dipole moments field 3'
        CALL OUTPUT(WRK(KIND3),1,NDIM,1,1,NDIM,1,1,LUPRI)
        WRITE(LUPRI,*)
      ENDIF

C     Now we need the integrals again transformed to the MO basis 
C     and also one-index transformed integrals. 
C     These could have been stored
C     before, but we choose to do all of this integral-direct.

      IF (MMTIME) THEN 
        DTIME = SECOND() - DTIME
        TMMCRO2 = TMMCRO2 + DTIME
      ENDIF

      IF (MMTIME) DTIME = SECOND()

#if defined(VAR_MPI)
      IF (NODTOT .GE. 1) THEN
         CALL QMMMCRO_M2(WRK(KIND1),WRK(KIND2),WRK(KIND3),WRK(KUCMO),
     &        WRK(KTRES),ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM2,ZYM3,
     &        WRK(KLAST),LWRK2)
      ELSE
#endif
         LRI = 1                ! group-index in vector
         DO 101 J = 1,MMCENT 
            IF (ZEROAL(J) .EQ. -1) GOTO 101
            CALL MMCRO_ITER2(J,WRK(KIND1+LRI-1),WRK(KIND2+LRI-1),
     *                       WRK(KIND3+LRI-1),WRK(KTRES),WRK(KUCMO),
     *                       ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM2,ZYM3,
     *                       WRK(KLAST),LWRK2,IPRRSP)
            LRI = LRI + 3
 101     CONTINUE
#if defined(VAR_MPI)
      ENDIF
#endif
      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMCRO3 = TMMCRO3 + DTIME
      ENDIF

      IF (MMTIME) DTIME = SECOND()

C       Make the gradient
C
C     / <0| [qj ,TRES] |0> \
C     |          0         |
C     | <0| [qj+,TRES] |0> |
C      \         0         /
C
      ISYMDN = 1
      OVLAP  = ONE
      JSPIN = 0
      ISYMV  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      NZYVEC = NCREF
      NZCVEC = NCREF

      CALL RSP1GR(NSIM,KZYVR,IDUM,JSPIN,IGRSYM,JSPIN,ISYMV,ETRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,WRK(KTRES),
     *            XINDX,MJWOP,WRK(KLAST),LWRK2,LORB,LCON,LREF)

      IF (MMTIME) THEN
        DTIME   = SECOND() - DTIME
        TMMCRO4 = TMMCRO4  + DTIME
        BTIME   = SECOND() - BTIME
        TMMCRO  = TMMCRO   + BTIME
        TMMRSP  = TMMRSP   + BTIME
      ENDIF

      CALL QEXIT('QMMMCRO')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck mmlno_iter1 */
      SUBROUTINE MMLNO_ITER1(J,EEX,EEY,EEZ,UDV,UDVTR,UCMO,UBO,
     *                       WRK,LWRK,IPRTMP)
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
#include "infdim.h"

      DIMENSION UDV(NASHDI,NASHDI),UDVTR(N2ASHX)
      DIMENSION UCMO(NORBT*NBAST),UBO(N2ORBX),WRK(LWRK)
 
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('MMLNO_ITER1')

      EEX = 0.0D0
      EEY = 0.0D0
      EEZ = 0.0D0

      DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *        (MMCORD(2,J)-QMCOM(2))**2 +
     *        (MMCORD(3,J)-QMCOM(3))**2
      DIST = SQRT(DIST2)

      IF (DIST .GT. RCUTMM) THEN
         CALL QEXIT('MMLNO_ITER1')
         RETURN
      ENDIF

C     Backup diporg. We use diporg to transfer coordinates to int.
C     program.

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      KMAT    = 1
      KURXAC  = KMAT  + 3*NNBASX
      KTRMO   = KURXAC  + N2ASHX
      KUTR    = KTRMO   + NNORBX
      KUTRX   = KUTR    + N2ORBX
      KLAST   = KUTRX   + N2ORBX

      LWRK2 = LWRK - KLAST + 1
      IF (LWRK2 .LT. 0) CALL ERRWRK('MMLNO_ITER1',-KLAST,LWRK2)

      CALL DZERO(WRK(KMAT),3*NNBASX)
      
      KPATOM = 0
      NOCOMP = 3
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      DIPORG(1) = MMCORD(1,J)
      DIPORG(2) = MMCORD(2,J)
      DIPORG(3) = MMCORD(3,J)

      RUNQM3 = .TRUE.
      CALL GET1IN(WRK(KMAT),'NEFIELD',NOCOMP,WRK(KLAST),
     &               LWRK2,LABINT,INTREP,INTADR,J,TOFILE,
     &               KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
      RUNQM3 = .FALSE.

      IF (QMDAMP) THEN
         IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
         ENDIF
         IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
            DIST = 9.99D+99
            MHIT = 0
            DO 123 M=1,NUCIND
               DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &                 (DIPORG(2)-CORD(2,M))**2 +
     &                 (DIPORG(3)-CORD(3,M))**2
               IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
               ENDIF
 123        CONTINUE
         ELSE IF (IDAMP .EQ. 2) THEN
            DIST = (DIPORG(1)-QMCOM(1))**2 +
     &             (DIPORG(2)-QMCOM(2))**2 +
     &             (DIPORG(3)-QMCOM(3))**2
         ENDIF
         DIST = SQRT(DIST)

         IF (IDAMP .EQ. 3) THEN
            IF (IPOLTP .EQ. 2) THEN
               TEMPI =  D3I*(POLMM(1,J)+POLMM(4,J)+POLMM(6,J))
            ELSE IF (IPOLTP .EQ. 1) THEN
               IF (IPOLTP .EQ. 1) TEMPI = POLIMM(J)
            ENDIF
            TEMP = (TEMPI*QMPOL(MHIT))**(D6I)
            SIJ = 2.1304*DIST/TEMP
            DFACT = 1.0D0 - (1.0D0+SIJ+0.50D0*SIJ*SIJ)*exp(-SIJ)
         ELSE
            DFACT = (1-exp(-ADAMP*DIST))**3
         ENDIF
         CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
      ENDIF

C     x-component of (induced) electric field
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      CALL DZERO(WRK(KUTRX),N2ORBX)
      CALL DZERO(WRK(KURXAC),N2ASHX)

      CALL UTHU(WRK(KMAT),WRK(KTRMO),UCMO,WRK(KLAST),
     &             NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
      CALL ONEXH1(UBO,WRK(KUTR),WRK(KUTRX))

      IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),WRK(KURXAC))
      IF (TRPLET) THEN
         EEX = SLVTLM(UDVTR,WRK(KURXAC),WRK(KUTRX),TAC)
      ELSE
         EEX = SLVQLM(UDV,WRK(KURXAC),WRK(KUTRX),TAC)
      ENDIF

C     y-component of (induced) electric field
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      CALL DZERO(WRK(KUTRX),N2ORBX)
      CALL DZERO(WRK(KURXAC),N2ASHX)

      CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),UCMO,
     &             WRK(KLAST),NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
      CALL ONEXH1(UBO,WRK(KUTR),WRK(KUTRX))
      IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),WRK(KURXAC))
      IF (TRPLET) THEN
         EEY = SLVTLM(UDVTR,WRK(KURXAC),WRK(KUTRX),TAC)
      ELSE
         EEY = SLVQLM(UDV,WRK(KURXAC),WRK(KUTRX),TAC)
      ENDIF

C     z-component of (induced) electric field
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      CALL DZERO(WRK(KUTRX),N2ORBX)
      CALL DZERO(WRK(KURXAC),N2ASHX)

      CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),UCMO,
     &          WRK(KLAST),NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
      CALL ONEXH1(UBO,WRK(KUTR),WRK(KUTRX))
      IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),WRK(KURXAC))
      IF (TRPLET) THEN
         EEZ = SLVTLM(UDVTR,WRK(KURXAC),WRK(KUTRX),TAC)
      ELSE
         EEZ = SLVQLM(UDV,WRK(KURXAC),WRK(KUTRX),TAC)
      ENDIF

C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('MMLNO_ITER1')
      RETURN
      END


C-------------------------------------------------------------------------------
C  /* Deck mmlno_iter2 */
      SUBROUTINE MMLNO_ITER2(J,XIND,YIND,ZIND,UCMO,RXYO,
     *                       WRK,LWRK,IPRTMP)
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
#include "infdim.h"

      DIMENSION RXYO(*),UCMO(*),WRK(LWRK)
 
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('MMLNO_ITER2')

      DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *        (MMCORD(2,J)-QMCOM(2))**2 +
     *        (MMCORD(3,J)-QMCOM(3))**2
      DIST = SQRT(DIST2)

      IF (DIST .GT. RCUTMM) THEN
         CALL QEXIT('MMLNO_ITER2')
         RETURN
      ENDIF

C     Backup the dipole origin

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      KTRMO   = 1
      KUTR    = KTRMO   + NNORBX
      KMAT    = KUTR    + N2ORBX
      KLAST   = KMAT    + 3*NNBASX
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('MMLNO_ITER2',-KLAST,LWRK2)

      CALL DZERO(WRK(KMAT),3*NNBASX)

      KPATOM = 0
      NOCOMP = 3
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      DIPORG(1) = MMCORD(1,J)
      DIPORG(2) = MMCORD(2,J)
      DIPORG(3) = MMCORD(3,J)

      RUNQM3 = .TRUE.
      CALL GET1IN(WRK(KMAT),'NEFIELD',NOCOMP,WRK(KLAST),
     &               LWRK2,LABINT,INTREP,INTADR,J,TOFILE,
     &               KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
      RUNQM3 = .FALSE.

      IF (QMDAMP) THEN
         IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
         ENDIF
         IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) )THEN
            DIST = 9.99D+99
            MHIT = 0
            DO 124 M=1,NUCIND
               DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &                 (DIPORG(2)-CORD(2,M))**2 +
     &                 (DIPORG(3)-CORD(3,M))**2
               IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
               ENDIF
 124        CONTINUE
         ELSE IF (IDAMP .EQ. 2) THEN
            DIST = (DIPORG(1)-QMCOM(1))**2 +
     &             (DIPORG(2)-QMCOM(2))**2 +
     &             (DIPORG(3)-QMCOM(3))**2
         ENDIF
         DIST = SQRT(DIST)

         IF (IDAMP .EQ. 3) THEN
            IF (IPOLTP .EQ. 2) THEN
               TEMPI =  D3I*(POLMM(1,J)+POLMM(4,J)+POLMM(6,J))
            ELSE IF (IPOLTP .EQ. 1) THEN
               IF (IPOLTP .EQ. 1) TEMPI = POLIMM(J)
            ENDIF
            TEMP = (TEMPI*QMPOL(MHIT))**(D6I)
            SIJ = 2.1304*DIST/TEMP
            DFACT = 1.0D0 - (1.0D0+SIJ+0.50D0*SIJ*SIJ)*exp(-SIJ)
         ELSE
            DFACT = (1-exp(-ADAMP*DIST))**3
         ENDIF
         CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
      ENDIF
      
C     x-component of electric field 
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      CALL UTHU(WRK(KMAT),WRK(KTRMO),UCMO,WRK(KLAST),
     &             NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
      FACx = -XIND
      CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,RXYO,1)

C     y-component of electric field
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),UCMO,
     &             WRK(KLAST),NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
      FACy = -YIND
      CALL DAXPY(N2ORBX,FACy,WRK(KUTR),1,RXYO,1)

C     z-component of (induced) electric field
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),UCMO,
     &             WRK(KLAST),NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
      FACz = -ZIND
      CALL DAXPY(N2ORBX,FACz,WRK(KUTR),1,RXYO,1)

C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('MMLNO_ITER2')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck mmqro_iter1 */
      SUBROUTINE MMQRO_ITER1(J,EF1,EF2,UDV,UCMO,ISYMT,ISYMV1,ISYMV2,
     &                       ZYM1,ZYM2,WRK,LWRK,IPRTMP)
C
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "infdim.h"

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)

      DIMENSION WRK(LWRK), EF1(3), EF2(3)
      DIMENSION UDV(NASHDI,NASHDI),UCMO(*),ZYM1(*),ZYM2(*)

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('MMQRO_ITER1')

      DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *        (MMCORD(2,J)-QMCOM(2))**2 +
     *        (MMCORD(3,J)-QMCOM(3))**2
      DIST = SQRT(DIST2)

      IF (DIST .GT. RCUTMM) THEN
         CALL QEXIT('MMQRO_ITER1')
         RETURN
      ENDIF

C     Backup the dipole origin

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      KTRMO   = 1
      KUTR    = KTRMO   + NNORBX
      KMAT    = KUTR    + N2ORBX
      KTLMA   = KMAT    + 3*NNBASX
      KTLMB   = KTLMA   + N2ORBX
      KLAST   = KTLMB   + N2ORBX
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('MMQRO_ITER1',-KLAST,LWRK2)

      CALL DZERO(WRK(KMAT),3*NNBASX)

      KPATOM = 0
      NOCOMP = 3
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      DIPORG(1) = MMCORD(1,J)
      DIPORG(2) = MMCORD(2,J)
      DIPORG(3) = MMCORD(3,J)

      RUNQM3 = .TRUE.
      CALL GET1IN(WRK(KMAT),'NEFIELD',NOCOMP,WRK(KLAST),
     &                LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
      RUNQM3 = .FALSE.

      IF (QMDAMP) THEN
         IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
         ENDIF
         IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
            DIST = 9.99D+99
            MHIT = 0
            DO 125 M=1,NUCIND
               DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &              (DIPORG(2)-CORD(2,M))**2 +
     &              (DIPORG(3)-CORD(3,M))**2
               IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
               ENDIF
 125        CONTINUE
         ELSE IF (IDAMP .EQ. 2) THEN
            DIST = (DIPORG(1)-QMCOM(1))**2 +
     &           (DIPORG(2)-QMCOM(2))**2 +
     &           (DIPORG(3)-QMCOM(3))**2
         ENDIF
         DIST = SQRT(DIST)
         
         IF (IDAMP .EQ. 3) THEN
            IF (IPOLTP .EQ. 2) THEN
               TEMPI =  D3I*(POLMM(1,J)+POLMM(4,J)+POLMM(6,J))
            ELSE IF (IPOLTP .EQ. 1) THEN
               IF (IPOLTP .EQ. 1) TEMPI = POLIMM(J)
            ENDIF
            TEMP = (TEMPI*QMPOL(MHIT))**(D6I)
            SIJ = 2.1304*DIST/TEMP
            DFACT = 1.0D0 - (1.0D0+SIJ+0.50D0*SIJ*SIJ)*exp(-SIJ)
         ELSE
            DFACT = (1-exp(-ADAMP*DIST))**3
         ENDIF
         CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
      ENDIF

C     x-component
      F1=ZERO
      F2=ZERO
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      CALL UTHU(WRK(KMAT),WRK(KTRMO),UCMO,WRK(KLAST),
     &              NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
C
      IF (ISYMT.EQ.ISYMV1) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMQRO')
         EF1(1) = F1
      ELSE
         EF1(1) = ZERO
      ENDIF 

      IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
         CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMQRO')
         EF2(1) = F2
      ELSE
         EF2(1) = ZERO
      ENDIF

C     y-component
      F1=ZERO
      F2=ZERO
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      
      CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),UCMO,
     &              WRK(KLAST),NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

      IF (ISYMT.EQ.ISYMV1) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMQRO')
         EF1(2) = F1
      ELSE
         EF1(2) = ZERO
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
         CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMQRO')
         EF2(2) = F2
      ELSE
         EF2(2) = ZERO
      ENDIF

C         z-component
      F1=ZERO
      F2=ZERO
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),UCMO,
     &     WRK(KLAST),NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

      IF (ISYMT.EQ.ISYMV1) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMQRO')
         EF1(3) = F1
      ELSE
         EF1(3) = ZERO
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
         CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMQRO')
         EF2(3) = F2
      ELSE
         EF2(3) = ZERO
      ENDIF

C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('MMQRO_ITER1')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck mmqro_iter2 */
      SUBROUTINE MMQRO_ITER2(J,EIND1,EIND2,TRES,UCMO,
     &                      ISYMT,ISYMV2,ZYM2,WRK,LWRK,IPRTMP)
C
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)

      DIMENSION WRK(LWRK), EIND1(3), EIND2(3)
      DIMENSION UCMO(*),ZYM2(*),TRES(*)

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('MMQRO_ITER2')

      DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *        (MMCORD(2,J)-QMCOM(2))**2 +
     *        (MMCORD(3,J)-QMCOM(3))**2
      DIST = SQRT(DIST2)

      IF (DIST .GT. RCUTMM) THEN
         CALL QEXIT('MMQRO_ITER2')
         RETURN
      ENDIF

C     Backup the dipole origin

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      KTRMO   = 1
      KUTR    = KTRMO   + NNORBX
      KMAT    = KUTR    + N2ORBX
      KTLMA   = KMAT    + 3*NNBASX
      KLAST   = KTLMA   + N2ORBX
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('MMQRO_ITER2',-KLAST,LWRK2)

      CALL DZERO(WRK(KMAT),3*NNBASX)

      KPATOM = 0
      NOCOMP = 3
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      DIPORG(1) = MMCORD(1,J)
      DIPORG(2) = MMCORD(2,J)
      DIPORG(3) = MMCORD(3,J)

      RUNQM3 = .TRUE.
      CALL GET1IN(WRK(KMAT),'NEFIELD',NOCOMP,WRK(KLAST),
     &                LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
      RUNQM3 = .FALSE.

      IF (QMDAMP) THEN
         IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
         ENDIF
         IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) )THEN
            DIST = 9.99D+99
            MHIT = 0
            DO 126 M=1,NUCIND
               DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &                  (DIPORG(2)-CORD(2,M))**2 +
     &                  (DIPORG(3)-CORD(3,M))**2
               IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
               ENDIF
 126        CONTINUE
         ELSE IF (IDAMP .EQ. 2) THEN
            DIST = (DIPORG(1)-QMCOM(1))**2 +
     &               (DIPORG(2)-QMCOM(2))**2 +
     &               (DIPORG(3)-QMCOM(3))**2
         ENDIF
         DIST = SQRT(DIST)
            
         IF (IDAMP .EQ. 3) THEN
            IF (IPOLTP .EQ. 2) THEN
               TEMPI =  D3I*(POLMM(1,J)+POLMM(4,J)+POLMM(6,J))
            ELSE IF (IPOLTP .EQ. 1) THEN
               IF (IPOLTP .EQ. 1) TEMPI = POLIMM(J)
            ENDIF
            TEMP = (TEMPI*QMPOL(MHIT))**(D6I)
            SIJ = 2.1304*DIST/TEMP
            DFACT = 1.0D0 - (1.0D0+SIJ+0.50D0*SIJ*SIJ)*exp(-SIJ)
         ELSE
            DFACT = (1-exp(-ADAMP*DIST))**3
         ENDIF
         CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
      ENDIF

C     x-component
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      CALL UTHU(WRK(KMAT),WRK(KTRMO),UCMO,WRK(KLAST),
     &              NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
C
      F1 = -EIND1(1)
      CALL DZERO(WRK(KTLMA),N2ORBX)
      CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
      CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,TRES,1)

      F2 = -0.50D0*EIND2(1)
      CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,TRES,1)

C     y-component
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),UCMO,
     &          WRK(KLAST),NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

      F1 = -EIND1(2)
      CALL DZERO(WRK(KTLMA),N2ORBX)
      CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
      CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,TRES,1)

      F2 = -0.50D0*EIND2(2)
      CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,TRES,1)

C     z-component
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),UCMO,
     &              WRK(KLAST),NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

      F1 = -EIND1(3)
      CALL DZERO(WRK(KTLMA),N2ORBX)
      CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
      CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,TRES,1)   

      F2 = -0.50D0*EIND2(3)
      CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,TRES,1)

C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('MMQRO_ITER2')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck mmcro_iter1 */
      SUBROUTINE MMCRO_ITER1(J,EF1,EF2,EF3,UDV,UCMO,
     &                       ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM1,ZYM2,ZYM3,
     &                       WRK,LWRK,IPRTMP)
C
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "infdim.h"

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)

      DIMENSION WRK(LWRK), EF1(3), EF2(3), EF3(3)
      DIMENSION UDV(NASHDI,NASHDI),UCMO(*),ZYM1(*),ZYM2(*),ZYM3(*)

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('MMCRO_ITER1')

      DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *        (MMCORD(2,J)-QMCOM(2))**2 +
     *        (MMCORD(3,J)-QMCOM(3))**2
      DIST = SQRT(DIST2)

      IF (DIST .GT. RCUTMM) THEN
         CALL QEXIT('MMCRO_ITER1')
         RETURN
      ENDIF

C     Backup the dipole origin

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      KTRMO   = 1
      KUTR    = KTRMO   + NNORBX
      KMAT    = KUTR    + N2ORBX
      KTLMA   = KMAT    + 3*NNBASX
      KTLMB   = KTLMA   + N2ORBX
      KLAST   = KTLMB   + N2ORBX
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('MMCRO_ITER1',-KLAST,LWRK2)

      CALL DZERO(WRK(KMAT),3*NNBASX)

      KPATOM = 0
      NOCOMP = 3
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      DIPORG(1) = MMCORD(1,J)
      DIPORG(2) = MMCORD(2,J)
      DIPORG(3) = MMCORD(3,J)

C     Get 1-electron integrals

      RUNQM3 = .TRUE.
      CALL GET1IN(WRK(KMAT),'NEFIELD',NOCOMP,WRK(KLAST),
     &                LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
      RUNQM3 = .FALSE.

      IF (QMDAMP) THEN
         IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
         ENDIF
         IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
            DIST = 9.99D+99
            MHIT = 0
            DO 123 M=1,NUCIND
               DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &                 (DIPORG(2)-CORD(2,M))**2 +
     &                 (DIPORG(3)-CORD(3,M))**2
               IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
               ENDIF
 123        CONTINUE
         ELSE IF (IDAMP .EQ. 2) THEN
            DIST = (DIPORG(1)-QMCOM(1))**2 +
     &             (DIPORG(2)-QMCOM(2))**2 +
     &             (DIPORG(3)-QMCOM(3))**2
         ENDIF
         DIST = SQRT(DIST)

         IF (IDAMP .EQ. 3) THEN
            IF (IPOLTP .EQ. 2) THEN
               TEMPI =  D3I*(POLMM(1,J)+POLMM(4,J)+POLMM(6,J))
            ELSE IF (IPOLTP .EQ. 1) THEN
               IF (IPOLTP .EQ. 1) TEMPI = POLIMM(J)
            ENDIF
            TEMP = (TEMPI*QMPOL(MHIT))**(D6I)
            SIJ = 2.1304*DIST/TEMP
            DFACT = 1.0D0 - (1.0D0+SIJ+0.50D0*SIJ*SIJ)*exp(-SIJ)
         ELSE
            DFACT = (1-exp(-ADAMP*DIST))**3
         ENDIF
         CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
      ENDIF

C     x-component
      F1=ZERO
      F2=ZERO
      F3=ZERO
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
C 
      CALL UTHU(WRK(KMAT),WRK(KTRMO),UCMO,WRK(KLAST),
     &          NBAST,NORBT)

C     Unpack

      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

      IF (ISYMT.EQ.ISYMV1) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)

C     1-index transformation

         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)

C     Evaluate expectation value

         CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMCRO')
         EF1(1) = F1
      ELSE
         EF1(1) = ZERO
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
         CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMCRO')
         EF2(1) = F2
      ELSE
         EF2(1) = ZERO
      END IF

      IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &              MULD2H(ISYMT,ISYMV1))
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV3,ZYM3,WRK(KTLMB),WRK(KTLMA),ISYMV3)
         CALL MELONE(WRK(KTLMA),1,UDV,ONE,F3,200,'QMMMCRO')
         EF3(1) = F3
      ELSE
         EF3(1) = ZERO
      END IF

C     y-component
      F1=ZERO
      F2=ZERO
      F3=ZERO
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      
C     AO to MO transformation

      CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),UCMO,WRK(KLAST),
     &          NBAST,NORBT)

C     Unpack from tri to full

      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

      IF (ISYMT.EQ.ISYMV1) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMCRO')
         EF1(2) = F1
      ELSE
         EF1(2) = ZERO
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
         CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMCRO')
         EF2(2) = F2
      ELSE
         EF2(2) = ZERO
      END IF

      IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &              MULD2H(ISYMT,ISYMV1))
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV3,ZYM3,WRK(KTLMB),WRK(KTLMA),ISYMV3)
         CALL MELONE(WRK(KTLMA),1,UDV,ONE,F3,200,'QMMMCRO')
         EF3(2) = F3
      ELSE
         EF3(2) = ZERO
      END IF

C     z-component
      F1=ZERO
      F2=ZERO
      F3=ZERO
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),UCMO,WRK(KLAST),
     &              NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

      IF (ISYMT.EQ.ISYMV1) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMCRO')
         EF1(3) = F1
      ELSE
         EF1(3) = ZERO
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
         CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMCRO')
         EF2(3) = F2
      ELSE
         EF2(3) = ZERO
      END IF
      
      IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &              MULD2H(ISYMT,ISYMV1))
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV3,ZYM3,WRK(KTLMB),WRK(KTLMA),ISYMV3)
         CALL MELONE(WRK(KTLMA),1,UDV,ONE,F3,200,'QMMMCRO')
         EF3(3) = F3
      ELSE
         EF3(3) = ZERO
      END IF

C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('MMCRO_ITER1')
      RETURN
      END


C-------------------------------------------------------------------------------
C  /* Deck mmcro_iter2 */
      SUBROUTINE MMCRO_ITER2(J,EIND1,EIND2,EIND3,TRES,UCMO,
     &                      ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM2,ZYM3,
     &                      WRK,LWRK,IPRTMP)
C
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"

      DIMENSION WRK(LWRK), EIND1(3), EIND2(3), EIND3(3)
      DIMENSION UCMO(*),ZYM2(*),ZYM3(*),TRES(*)

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('MMCRO_ITER2')

      DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *        (MMCORD(2,J)-QMCOM(2))**2 +
     *        (MMCORD(3,J)-QMCOM(3))**2
      DIST = SQRT(DIST2)

      IF (DIST .GT. RCUTMM) THEN
         CALL QEXIT('MMCRO_ITER2')
         RETURN
      ENDIF

C     Backup the dipole origin

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      KTRMO   = 1
      KUTR    = KTRMO   + NNORBX
      KMAT    = KUTR    + N2ORBX
      KTLMA   = KMAT    + 3*NNBASX
      KTLMB   = KTLMA   + N2ORBX
      KLAST   = KTLMB   + N2ORBX
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('MMCRO_ITER2',-KLAST,LWRK2)

      CALL DZERO(WRK(KMAT),3*NNBASX)

      KPATOM = 0
      NOCOMP = 3
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      DIPORG(1) = MMCORD(1,J)
      DIPORG(2) = MMCORD(2,J)
      DIPORG(3) = MMCORD(3,J)

      RUNQM3 = .TRUE.
      CALL GET1IN(WRK(KMAT),'NEFIELD',NOCOMP,WRK(KLAST),
     &            LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &            KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
      RUNQM3 = .FALSE.

      IF (QMDAMP) THEN
         IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
         ENDIF
         IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
            DIST = 9.99D+99
            MHIT = 0
            DO 124 M=1,NUCIND
               DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &                 (DIPORG(2)-CORD(2,M))**2 +
     &                 (DIPORG(3)-CORD(3,M))**2
               IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
               ENDIF
 124        CONTINUE
         ELSE IF (IDAMP .EQ. 2) THEN
            DIST = (DIPORG(1)-QMCOM(1))**2 +
     &             (DIPORG(2)-QMCOM(2))**2 +
     &             (DIPORG(3)-QMCOM(3))**2
         ENDIF
         DIST = SQRT(DIST)

         IF (IDAMP .EQ. 3) THEN
            IF (IPOLTP .EQ. 2) THEN
               TEMPI =  D3I*(POLMM(1,J)+POLMM(4,J)+POLMM(6,J))
            ELSE IF (IPOLTP .EQ. 1) THEN
               IF (IPOLTP .EQ. 1) TEMPI = POLIMM(J)
            ENDIF
            TEMP = (TEMPI*QMPOL(MHIT))**(D6I)
            SIJ = 2.1304*DIST/TEMP
            DFACT = 1.0D0 - (1.0D0+SIJ+0.50D0*SIJ*SIJ)*exp(-SIJ)
         ELSE
            DFACT = (1-exp(-ADAMP*DIST))**3
         ENDIF
         CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
      ENDIF
      
C     x-component
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      
      CALL UTHU(WRK(KMAT),WRK(KTRMO),UCMO,WRK(KLAST),
     &          NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

      IF (ISYMT.EQ.ISYMV1) THEN 
         F1 = -(1.0D0/2.0D0)*EIND1(1)
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV3,ZYM3,WRK(KTLMA),WRK(KTLMB),
     &              MULD2H(ISYMT,ISYMV2))
         CALL DAXPY(N2ORBX,F1,WRK(KTLMB),1,TRES,1)
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
         F2 = -(1.0D0/2.0D0)*EIND2(1)
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV3,ZYM3,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL DAXPY(N2ORBX,F2,WRK(KTLMA),1,TRES,1)
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
         F3 = -(1.0D0/6.0D0)*EIND3(1)
         CALL DAXPY(N2ORBX,F3,WRK(KUTR),1,TRES,1)
      ENDIF

C     y-component
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),UCMO,WRK(KLAST),
     &          NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

      IF (ISYMT.EQ.ISYMV1) THEN
         F1 = -(1.0D0/2.0D0)*EIND1(2)
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV3,ZYM3,WRK(KTLMA),WRK(KTLMB),
     &              MULD2H(ISYMT,ISYMV2))
         CALL DAXPY(N2ORBX,F1,WRK(KTLMB),1,TRES,1)
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
         F2 = -(1.0D0/2.0D0)*EIND2(2)
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV3,ZYM3,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL DAXPY(N2ORBX,F2,WRK(KTLMA),1,TRES,1)
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
         F3 = -(1.0D0/6.0D0)*EIND3(2)
         CALL DAXPY(N2ORBX,F3,WRK(KUTR),1,TRES,1)
      ENDIF

C     z-component
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      
      CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),UCMO,WRK(KLAST),
     &          NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
      
      IF (ISYMT.EQ.ISYMV1) THEN
         F1 = -(1.0D0/2.0D0)*EIND1(3)
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL DZERO(WRK(KTLMB),N2ORBX)
         CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL OITH1(ISYMV3,ZYM3,WRK(KTLMA),WRK(KTLMB),
     &              MULD2H(ISYMT,ISYMV2))
         CALL DAXPY(N2ORBX,F1,WRK(KTLMB),1,TRES,1)
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
         F2 = -(1.0D0/2.0D0)*EIND2(3)
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV3,ZYM3,WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL DAXPY(N2ORBX,F2,WRK(KTLMA),1,TRES,1)
      ENDIF

      IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
         F3 = -(1.0D0/6.0D0)*EIND3(3)
         CALL DAXPY(N2ORBX,F3,WRK(KUTR),1,TRES,1)
      ENDIF

C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('MMCRO_ITER2')
      RETURN
      END
#if defined(VAR_MPI)
C ************************************************
C Parallel QM/MM response routines (AHS 09-10)
C ************************************************
C  /* Deck qmmmlno_m1 */
      SUBROUTINE QMMMLNO_M1(UDV,UDVTR,CEFIELD,UCMO,
     *                      UBO,WRK,LWRK)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "gnrinf.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "infdim.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WRK(LWRK), UBO(*), UCMO(*)
      DIMENSION UDVTR(*), UDV(NASHDI,NASHDI), CEFIELD(3*NNZAL)

      CALL QENTER('QMMMLNO_M1')

      KEFIELD = 1
      KLAST   = KEFIELD + 3*NNZAL
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMLNO_M1',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = QMMMLNO_1_WORK
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

C     Send info
C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
      ENDIF

      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,NORB,'INTEGER',MASTER)

      IF (TRPLET) THEN
         CALL MPIXBCAST(UDVTR,N2ASHX,'DOUBLE',MASTER)
      ELSE
         CALL MPIXBCAST(UDV,NASHDI*NASHDI,'DOUBLE',MASTER)
      ENDIF

      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(UBO,N2ORBX,'DOUBLE',MASTER)

C     Do the work
      LRI = 1
      DO 20  J=1,MMCENT
         IWHO = -1
         IF (ZEROAL(J) .EQ. -1) GOTO 20
         CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
         CALL MPIXSEND(J, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
         LRI = LRI + 3
 20   CONTINUE

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data
      CALL DZERO(WRK(KEFIELD),3*NNZAL)
      CALL MPI_REDUCE(WRK(KEFIELD),CEFIELD,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMLNO_M1')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmlno_s1 */
      SUBROUTINE QMMMLNO_S1(WRK,LWRK,IPRTMP)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
#include "infdim.h"

      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMLNO_S1')

      QMMM = .TRUE.

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
      ENDIF

      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,NORB,'INTEGER',MASTER)
C    
      IF (TRPLET) THEN
         NFACTOR = N2ASHX
      ELSE
         NFACTOR = NASHDI*NASHDI
      ENDIF
C
      KUCMO   = 1
      JUBO    = KUCMO   + NORBT*NBAST
      KEFIELD = JUBO    + N2ORBX
      KUDV    = KEFIELD + 3*NNZAL 
      KLAST   = KUDV    + NFACTOR
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMLNO_S1',-KLAST,LWRK2)

      CALL MPIXBCAST(WRK(KUDV),NFACTOR,'DOUBLE',MASTER)      
      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(JUBO),N2ORBX,'DOUBLE',MASTER)

      CALL DZERO(WRK(KEFIELD),3*NNZAL)

 201  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN
         CALL MMLNO_ITER1(J,WRK(KEFIELD+LRI-1+0),
     *                    WRK(KEFIELD+LRI-1+1),WRK(KEFIELD+LRI-1+2),
     *                    WRK(KUDV),WRK(KUDV),WRK(KUCMO),WRK(JUBO),
     *                    WRK(KLAST),LWRK2,IPRRSP)
         GOTO 201
      ENDIF

C     No more integrals to calculate

      CALL MPI_REDUCE(WRK(KEFIELD),MPI_IN_PLACE,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      CALL QEXIT('QMMMLNO_S1')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmlno_m2 */
      SUBROUTINE QMMMLNO_M2(RXYO,UCMO,XINDMOM,WRK,LWRK)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "infdim.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WRK(LWRK)
      DIMENSION UCMO(*), RXYO(N2ORBX), XINDMOM(3*NNZAL)

      CALL QENTER('QMMMLNO_M2')

      KRXYO = 1
      KLAST = KRXYO + N2ORBX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMLNO_M2',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = QMMMLNO_2_WORK

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(XINDMOM,3*NNZAL,'DOUBLE',MASTER)

C     Do the work
      LRI = 1
      DO 20  J=1,MMCENT
         IWHO = -1
         IF (ZEROAL(J) .EQ. -1) GOTO 20
         CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
         CALL MPIXSEND(J, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
         LRI = LRI + 3
 20   CONTINUE

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data
      CALL DZERO(WRK(KRXYO),N2ORBX)
      CALL MPI_REDUCE(WRK(KRXYO),RXYO,N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMLNO_M2')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmlno_s2 */
      SUBROUTINE QMMMLNO_S2(WRK,LWRK,IPRTMP)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
#include "infdim.h"

      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMLNO_S2')

      QMMM = .TRUE.
      
      KINDMOM = 1
      KUCMO   = KINDMOM + 3*NNZAL
      JRXYO   = KUCMO   + NORBT*NBAST
      KLAST   = JRXYO   + N2ORBX
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMLNO_S2',-KLAST,LWRK2)

      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KINDMOM),3*NNZAL,'DOUBLE',MASTER)

      CALL DZERO(WRK(JRXYO),N2ORBX)

 201  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN
         CALL MMLNO_ITER2(J,WRK(KINDMOM+LRI-1),
     *                    WRK(KINDMOM+LRI),WRK(KINDMOM+LRI+1),
     *                    WRK(KUCMO),WRK(JRXYO),
     *                    WRK(KLAST),LWRK2,IPRTMP)
         GOTO 201
      ENDIF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(JRXYO),MPI_IN_PLACE,N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      CALL QEXIT('QMMMLNO_S2')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmqro_m1 */
      SUBROUTINE QMMMQRO_M1(EF1,EF2,UDV,UCMO,
     &                      ISYMT,ISYMV1,ISYMV2,ZYM1,ZYM2,WRK,LWRK)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "infdim.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WRK(LWRK)
      DIMENSION UDV(NASHDI,NASHDI),UCMO(*),ZYM1(*),ZYM2(*),EF1(*),EF2(*)

      CALL QENTER('QMMMQRO_M1')

      KEF1  = 1
      KEF2  = KEF1 + 3*NNZAL
      KLAST = KEF2 + 3*NNZAL
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMQRO_M1',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = QMMMQRO_1_WORK
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

C     Send info

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
      ENDIF

      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(UDV,NASHDI*NASHDI,'DOUBLE',MASTER)

C     QRO
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ZYM1,NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZYM2,NORBT*NORBT,'DOUBLE',MASTER)

C     Do the work
      LRI = 1
      DO 20  J=1,MMCENT
         IWHO = -1
         IF (ZEROAL(J) .EQ. -1) GOTO 20
         CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
         CALL MPIXSEND(J, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
         LRI = LRI + 3
 20   CONTINUE

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data
      CALL DZERO(WRK(KEF1),3*NNZAL)
      CALL MPI_REDUCE(WRK(KEF1),EF1,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
C
      CALL DZERO(WRK(KEF2),3*NNZAL)
      CALL MPI_REDUCE(WRK(KEF2),EF2,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMQRO_M1')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmqro_s1 */
      SUBROUTINE QMMMQRO_S1(WRK,LWRK,IPRTMP)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
#include "infdim.h"

      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMQRO_S1')

      QMMM = .TRUE.

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
      ENDIF
C
C
      KUDV    = 1
      KUCMO   = KUDV    + NASHDI*NASHDI
      KEF1    = KUCMO   + NORBT*NBAST
      KEF2    = KEF1    + 3*NNZAL
      KZYM1   = KEF2    + 3*NNZAL
      KZYM2   = KZYM1   + NORBT*NORBT
      KLAST   = KZYM2   + NORBT*NORBT
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMQRO_S1',-KLAST,LWRK2)

      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KUDV),NASHDI*NASHDI,'DOUBLE',MASTER)

C     QRO
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(WRK(KZYM1),NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KZYM2),NORBT*NORBT,'DOUBLE',MASTER)

      CALL DZERO(WRK(KEF1),3*NNZAL)
      CALL DZERO(WRK(KEF2),3*NNZAL)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN
         CALL MMQRO_ITER1(J,WRK(KEF1+LRI-1),WRK(KEF2+LRI-1),WRK(KUDV),
     *              WRK(KUCMO),ISYMT,ISYMV1,ISYMV2,WRK(KZYM1),
     *              WRK(KZYM2),WRK(KLAST),LWRK2,IPRTMP)
         GOTO 100
      ENDIF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KEF1),MPI_IN_PLACE,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(WRK(KEF2),MPI_IN_PLACE,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      CALL QEXIT('QMMMQRO_S1')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmqro_m2 */
      SUBROUTINE QMMMQRO_M2(ZIND1,ZIND2,UCMO,TRES,
     &                      ISYMT,ISYMV2,ZYM2,WRK,LWRK)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WRK(LWRK)
      DIMENSION UCMO(*),ZYM2(*),TRES(*),ZIND1(3*NNZAL),ZIND2(3*NNZAL)

      CALL QENTER('QMMMQRO_M2')

      KTRES = 1
      KLAST = KTRES + N2ORBX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMQRO_M2',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = QMMMQRO_2_WORK
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

C     Send info
      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)

C     QRO
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ZYM2,NORBT*NORBT,'DOUBLE',MASTER)

      CALL MPIXBCAST(ZIND1,3*NNZAL,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZIND2,3*NNZAL,'DOUBLE',MASTER)

C     Do the work
      LRI = 1
      DO 20  J=1,MMCENT
         IWHO = -1
         IF (ZEROAL(J) .EQ. -1) GOTO 20
         CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
         CALL MPIXSEND(J, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
         LRI = LRI + 3
 20   CONTINUE

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data
      CALL DZERO(WRK(KTRES),N2ORBX)
      CALL MPI_REDUCE(WRK(KTRES),TRES,N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMQRO_M2')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmqro_s2 */
      SUBROUTINE QMMMQRO_S2(WRK,LWRK,IPRTMP)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
#include "infdim.h"

      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMQRO_S2')

      QMMM = .TRUE.

      KUCMO   = 1
      KTRES   = KUCMO   + NORBT*NBAST
      KZYM2   = KTRES   + N2ORBX
      KIND1   = KZYM2   + NORBT*NORBT
      KIND2   = KIND1   + 3*NNZAL
      KLAST   = KIND2   + 3*NNZAL
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMQRO_S2',-KLAST,LWRK2)

      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)

C     QRO
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(WRK(KZYM2),NORBT*NORBT,'DOUBLE',MASTER)

      CALL MPIXBCAST(WRK(KIND1),3*NNZAL,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KIND2),3*NNZAL,'DOUBLE',MASTER)

      CALL DZERO(WRK(KTRES),N2ORBX)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN
         CALL MMQRO_ITER2(J,WRK(KIND1+LRI-1),WRK(KIND2+LRI-1),
     &                   WRK(KTRES),WRK(KUCMO),ISYMT,ISYMV2,WRK(KZYM2),
     &                   WRK(KLAST),LWRK2,IPRTMP)
         GOTO 100
      ENDIF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KTRES),MPI_IN_PLACE,N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      CALL QEXIT('QMMMQRO_S2')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmcro_m1 */
      SUBROUTINE QMMMCRO_M1(EF1,EF2,EF3,UDV,UCMO,
     &                      ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM1,ZYM2,ZYM3,
     &                      WRK,LWRK)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "infdim.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WRK(LWRK)
      DIMENSION UCMO(*),UDV(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*),EF1(*),EF2(*),EF3(*)

      CALL QENTER('QMMMCRO_M1')

      KEF1  = 1
      KEF2  = KEF1 + 3*NNZAL
      KEF3  = KEF2 + 3*NNZAL
      KLAST = KEF3 + 3*NNZAL
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMCRO_M1',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = QMMMCRO_1_WORK
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

C     Send info

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
      ENDIF

      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(UDV,NASHDI*NASHDI,'DOUBLE',MASTER)

C     CRO
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV3,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ZYM1,NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZYM2,NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZYM3,NORBT*NORBT,'DOUBLE',MASTER)

C     Do the work
      LRI = 1
      DO 20  J=1,MMCENT
         IWHO = -1
         IF (ZEROAL(J) .EQ. -1) GOTO 20
         CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
         CALL MPIXSEND(J, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
         LRI = LRI + 3
 20   CONTINUE

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data
      CALL DZERO(WRK(KEF1),3*NNZAL)
      CALL MPI_REDUCE(WRK(KEF1),EF1,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
C
      CALL DZERO(WRK(KEF2),3*NNZAL)
      CALL MPI_REDUCE(WRK(KEF2),EF2,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
C
      CALL DZERO(WRK(KEF3),3*NNZAL)
      CALL MPI_REDUCE(WRK(KEF3),EF3,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMCRO_M1')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmcro_s1 */
      SUBROUTINE QMMMCRO_S1(WRK,LWRK,IPRTMP)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
#include "infdim.h"

      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMCRO_S1')

      QMMM = .TRUE.

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
      ENDIF
C
      KUDV    = 1
      KUCMO   = KUDV    + NASHDI*NASHDI
      KEF1    = KUCMO   + NORBT*NBAST
      KEF2    = KEF1    + 3*NNZAL
      KEF3    = KEF2    + 3*NNZAL
      KZYM1   = KEF3    + 3*NNZAL
      KZYM2   = KZYM1   + NORBT*NORBT
      KZYM3   = KZYM2   + NORBT*NORBT
      KLAST   = KZYM3   + NORBT*NORBT
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMCRO_S1',-KLAST,LWRK2)

      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KUDV),NASHDI*NASHDI,'DOUBLE',MASTER)

C     CRO
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV3,1,'INTEGER',MASTER)
      CALL MPIXBCAST(WRK(KZYM1),NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KZYM2),NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KZYM3),NORBT*NORBT,'DOUBLE',MASTER)

      CALL DZERO(WRK(KEF1),3*NNZAL)
      CALL DZERO(WRK(KEF2),3*NNZAL)
      CALL DZERO(WRK(KEF3),3*NNZAL)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN
         CALL MMCRO_ITER1(J,WRK(KEF1+LRI-1),WRK(KEF2+LRI-1),
     *                   WRK(KEF3+LRI-1),WRK(KUDV),WRK(KUCMO),
     *                   ISYMT,ISYMV1,ISYMV2,ISYMV3,
     *                   WRK(KZYM1),WRK(KZYM2),WRK(KZYM3),
     *                   WRK(KLAST),LWRK2,IPRTMP)
         GOTO 100
      ENDIF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KEF1),MPI_IN_PLACE,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL MPI_REDUCE(WRK(KEF2),MPI_IN_PLACE,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL MPI_REDUCE(WRK(KEF3),MPI_IN_PLACE,3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMCRO_S1')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmcro_m2 */
      SUBROUTINE QMMMCRO_M2(ZIND1,ZIND2,ZIND3,UCMO,TRES,
     &                      ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM2,ZYM3,
     &                      WRK,LWRK)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WRK(LWRK)
      DIMENSION ZIND1(3*NNZAL),ZIND2(3*NNZAL),ZIND3(3*NNZAL)
      DIMENSION UCMO(*),TRES(*),ZYM2(*),ZYM3(*)

      CALL QENTER('QMMMCRO_M2')

      KTRES = 1
      KLAST = KTRES + N2ORBX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMCRO_M2',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = QMMMCRO_2_WORK
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

C     Send info

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
      ENDIF

      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)

C     CRO
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV3,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ZYM2,NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZYM3,NORBT*NORBT,'DOUBLE',MASTER)

      CALL MPIXBCAST(ZIND1,3*NNZAL,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZIND2,3*NNZAL,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZIND3,3*NNZAL,'DOUBLE',MASTER)

C     Do the work
      LRI = 1
      DO 20  J=1,MMCENT
         IWHO = -1
         IF (ZEROAL(J) .EQ. -1) GOTO 20
         CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
         CALL MPIXSEND(J, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
         LRI = LRI + 3
 20   CONTINUE

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(LRI, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data
      CALL DZERO(WRK(KTRES),N2ORBX)
      CALL MPI_REDUCE(WRK(KTRES),TRES,N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMCRO_M2')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmcro_s2 */
      SUBROUTINE QMMMCRO_S2(WRK,LWRK,IPRTMP)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include "implicit.h"
#include "maxorb.h"
#include "inflin.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "maxash.h"
#include "infind.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "ccinftap.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
#include "infdim.h"

      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMCRO_S2')

      QMMM = .TRUE.

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
      ENDIF
C    
      KUCMO   = 1
      KTRES   = KUCMO   + NORBT*NBAST
      KZYM2   = KTRES   + N2ORBX
      KZYM3   = KZYM2   + NORBT*NORBT
      KIND1   = KZYM3   + NORBT*NORBT
      KIND2   = KIND1   + 3*NNZAL
      KIND3   = KIND2   + 3*NNZAL
      KLAST   = KIND3   + 3*NNZAL
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMCRO_S2',-KLAST,LWRK2)

      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)

C     CRO
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV3,1,'INTEGER',MASTER)
      CALL MPIXBCAST(WRK(KZYM2),NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KZYM3),NORBT*NORBT,'DOUBLE',MASTER)

      CALL MPIXBCAST(WRK(KIND1),3*NNZAL,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KIND2),3*NNZAL,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KIND3),3*NNZAL,'DOUBLE',MASTER)

      CALL DZERO(WRK(KTRES),N2ORBX)

 101  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN
         CALL MMCRO_ITER2(J,WRK(KIND1+LRI-1),WRK(KIND2+LRI-1),
     *                    WRK(KIND3+LRI-1),WRK(KTRES),WRK(KUCMO),
     *                    ISYMT,ISYMV1,ISYMV2,ISYMV3,WRK(KZYM2),
     *                    WRK(KZYM3),WRK(KLAST),LWRK2,IPRTMP)
         GOTO 101
      ENDIF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KTRES),MPI_IN_PLACE,N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      CALL QEXIT('QMMMCRO_S2')
      RETURN
      END
#endif




