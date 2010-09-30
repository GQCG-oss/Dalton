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
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>
#include <infdim.h>

      DIMENSION BOVECS(*), CMO(*), XINDX(*), UDV(NASHDI,NASHDI)
      DIMENSION UDVTR(N2ASHX), DVTR(*), EVECS(KZYVAR,*)
      DIMENSION WRK(*), DV(*), CREF(*)

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LOCDEB
      LOGICAL FNDLAB
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      LOCDEB = .FALSE.

      CALL QENTER('QMMMLNO')

      IF (NNZAL .EQ. 0) THEN ! no pol. and nothing to do here
        CALL QEXIT('QMMMLNO')
        RETURN
      ENDIF

      IF (LGSPOL) THEN
        CALL QEXIT('QMMMLNO')
        RETURN
      ENDIF

      IF  ( ((NASHT.EQ.0) .AND. (TRPLET))  ) THEN ! no pol. contr.
C                                                   can arrise in that case
        CALL QEXIT('QMMMLNO')
        RETURN
      ENDIF


      KUCMO   = 1
      KUBO    = KUCMO   + NORBT*NBAST
      KFSOLMO = KUBO    + NOSIM*N2ORBX
      KUFSOL  = KFSOLMO + NNORBX
      KRXYO   = KUFSOL  + N2ORBX
      IF (TRPLET) THEN
        KRXYOT  = KRXYO   + NOSIM*N2ORBX
        KURXAC  = KRXYOT  + NOSIM*N2ORBX
      ELSE
        KURXAC  = KRXYO   + NOSIM*N2ORBX
      ENDIF
      KTRMO   = KURXAC  + N2ASHX
      KUTR    = KTRMO   + NNORBX
      KUTRX   = KUTR    + N2ORBX
      KEFIELD = KUTRX   + N2ORBX
      KINDMOM = KEFIELD + 3*NNZAL
      IF (MMMAT) THEN
        KINVAMT = KINDMOM + 3*NNZAL
        KWRK    = KINVAMT + (3*NNZAL)*(3*NNZAL)
      ELSE IF (MMITER) THEN
        KWRK    = KINDMOM + 3*NNZAL
      ENDIF
      LWRK1   = LWRK    - KWRK      


      IF (LWRK1 .LT. 0) CALL ERRWRK('QMMMLNO',-KWRK,LWRK)

      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KUBO),NOSIM*N2ORBX)
      CALL DZERO(WRK(KFSOLMO),NNORBX)
      CALL DZERO(WRK(KUFSOL),N2ORBX)
      CALL DZERO(WRK(KRXYO),NOSIM*N2ORBX)
      IF (TRPLET) CALL DZERO(WRK(KRXYOT),NOSIM*N2ORBX)
      CALL DZERO(WRK(KURXAC),N2ASHX)
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      CALL DZERO(WRK(KUTRX),N2ORBX)
      CALL DZERO(WRK(KEFIELD),3*NNZAL)
      CALL DZERO(WRK(KINDMOM),3*NNZAL)
      IF (MMMAT) CALL DZERO(WRK(KINVAMT),(3*NNZAL)*(3*NNZAL))

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
          CALL READT(LUQMMM,N*N,WRK(KINVAMT))
        ELSE
          CALL QUIT('Problem reading the matrix from the QMMMIM file.')
        ENDIF

        CALL GPCLOSE(LUQMMM,'KEEP')

        IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
          WRITE(LUPRI,*) ' The Relay matrix is read from file'
          CALL OUTPUT(WRK(KINVAMT),1,N,1,N,N,N,1,LUPRI)
        ENDIF
      ENDIF

      CALL UPKCMO(CMO,WRK(KUCMO))

      IF (NOSIM.GT.0) THEN
         CALL RSPZYM(NOSIM,BOVECS,WRK(KUBO))
         CALL DSCAL(NOSIM*N2ORBX,-1.0D0,WRK(KUBO),1)
      END IF

C     In contrast to the old code, the zero-order contribution, i.e.
C     - sum_a mu_ind_a t_a_pq has been included since we in 
C     sirirus do not subratct this from the fock/ks-operator. We here only 
C     ned to calculate the explicit electronic contribution, i.e.
C     - sum_a mu_ind_a[delta rho] t_q_pa, where delta rho is the first order
C     response/change in the density due to the applied perturbation. 
C

C     Backup diporg. We use diporg to transfer coordinates to int.
C     program.

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      KMAT = KWRK
      KLAST = KMAT + 3*NNBASX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMLNO',-KLAST,LWRK2)

      DO 200 IOSIM = 1,NOSIM

         JRXYO = KRXYO + (IOSIM-1)*N2ORBX
         IF (TRPLET) JRXYOT = KRXYOT + (IOSIM-1)*N2ORBX
         JUBO  = KUBO  + (IOSIM-1)*N2ORBX

C-------------------------------------------------------
C 1. Construct vector containing expectation value of electric 
C    field due to first order change in the density. Length of vector 
C    is NNZAL = no. of pol. sites.
C-------------------------------------------------------

         CALL DZERO(WRK(KEFIELD),3*NNZAL)

         IF (NODTOT .GE. 1) THEN
            CALL QMMMLNO_M1(IOSIM,UDV,UDVTR,WRK(KEFIELD),WRK(KUCMO),
     *                      WRK(JUBO),WRK(KLAST),LWRK2)
         ELSE
           LRI = 1             ! group-index in vector
           DO 201 J=1,MMCENT

             IF (ZEROAL(J) .EQ. -1) GOTO 201

             DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *               (MMCORD(2,J)-QMCOM(2))**2 +
     *               (MMCORD(3,J)-QMCOM(3))**2
             DIST = SQRT(DIST2)

             IF (DIST .GT. RCUTMM) THEN
               LRI = LRI + 3
               IF (LOCDEB) THEN
                 WRITE(LUPRI,*) 'Skipping site in F vector rsp 1', J
               ENDIF
               GOTO 201
             ENDIF

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
     &                   LWRK2,LABINT,INTREP,INTADR,J,TOFILE,
     &                   KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
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
     &                     (DIPORG(2)-CORD(2,M))**2 +
     &                     (DIPORG(3)-CORD(3,M))**2
                   IF (DISTC .LE. DIST) THEN
                     DIST = DISTC
                     MHIT = M
                   ENDIF
  123            CONTINUE
               ELSE IF (IDAMP .EQ. 2) THEN
                 DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                  (DIPORG(2)-QMCOM(2))**2 +
     &                  (DIPORG(3)-QMCOM(3))**2
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

C            x-component of (induced) electric field
             CALL DZERO(WRK(KTRMO),NNORBX)
             CALL DZERO(WRK(KUTR),N2ORBX)
             CALL DZERO(WRK(KUTRX),N2ORBX)
             CALL DZERO(WRK(KURXAC),N2ASHX)

             CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &                 NBAST,NORBT)
             CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
             CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))

             IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),WRK(KURXAC))
             IF (TRPLET) THEN
               EEX = SLVTLM(UDVTR,WRK(KURXAC),WRK(KUTRX),TAC)
             ELSE
               EEX = SLVQLM(UDV,WRK(KURXAC),WRK(KUTRX),TAC)
             ENDIF
             WRK(KEFIELD+LRI-1+0) = EEX

C            y-component of (induced) electric field
             CALL DZERO(WRK(KTRMO),NNORBX)
             CALL DZERO(WRK(KUTR),N2ORBX)
             CALL DZERO(WRK(KUTRX),N2ORBX)
             CALL DZERO(WRK(KURXAC),N2ASHX)

             CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                 WRK(KLAST),NBAST,NORBT)
             CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
             CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
             IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),WRK(KURXAC))
             IF (TRPLET) THEN
                EEY = SLVTLM(UDVTR,WRK(KURXAC),WRK(KUTRX),TAC)
             ELSE
                EEY = SLVQLM(UDV,WRK(KURXAC),WRK(KUTRX),TAC)
             ENDIF

             WRK(KEFIELD+LRI-1+1) = EEY

C            z-component of (induced) electric field
             CALL DZERO(WRK(KTRMO),NNORBX)
             CALL DZERO(WRK(KUTR),N2ORBX)
             CALL DZERO(WRK(KUTRX),N2ORBX)
             CALL DZERO(WRK(KURXAC),N2ASHX)

             CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                 WRK(KLAST),NBAST,NORBT)
             CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
             CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
             IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),WRK(KURXAC))
             IF (TRPLET) THEN
               EEZ = SLVTLM(UDVTR,WRK(KURXAC),WRK(KUTRX),TAC)
             ELSE
               EEZ = SLVQLM(UDV,WRK(KURXAC),WRK(KUTRX),TAC)
             ENDIF

             WRK(KEFIELD+LRI-1+2) = EEZ

             LRI = LRI + 3

  201      CONTINUE        
         ENDIF

         NDIM = 3*NNZAL
         NTOTI = MAX(NDIM,1)

         IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
           WRITE(LUPRI,*) 'Transformed electric fields'
           CALL OUTPUT(WRK(KEFIELD),1,NDIM,1,1,NDIM,1,1,LUPRI)
           WRITE(LUPRI,*)
         ENDIF

C        Calculate the induced dipoles corresponding to this field

         IF (MMMAT) THEN
           CALL DGEMV('N',NDIM,NDIM,ONE,WRK(KINVAMT),
     *                NTOTI,WRK(KEFIELD),1,ZERO,WRK(KINDMOM),1)
         ELSE IF (MMITER) THEN
           IOPT = 2 ! Do not read from file any previuos induced moments
           CALL F2QMMM(WRK(KEFIELD),NNZAL,WRK(KINDMOM),WRK(KLAST),LWRK2,
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

         IF (NODTOT .GE. 1) THEN
            IF (TRPLET) THEN
               CALL QMMMLNO_M2(WRK(JRXYOT),WRK(KUCMO),WRK(KINDMOM),
     *                      WRK(KWRK),LWRK1)
            ELSE
               CALL QMMMLNO_M2(WRK(JRXYO),WRK(KUCMO),WRK(KINDMOM),
     *                      WRK(KWRK),LWRK1)
            ENDIF
         ELSE
           LRI = 1

           DO 202 J=1,MMCENT

             IF (ZEROAL(J) .EQ. -1) GOTO 202

             DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *              (MMCORD(2,J)-QMCOM(2))**2 +
     *              (MMCORD(3,J)-QMCOM(3))**2
             DIST = SQRT(DIST2)

             IF (DIST .GT. RCUTMM) THEN
               LRI = LRI + 3
               IF (LOCDEB) THEN
                 WRITE(LUPRI,*) 'Skipping site in F vector rsp 2', J
               ENDIF
               GOTO 202
             ENDIF

             CALL DZERO(WRK(KTRMO),NNORBX)
             CALL DZERO(WRK(KUTR),N2ORBX)
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
     &                   LWRK2,LABINT,INTREP,INTADR,J,TOFILE,
     &                   KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
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
     &                     (DIPORG(2)-CORD(2,M))**2 +
     &                     (DIPORG(3)-CORD(3,M))**2
                   IF (DISTC .LE. DIST) THEN
                     DIST = DISTC
                     MHIT = M
                   ENDIF
  124            CONTINUE
               ELSE IF (IDAMP .EQ. 2) THEN
                 DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                  (DIPORG(2)-QMCOM(2))**2 +
     &                  (DIPORG(3)-QMCOM(3))**2
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

C            x-component of electric field 
             CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &                 NBAST,NORBT)
             CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
             FACx = -WRK(KINDMOM+LRI-1+0)   
             IF (TRPLET) THEN 
               CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,WRK(JRXYOT),1)
             ELSE
               CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,WRK(JRXYO),1)
             ENDIF

C            y-component of electric field
             CALL DZERO(WRK(KTRMO),NNORBX)
             CALL DZERO(WRK(KUTR),N2ORBX)

             CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                 WRK(KLAST),NBAST,NORBT)
             CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
             FACy = -WRK(KINDMOM+LRI-1+1)
             IF (TRPLET) THEN
               CALL DAXPY(N2ORBX,FACy,WRK(KUTR),1,WRK(JRXYOT),1)
             ELSE
               CALL DAXPY(N2ORBX,FACy,WRK(KUTR),1,WRK(JRXYO),1)
             ENDIF

C            z-component of (induced) electric field
             CALL DZERO(WRK(KTRMO),NNORBX)
             CALL DZERO(WRK(KUTR),N2ORBX)

             CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                 WRK(KLAST),NBAST,NORBT)
             CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
             FACz = -WRK(KINDMOM+LRI-1+2)
             IF (TRPLET) THEN
               CALL DAXPY(N2ORBX,FACz,WRK(KUTR),1,WRK(JRXYOT),1)
             ELSE
               CALL DAXPY(N2ORBX,FACz,WRK(KUTR),1,WRK(JRXYO),1)
             ENDIF

             LRI = LRI + 3

  202      CONTINUE 
         ENDIF
  200 CONTINUE

C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

C--------------------------------------------
C 3.  Add effective operators to response
C--------------------------------------------


      IF (TRPLET) THEN
         CALL SLVSOR(.TRUE.,.FALSE.,NOSIM,UDVTR,EVECS(1,1),WRK(KRXYO))
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYOT))
      ELSE
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYO))
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
#include <implicit.h>
#include <dummy.h>
#include <maxorb.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <infvar.h>
#include <infrsp.h>
#include <infpri.h>
#include <rspprp.h>
#include <infcr.h>
#include <inftap.h>
#include <qrinf.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <ccinftap.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>

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

      IF (IPQMMM .GT. 1) CALL TIMER('START ',TIMSTR,TIMEND)

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
        KINVAMT = KUTR    + N2ORBX
        KEF1    = KINVAMT + (3*NNZAL)*(3*NNZAL)
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
      IF (MMMAT) CALL DZERO(WRK(KINVAMT),3*NNZAL*3*NNZAL)
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

C       Read the relay matrix from file

        LUQMMM = -1
        IF (LUQMMM .LT. 0) THEN
          CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
        ENDIF
        REWIND(LUQMMM)

        IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
          CALL READT(LUQMMM,N*N,WRK(KINVAMT))
        ELSE
          CALL QUIT('Problem reading the matrix from the QMMMIM file.')
        ENDIF

        CALL GPCLOSE(LUQMMM,'KEEP')

        IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
          WRITE(LUPRI,*) ' The Relay matrix is read from file'
          CALL OUTPUT(WRK(KINVAMT),1,N,1,N,N,N,1,LUPRI)
        ENDIF

      ENDIF ! MMMAT

C     Backup the dipole origin

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      KMAT = KWRK
      KLAST = KMAT + 3*NNBASX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMQRO',-KLAST,LWRK2)

      IF (NODTOT .GE. 1) THEN
         CALL QMMMQRO_M1(WRK(KEF1),WRK(KEF2),WRK(KUCMO),
     &                   ISYMT,ISYMV1,ISYMV2,ZYM1,ZYM2,WRK(KLAST),LWRK2)
      ELSE
        LRI = 1                ! group-index in vector
        DO 100 J = 1,MMCENT 

          IF (ZEROAL(J) .EQ. -1) GOTO 100

          DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *            (MMCORD(2,J)-QMCOM(2))**2 +
     *            (MMCORD(3,J)-QMCOM(3))**2
          DIST = SQRT(DIST2)

          IF (DIST .GT. RCUTMM) THEN
            LRI = LRI + 3
            IF (LOCDEB) THEN
              WRITE(LUPRI,*) 'Skipping site in F vector rsp 3', J
            ENDIF
            GOTO 100
          ENDIF

C         x-component
          F1=ZERO
          F2=ZERO
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)
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
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
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
     &                  (DIPORG(2)-CORD(2,M))**2 +
     &                  (DIPORG(3)-CORD(3,M))**2
                IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
                ENDIF
  125         CONTINUE
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

          CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
C
          IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMQRO')
            WRK(KEF1+LRI-1+0) = F1
          ELSE
            WRK(KEF1+LRI-1+0) = ZERO
          ENDIF 

          IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMQRO')
            WRK(KEF2+LRI-1+0) = F2
          ELSE
            WRK(KEF2+LRI-1+0) = ZERO
          ENDIF

C         y-component
          F1=ZERO
          F2=ZERO
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)

          CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &              WRK(KLAST),NBAST,NORBT)
          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

          IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QM3MMRO')
            WRK(KEF1+LRI-1+1) = F1
          ELSE
            WRK(KEF1+LRI-1+1) = ZERO
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMQRO')
            WRK(KEF2+LRI-1+1) = F2
          ELSE
            WRK(KEF2+LRI-1+1) = ZERO
          ENDIF

C         z-component
          F1=ZERO
          F2=ZERO
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)

          CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &              WRK(KLAST),NBAST,NORBT)
          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

          IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMQRO')
            WRK(KEF1+LRI-1+2) = F1
          ELSE
            WRK(KEF1+LRI-1+2) = ZERO
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMQRO')
            WRK(KEF2+LRI-1+2) = F2
          ELSE
            WRK(KEF2+LRI-1+2) = ZERO
          ENDIF

          LRI = LRI + 3

  100   CONTINUE
      ENDIF

C     Calculate the induced dipoles corresponding to the 
C     (transformed) fields 1 and 2

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
        CALL DGEMV('N',NDIM,NDIM,ONE,WRK(KINVAMT),
     *             NTOTI,WRK(KEF1),1,ZERO,WRK(KIND1),1)
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
         CALL DGEMV('N',NDIM,NDIM,ONE,WRK(KINVAMT),
     *             NTOTI,WRK(KEF2),1,ZERO,WRK(KIND2),1)
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

C     Now we need the integrals again transformed to the MO basis 
C     and also one-index transformed integrals. 
C     These could have been stored
C     before, but we choose to do all of this integral-direct.

      IF (NODTOT .GE. 1) THEN
         CALL QMMMQRO_M2(WRK(KIND1),WRK(KIND2),WRK(KUCMO),
     &        WRK(KTRES),ISYMT,ISYMV2,ZYM2,WRK(KLAST),LWRK2)
      ELSE
        LRI = 1                ! group-index in vector
        DO 101 J = 1,MMCENT 

          IF (ZEROAL(J) .EQ. -1) GOTO 101

          DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *            (MMCORD(2,J)-QMCOM(2))**2 +
     *            (MMCORD(3,J)-QMCOM(3))**2
          DIST = SQRT(DIST2)

          IF (DIST .GT. RCUTMM) THEN
            LRI = LRI + 3
            IF (LOCDEB) THEN
              WRITE(LUPRI,*) 'Skipping site in F vector rsp 4', J
            ENDIF
            GOTO 101
          ENDIF

C         x-component
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)
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
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
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
  126         CONTINUE
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

          CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
C
          F1 = -WRK(KIND1+LRI-1+0) 
          CALL DZERO(WRK(KTLMA),N2ORBX)
          CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
          CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)

          F2 = -0.50D0*WRK(KIND2+LRI-1+0)
          CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)

C         y-component
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)

          CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &              WRK(KLAST),NBAST,NORBT)
          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

          F1 = -WRK(KIND1+LRI-1+1)
          CALL DZERO(WRK(KTLMA),N2ORBX)
          CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
          CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)

          F2 = -0.50D0*WRK(KIND2+LRI-1+1)
          CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)

C         z-component
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)

          CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &              WRK(KLAST),NBAST,NORBT)
          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

          F1 = -WRK(KIND1+LRI-1+2)
          CALL DZERO(WRK(KTLMA),N2ORBX)
          CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
          CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)   

          F2 = -0.50D0*WRK(KIND2+LRI-1+2)
          CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)

          LRI = LRI + 3

  101   CONTINUE
      ENDIF
C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

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
      IF (IPQMMM .GT. 1) CALL TIMER('QMMMQRO',TIMSTR,TIMEND)
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
#include <implicit.h>
#include <maxorb.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <infvar.h>
#include <infrsp.h>
#include <infpri.h>
#include <rspprp.h>
#include <infcr.h>
#include <inftap.h>
#include <qrinf.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <ccinftap.h>
#include <nuclei.h>
#include <dummy.h>
#include <infpar.h>

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

      IF (IPQMMM .GT. 1) CALL TIMER('START ',TIMSTR,TIMEND)

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
        KINVAMT = KUTR    + N2ORBX
        KEF1    = KINVAMT + (3*NNZAL)*(3*NNZAL)
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
      IF (MMMAT) CALL DZERO(WRK(KINVAMT),3*NNZAL*3*NNZAL)
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

C       Read the relay matrix from file        

        LUQMMM = -1
        IF (LUQMMM .LT. 0) THEN
          CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
        ENDIF
        REWIND(LUQMMM)

        IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
          CALL READT(LUQMMM,N*N,WRK(KINVAMT))
        ELSE
          CALL QUIT('Problem reading the matrix from the QMMMIM file.')
        ENDIF

        CALL GPCLOSE(LUQMMM,'KEEP')

        IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) .OR. (IPQMMM .GE. 5) )THEN
          WRITE(LUPRI,*) ' The Relay matrix is read from file'
          CALL OUTPUT(WRK(KINVAMT),1,N,1,N,N,N,1,LUPRI)
        ENDIF

      ENDIF 

C     Backup the dipole origin

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      KMAT = KWRK
      KLAST = KMAT + 3*NNBASX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMCRO',-KLAST,LWRK2)

      IF (NODTOT .GE. 1) THEN
         CALL QMMMCRO_M1(WRK(KEF1),WRK(KEF2),WRK(KEF3),WRK(KUCMO),
     &                   ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM1,ZYM2,ZYM3,
     &                   WRK(KLAST),LWRK2)
      ELSE
        LRI = 1                ! group-index in vector
        DO 100 J = 1,MMCENT 

          IF (ZEROAL(J) .EQ. -1) GOTO 100

          DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *            (MMCORD(2,J)-QMCOM(2))**2 +
     *            (MMCORD(3,J)-QMCOM(3))**2
          DIST = SQRT(DIST2)

          IF (DIST .GT. RCUTMM) THEN
            LRI = LRI + 3
            IF (LOCDEB) THEN
              WRITE(LUPRI,*) 'Skipping site in F vector rsp 5', J
            ENDIF
            GOTO 100
          ENDIF


          CALL DZERO(WRK(KMAT),3*NNBASX)

          KPATOM = 0
          NOCOMP = 3
          TOFILE = .FALSE.
          TRIMAT = .TRUE.
          EXP1VL = .FALSE.
          DIPORG(1) = MMCORD(1,J)
          DIPORG(2) = MMCORD(2,J)
          DIPORG(3) = MMCORD(3,J)

C         Get 1-electron integrals

          RUNQM3 = .TRUE.
          CALL GET1IN(WRK(KMAT),'NEFIELD',NOCOMP,WRK(KLAST),
     &                LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
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
     &                  (DIPORG(2)-CORD(2,M))**2 +
     &                  (DIPORG(3)-CORD(3,M))**2
                IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
                ENDIF
  123         CONTINUE
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

C         x-component
          F1=ZERO
          F2=ZERO
          F3=ZERO
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)

C 
          CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)

C         Unpack

          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

          IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)

C           1-index transformation

            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)

C           Evaluate expectation value

            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMCRO')
            WRK(KEF1+LRI-1+0) = F1
          ELSE
            WRK(KEF1+LRI-1+0) = ZERO
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMCRO')
            WRK(KEF2+LRI-1+0) = F2
          ELSE
            WRK(KEF2+LRI-1+0) = ZERO
          END IF

          IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV1))
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTLMB),WRK(KTLMA),ISYMV3)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F3,200,'QMMMCRO')
            WRK(KEF3+LRI-1+0) = F3
          ELSE
            WRK(KEF3+LRI-1+0) = ZERO
          END IF

C         y-component
          F1=ZERO
          F2=ZERO
          F3=ZERO
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)

C         AO to MO transformation

          CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)

C         Unpack from tri to full

          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

          IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMCRO')
            WRK(KEF1+LRI-1+1) = F1
          ELSE
             WRK(KEF1+LRI-1+1) = ZERO
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMCRO')
            WRK(KEF2+LRI-1+1) = F2
          ELSE
            WRK(KEF2+LRI-1+1) = ZERO
          END IF

          IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV1))
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTLMB),WRK(KTLMA),ISYMV3)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F3,200,'QMMMCRO')
            WRK(KEF3+LRI-1+1) = F3
          ELSE
            WRK(KEF3+LRI-1+1) = ZERO
          END IF

C         z-component
          F1=ZERO
          F2=ZERO
          F3=ZERO
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)

          CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

          IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMCRO')
            WRK(KEF1+LRI-1+2) = F1
          ELSE
            WRK(KEF1+LRI-1+2) = ZERO
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMCRO')
            WRK(KEF2+LRI-1+2) = F2
          ELSE
            WRK(KEF2+LRI-1+2) = ZERO
          END IF

          IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV1))
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTLMB),WRK(KTLMA),ISYMV3)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F3,200,'QMMMCRO')
            WRK(KEF3+LRI-1+2) = F3
          ELSE
            WRK(KEF3+LRI-1+2) = ZERO
          END IF

          LRI = LRI + 3

  100   CONTINUE
      ENDIF

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
        CALL DGEMV('N',NDIM,NDIM,ONE,WRK(KINVAMT),
     *             NTOTI,WRK(KEF1),1,ZERO,WRK(KIND1),1)
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
        CALL DGEMV('N',NDIM,NDIM,ONE,WRK(KINVAMT),
     *             NTOTI,WRK(KEF2),1,ZERO,WRK(KIND2),1)
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
        CALL DGEMV('N',NDIM,NDIM,ONE,WRK(KINVAMT),
     *            NTOTI,WRK(KEF3),1,ZERO,WRK(KIND3),1)
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

      IF (NODTOT .GE. 1) THEN
         CALL QMMMCRO_M2(WRK(KIND1),WRK(KIND2),WRK(KIND3),WRK(KUCMO),
     &        WRK(KTRES),ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM2,ZYM3,
     &        WRK(KLAST),LWRK2)
      ELSE
        LRI = 1                ! group-index in vector
        DO 101 J = 1,MMCENT 

          IF (ZEROAL(J) .EQ. -1) GOTO 101

          DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *            (MMCORD(2,J)-QMCOM(2))**2 +
     *            (MMCORD(3,J)-QMCOM(3))**2
          DIST = SQRT(DIST2)

          IF (DIST .GT. RCUTMM) THEN
            LRI = LRI + 3
            IF (LOCDEB) THEN
              WRITE(LUPRI,*) 'Skipping site in F vector rsp 6', J
            ENDIF
            GOTO 101
          ENDIF

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
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
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
     &                  (DIPORG(2)-CORD(2,M))**2 +
     &                  (DIPORG(3)-CORD(3,M))**2
                IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
                ENDIF
  124         CONTINUE
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

C         x-component
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)
          
          CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

          IF (ISYMT.EQ.ISYMV1) THEN 
            F1 = -(1.0D0/2.0D0)*WRK(KIND1+LRI-1+0)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV2))
            CALL DAXPY(N2ORBX,F1,WRK(KTLMB),1,WRK(KTRES),1)
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            F2 = -(1.0D0/2.0D0)*WRK(KIND2+LRI-1+0)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,ZYM3,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL DAXPY(N2ORBX,F2,WRK(KTLMA),1,WRK(KTRES),1)
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            F3 = -(1.0D0/6.0D0)*WRK(KIND3+LRI-1+0)
            CALL DAXPY(N2ORBX,F3,WRK(KUTR),1,WRK(KTRES),1)
          ENDIF

C         y-component
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)

          CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

          IF (ISYMT.EQ.ISYMV1) THEN
            F1 = -(1.0D0/2.0D0)*WRK(KIND1+LRI-1+1)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV2))
            CALL DAXPY(N2ORBX,F1,WRK(KTLMB),1,WRK(KTRES),1)
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            F2 = -(1.0D0/2.0D0)*WRK(KIND2+LRI-1+1)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,ZYM3,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL DAXPY(N2ORBX,F2,WRK(KTLMA),1,WRK(KTRES),1)
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            F3 = -(1.0D0/6.0D0)*WRK(KIND3+LRI-1+1)
            CALL DAXPY(N2ORBX,F3,WRK(KUTR),1,WRK(KTRES),1)
          ENDIF

C         z-component
          CALL DZERO(WRK(KTRMO),NNORBX)
          CALL DZERO(WRK(KUTR),N2ORBX)

          CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
          CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

          IF (ISYMT.EQ.ISYMV1) THEN
            F1 = -(1.0D0/2.0D0)*WRK(KIND1+LRI-1+2)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV2))
            CALL DAXPY(N2ORBX,F1,WRK(KTLMB),1,WRK(KTRES),1)
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            F2 = -(1.0D0/2.0D0)*WRK(KIND2+LRI-1+2)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,ZYM3,WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL DAXPY(N2ORBX,F2,WRK(KTLMA),1,WRK(KTRES),1)
          ENDIF

          IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            F3 = -(1.0D0/6.0D0)*WRK(KIND3+LRI-1+2)
            CALL DAXPY(N2ORBX,F3,WRK(KUTR),1,WRK(KTRES),1)
          ENDIF

          LRI = LRI + 3

  101   CONTINUE
      ENDIF
C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

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

      IF (IPQMMM .GT. 1) CALL TIMER('QMMMCRO',TIMSTR,TIMEND)

      CALL QEXIT('QMMMCRO')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmlno_m1 */
      SUBROUTINE QMMMLNO_M1(IOSIM,UDV,UDVTR,CEFIELD,UCMO,
     *                      UBO,WRK,LWRK)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>
#include <infdim.h>

      INTEGER IOSIM
      DIMENSION WRK(LWRK), UBO(*), UCMO(*)
      DIMENSION UDVTR(*), UDV(NASHDI,NASHDI), CEFIELD(3*NNZAL)

      CALL QENTER('QMMMLNO_M1')

      KEFIELD = 1
      KLAST   = KEFIELD + 3*NNZAL
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMLNO_M1',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = 28
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
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>
#include <gnrinf.h>
#include <infdim.h>

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMLNO_S1')

      QM3 = .TRUE.

C      print*,mmcent,natoms,nnbasx,nnzal,qmdamp,mmcord(1,1),
C     &       qmcom(1),rcutmm,nbast,isytp,nctot,cord(1,1),
C     &       charge(1),nucind,nucdep,luprop,npatom,ipatom(1),
C     &       trplet,nnorbx,n2orbx,n2ashx,nashdi,norbt,nisht,
C     &       norb,isx(1)
C      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)

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
C      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
C      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
C      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)
C
C      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
C      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NISHT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NORB,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,NORB,'INTEGER',MASTER)

C      print*,mmcent,natoms,nnbasx,nnzal,qmdamp,mmcord(1,1),
C     &       qmcom(1),rcutmm,nbast,isytp,nctot,cord(1,1),
C     &       charge(1),nucind,nucdep,luprop,npatom,ipatom(1),
C      CALL QUIT('DONE')
      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)
C    
      KUCMO  = 1
      JUBO   = KUCMO + NORBT*NBAST
      KURXAC = JUBO  + N2ORBX
      KTRMO  = KURXAC+ N2ASHX
      KUTR   = KTRMO + NNORBX
      KUTRX  = KUTR  + N2ORBX
      KEFIELD= KUTRX + N2ORBX

      IF (TRPLET) THEN
         KUDVTR  = KEFIELD + 3*NNZAL
         KWRK    = KUDVTR  + N2ASHX
         CALL MPIXBCAST(WRK(KUDVTR),N2ASHX,'DOUBLE',MASTER)
      ELSE
         NNASH   = NASHDI*NASHDI
         KUDV    = KEFIELD + 3*NNZAL 
         KWRK    = KUDV    + NASHDI*NASHDI
         CALL MPIXBCAST(WRK(KUDV),NASHDI*NASHDI,'DOUBLE',MASTER)
      ENDIF

      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(JUBO),N2ORBX,'DOUBLE',MASTER)

      CALL DZERO(WRK(KEFIELD),3*NNZAL)

      KMAT = KWRK
      KLAST = KMAT + 3*NNBASX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMLNO_S1',-KLAST,LWRK2)

 201  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN

         DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *           (MMCORD(2,J)-QMCOM(2))**2 +
     *           (MMCORD(3,J)-QMCOM(3))**2
         DIST = SQRT(DIST2)

         IF (DIST .GT. RCUTMM) THEN
            GOTO 201
         ENDIF

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
     &                    (DIPORG(2)-CORD(2,M))**2 +
     &                    (DIPORG(3)-CORD(3,M))**2
                  IF (DISTC .LE. DIST) THEN
                     DIST = DISTC
                     MHIT = M
                  ENDIF
 123           CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
               DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                (DIPORG(2)-QMCOM(2))**2 +
     &                (DIPORG(3)-QMCOM(3))**2
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

C        x-component of (induced) electric field
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)
         CALL DZERO(WRK(KUTRX),N2ORBX)
         CALL DZERO(WRK(KURXAC),N2ASHX)

         CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &             NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
         CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))

         IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),WRK(KURXAC))
         IF (TRPLET) THEN
            EEX = SLVTLM(WRK(KUDVTR),WRK(KURXAC),WRK(KUTRX),TAC)
         ELSE
            EEX = SLVQLM(WRK(KUDV),WRK(KURXAC),WRK(KUTRX),TAC)
         ENDIF
         WRK(KEFIELD+LRI-1+0) = EEX

C        y-component of (induced) electric field
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)
         CALL DZERO(WRK(KUTRX),N2ORBX)
         CALL DZERO(WRK(KURXAC),N2ASHX)

         CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &             WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
         CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
         IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),WRK(KURXAC))
         IF (TRPLET) THEN
            EEY = SLVTLM(WRK(KUDVTR),WRK(KURXAC),WRK(KUTRX),TAC)
         ELSE
            EEY = SLVQLM(WRK(KUDV),WRK(KURXAC),WRK(KUTRX),TAC)
         ENDIF
         WRK(KEFIELD+LRI-1+1) = EEY

C        z-component of (induced) electric field
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)
         CALL DZERO(WRK(KUTRX),N2ORBX)
         CALL DZERO(WRK(KURXAC),N2ASHX)

         CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &             WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
         CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
         IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),WRK(KURXAC))
         IF (TRPLET) THEN
            EEZ = SLVTLM(WRK(KUDVTR),WRK(KURXAC),WRK(KUTRX),TAC)
         ELSE
            EEZ = SLVQLM(WRK(KUDV),WRK(KURXAC),WRK(KUTRX),TAC)
         ENDIF

         WRK(KEFIELD+LRI-1+2) = EEZ

         GO TO 201
      ENDIF

C     No more integrals to calculate

      CALL MPI_REDUCE(WRK(KEFIELD),WRK(KEFIELD),3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

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
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>
#include <infdim.h>

      DIMENSION WRK(LWRK)
      DIMENSION UCMO(*), RXYO(N2ORBX), XINDMOM(3*NNZAL)

      CALL QENTER('QMMMLNO_M2')

      KRXYO = 1
      KLAST = KRXYO + N2ORBX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMLNO_M2',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = 29

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

C     Send info
C      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)

C     Damping
C      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
C      IF (QMDAMP) THEN
C         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
C         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
C         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
C         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
C      ENDIF
C
C      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
C      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
C      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)
C
C      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
C      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)

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
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>
#include <gnrinf.h>
#include <infdim.h>

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMLNO_S2')

      QM3 = .TRUE.
      
C      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)

C     Damping
C      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
C      IF (QMDAMP) THEN
C         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
C         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
C         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
C         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
C      ENDIF
C
C      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
C      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
C      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)
C
C      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
C      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)
C    
      KINDMOM = 1
      KUCMO   = KINDMOM + 3*NNZAL
      KTRMO   = KUCMO   + NORBT*NBAST
      KUTR    = KTRMO   + NNORBX
      KMAT    = KUTR    + N2ORBX
      JRXYO   = KMAT    + 3*NNBASX
      KLAST   = JRXYO   + N2ORBX
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMLNO_S2',-KLAST,LWRK2)

      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KINDMOM),3*NNZAL,'DOUBLE',MASTER)

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL DZERO(WRK(JRXYO),N2ORBX)

 201  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN

         DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *           (MMCORD(2,J)-QMCOM(2))**2 +
     *           (MMCORD(3,J)-QMCOM(3))**2
         DIST = SQRT(DIST2)

         IF (DIST .GT. RCUTMM) THEN
            GOTO 201
         ENDIF

         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)
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
     &                    (DIPORG(2)-CORD(2,M))**2 +
     &                    (DIPORG(3)-CORD(3,M))**2
                  IF (DISTC .LE. DIST) THEN
                     DIST = DISTC
                     MHIT = M
                  ENDIF
 124           CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
               DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                (DIPORG(2)-QMCOM(2))**2 +
     &                (DIPORG(3)-QMCOM(3))**2
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

C        x-component of electric field 
         CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &             NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
         FACx = -WRK(KINDMOM+LRI-1+0)   
         CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,WRK(JRXYO),1)

C        y-component of electric field
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &             WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
         FACy = -WRK(KINDMOM+LRI-1+1)
         CALL DAXPY(N2ORBX,FACy,WRK(KUTR),1,WRK(JRXYO),1)

C        z-component of (induced) electric field
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &             WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
         FACz = -WRK(KINDMOM+LRI-1+2)
         CALL DAXPY(N2ORBX,FACz,WRK(KUTR),1,WRK(JRXYO),1)

         GO TO 201
      ENDIF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(JRXYO),WRK(JRXYO),N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('QMMMLNO_S2')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmqro_m1 */
      SUBROUTINE QMMMQRO_M1(EF1,EF2,UCMO,
     &                      ISYMT,ISYMV1,ISYMV2,ZYM1,ZYM2,WRK,LWRK)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>

      DIMENSION WRK(LWRK)
      DIMENSION UCMO(*),ZYM1(*),ZYM2(*),EF1(*),EF2(*)

      CALL QENTER('QMMMQRO_M1')

      KEF1  = 1
      KEF2  = KEF1 + 3*NNZAL
      KLAST = KEF2 + 3*NNZAL
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMQRO_M1',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = 30
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

C     Send info

C     Integers for arrays/matrixes
C      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)

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

C      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
C      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
C      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)
C
C      CALL MPIXBCAST(IATNOW,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
C      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)

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
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>
#include <gnrinf.h>
#include <infdim.h>

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMQRO_S1')

      QM3 = .TRUE.

C      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)

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
C      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
C      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
C      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)
C
C      CALL MPIXBCAST(IATNOW,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
C      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
C
      KUCMO   = 1
      KTRMO   = KUCMO   + NORBT*NBAST
      KUTR    = KTRMO   + NNORBX
      KMAT    = KUTR    + N2ORBX
      KEF1    = KMAT    + 3*NNBASX
      KEF2    = KEF1    + 3*NNZAL
      KTLMA   = KEF2    + 3*NNZAL
      KTLMB   = KTLMA   + N2ORBX
      KZYM1   = KTLMB   + N2ORBX
      KZYM2   = KZYM1   + NORBT*NORBT
      KLAST   = KZYM2   + NORBT*NORBT
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMQRO_S1',-KLAST,LWRK2)

      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)

C     QRO
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(WRK(KZYM1),NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KZYM2),NORBT*NORBT,'DOUBLE',MASTER)

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL DZERO(WRK(KEF1),3*NNZAL)
      CALL DZERO(WRK(KEF2),3*NNZAL)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN

         DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *           (MMCORD(2,J)-QMCOM(2))**2 +
     *           (MMCORD(3,J)-QMCOM(3))**2
         DIST = SQRT(DIST2)

         IF (DIST .GT. RCUTMM) THEN
            GOTO 100
         ENDIF

C        x-component
         F1=ZERO
         F2=ZERO
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)
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
     &                LWRK2,LABINT,INTREP,INTADR,J,TOFILE,
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
     &                    (DIPORG(2)-CORD(2,M))**2 +
     &                    (DIPORG(3)-CORD(3,M))**2
                  IF (DISTC .LE. DIST) THEN
                     DIST = DISTC
                     MHIT = M
                  ENDIF
 125           CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
               DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                (DIPORG(2)-QMCOM(2))**2 +
     &                (DIPORG(3)-QMCOM(3))**2
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

         CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &             NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
C
         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMQRO')
            WRK(KEF1+LRI-1+0) = F1
         ELSE
            WRK(KEF1+LRI-1+0) = ZERO
         ENDIF 

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMQRO')
            WRK(KEF2+LRI-1+0) = F2
         ELSE
            WRK(KEF2+LRI-1+0) = ZERO
         ENDIF

C        y-component
         F1=ZERO
         F2=ZERO
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &             WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QM3MMRO')
            WRK(KEF1+LRI-1+1) = F1
         ELSE
            WRK(KEF1+LRI-1+1) = ZERO
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMQRO')
            WRK(KEF2+LRI-1+1) = F2
         ELSE
            WRK(KEF2+LRI-1+1) = ZERO
         ENDIF

C        z-component
         F1=ZERO
         F2=ZERO
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &             WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMQRO')
            WRK(KEF1+LRI-1+2) = F1
         ELSE
            WRK(KEF1+LRI-1+2) = ZERO
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMQRO')
            WRK(KEF2+LRI-1+2) = F2
         ELSE
            WRK(KEF2+LRI-1+2) = ZERO
         ENDIF
         GO TO 100
      ENDIF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KEF1),WRK(KEF1),3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(WRK(KEF2),WRK(KEF2),3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

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
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>

      DIMENSION WRK(LWRK)
      DIMENSION UCMO(*),ZYM2(*),TRES(*),ZIND1(3*NNZAL),ZIND2(3*NNZAL)

      CALL QENTER('QMMMQRO_M2')

      KTRES = 1
      KLAST = KTRES + N2ORBX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMQRO_M2',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = 31
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

C     Send info
C     Integers for arrays/matrixes
C      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)

C     Damping
C      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
C      IF (QMDAMP) THEN
C         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
C         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
C         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
C         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
C      ENDIF
C
C      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
C      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
C      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)
C
C      CALL MPIXBCAST(IATNOW,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
C      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)

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
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>
#include <gnrinf.h>
#include <infdim.h>

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMQRO_S2')

      QM3 = .TRUE.

C      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)

C     Damping
C      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
C      IF (QMDAMP) THEN
C         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
C         CALL MPIXBCAST(QMPOL,NATOMS,'DOUBLE',MASTER)
C         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
C         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
C         CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
C      ENDIF
C
C      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
C      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
C      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)
C
C      CALL MPIXBCAST(IATNOW,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
C      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
C
C      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
C      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
C      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
C    
      KUCMO   = 1
      KTRMO   = KUCMO   + NORBT*NBAST
      KUTR    = KTRMO   + NNORBX
      KMAT    = KUTR    + N2ORBX
      KTRES   = KMAT    + 3*NNBASX
      KTLMA   = KTRES   + N2ORBX
      KZYM2   = KTLMA   + N2ORBX
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

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL DZERO(WRK(KTRES),N2ORBX)


 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN

         DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *           (MMCORD(2,J)-QMCOM(2))**2 +
     *           (MMCORD(3,J)-QMCOM(3))**2
         DIST = SQRT(DIST2)

         IF (DIST .GT. RCUTMM) THEN
            GOTO 100
         ENDIF

C        x-component
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)
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
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
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
     &                 (DIPORG(2)-CORD(2,M))**2 +
     &                 (DIPORG(3)-CORD(3,M))**2
                  IF (DISTC .LE. DIST) THEN
                     DIST = DISTC
                     MHIT = M
                  ENDIF
 126           CONTINUE
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

         CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
C
         F1 = -WRK(KIND1+LRI-1+0) 
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)

         F2 = -0.50D0*WRK(KIND2+LRI-1+0)
         CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)

C        y-component
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &              WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         F1 = -WRK(KIND1+LRI-1+1)
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)

         F2 = -0.50D0*WRK(KIND2+LRI-1+1)
         CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)

C        z-component
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &              WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         F1 = -WRK(KIND1+LRI-1+2)
         CALL DZERO(WRK(KTLMA),N2ORBX)
         CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KUTR),WRK(KTLMA),ISYMT)
         CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)   

         F2 = -0.50D0*WRK(KIND2+LRI-1+2)
         CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)

         GO TO 100
      ENDIF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KTRES),WRK(KTRES),N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('QMMMQRO_S2')
      RETURN
      END
C-------------------------------------------------------------------------------
C  /* Deck qmmmcro_m1 */
      SUBROUTINE QMMMCRO_M1(EF1,EF2,EF3,UCMO,
     &                      ISYMT,ISYMV1,ISYMV2,ISYMV3,ZYM1,ZYM2,ZYM3,
     &                      WRK,LWRK)
C
C  JK, Nov.08
C  Parallel version. AHS Nov.09
C  
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>

      DIMENSION WRK(LWRK)
      DIMENSION UCMO(*),ZYM1(*),ZYM2(*),ZYM3(*),EF1(*),EF2(*),EF3(*)

      CALL QENTER('QMMMCRO_M1')

      KEF1  = 1
      KEF2  = KEF1 + 3*NNZAL
      KEF3  = KEF2 + 3*NNZAL
      KLAST = KEF3 + 3*NNZAL
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMCRO_M1',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = 32
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

C     Send info

C     Integers for arrays/matrixes
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)

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

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)

      CALL MPIXBCAST(IATNOW,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)

      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)

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
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>
#include <gnrinf.h>
#include <infdim.h>

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMCRO_S1')

      QM3 = .TRUE.

      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)

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
      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)

      CALL MPIXBCAST(IATNOW,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)

      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
C
      KUCMO   = 1
      KTRMO   = KUCMO   + NORBT*NBAST
      KUTR    = KTRMO   + NNORBX
      KMAT    = KUTR    + N2ORBX
      KEF1    = KMAT    + 3*NNBASX
      KEF2    = KEF1    + 3*NNZAL
      KEF3    = KEF2    + 3*NNZAL
      KTLMA   = KEF3    + 3*NNZAL
      KTLMB   = KTLMA   + N2ORBX
      KZYM1   = KTLMB   + N2ORBX
      KZYM2   = KZYM1   + NORBT*NORBT
      KZYM3   = KZYM2   + NORBT*NORBT
      KLAST   = KZYM3   + NORBT*NORBT
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMCRO_S1',-KLAST,LWRK2)

      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)

C     CRO
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV3,1,'INTEGER',MASTER)
      CALL MPIXBCAST(WRK(KZYM1),NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KZYM2),NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KZYM3),NORBT*NORBT,'DOUBLE',MASTER)

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL DZERO(WRK(KEF1),3*NNZAL)
      CALL DZERO(WRK(KEF2),3*NNZAL)
      CALL DZERO(WRK(KEF3),3*NNZAL)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN

         DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *           (MMCORD(2,J)-QMCOM(2))**2 +
     *           (MMCORD(3,J)-QMCOM(3))**2
         DIST = SQRT(DIST2)

         IF (DIST .GT. RCUTMM) THEN
            GOTO 100
         ENDIF

         CALL DZERO(WRK(KMAT),3*NNBASX)

         KPATOM = 0
         NOCOMP = 3
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         DIPORG(1) = MMCORD(1,J)
         DIPORG(2) = MMCORD(2,J)
         DIPORG(3) = MMCORD(3,J)

C        Get 1-electron integrals

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
     &                  (DIPORG(2)-CORD(2,M))**2 +
     &                  (DIPORG(3)-CORD(3,M))**2
                  IF (DISTC .LE. DIST) THEN
                     DIST = DISTC
                     MHIT = M
                  ENDIF
 123           CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
               DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                (DIPORG(2)-QMCOM(2))**2 +
     &                (DIPORG(3)-QMCOM(3))**2
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

C        x-component
         F1=ZERO
         F2=ZERO
         F3=ZERO
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

C 
         CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)

C        Unpack

         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)

C           1-index transformation

            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)

C           Evaluate expectation value

            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMCRO')
            WRK(KEF1+LRI-1+0) = F1
         ELSE
            WRK(KEF1+LRI-1+0) = ZERO
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMCRO')
            WRK(KEF2+LRI-1+0) = F2
         ELSE
            WRK(KEF2+LRI-1+0) = ZERO
         END IF

         IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV1))
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,WRK(KZYM3),WRK(KTLMB),WRK(KTLMA),ISYMV3)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F3,200,'QMMMCRO')
            WRK(KEF3+LRI-1+0) = F3
         ELSE
            WRK(KEF3+LRI-1+0) = ZERO
         END IF

C         y-component
         F1=ZERO
         F2=ZERO
         F3=ZERO
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

C         AO to MO transformation

         CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)

C         Unpack from tri to full

         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
          
         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMCRO')
            WRK(KEF1+LRI-1+1) = F1
         ELSE
            WRK(KEF1+LRI-1+1) = ZERO
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMCRO')
            WRK(KEF2+LRI-1+1) = F2
         ELSE
            WRK(KEF2+LRI-1+1) = ZERO
         END IF

         IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV1))
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,WRK(KZYM3),WRK(KTLMB),WRK(KTLMA),ISYMV3)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F3,200,'QMMMCRO')
            WRK(KEF3+LRI-1+1) = F3
         ELSE
            WRK(KEF3+LRI-1+1) = ZERO
         END IF

C         z-component
         F1=ZERO
         F2=ZERO
         F3=ZERO
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F1,200,'QMMMCRO')
            WRK(KEF1+LRI-1+2) = F1
         ELSE
            WRK(KEF1+LRI-1+2) = ZERO
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,ONE,F2,200,'QMMMCRO')
            WRK(KEF2+LRI-1+2) = F2
         ELSE
            WRK(KEF2+LRI-1+2) = ZERO
         END IF

         IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV1))
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,WRK(KZYM3),WRK(KTLMB),WRK(KTLMA),ISYMV3)
            CALL MELONE(WRK(KTLMA),1,UDV,ONE,F3,200,'QMMMCRO')
            WRK(KEF3+LRI-1+2) = F3
         ELSE
            WRK(KEF3+LRI-1+2) = ZERO
         END IF

         GO TO 100
      ENDIF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KEF1),WRK(KEF1),3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL MPI_REDUCE(WRK(KEF2),WRK(KEF2),3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL MPI_REDUCE(WRK(KEF3),WRK(KEF3),3*NNZAL,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ


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
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>

      DIMENSION WRK(LWRK)
      DIMENSION ZIND1(3*NNZAL),ZIND2(3*NNZAL),ZIND3(3*NNZAL)
      DIMENSION UCMO(*),TRES(*),ZYM2(*),ZYM3(*)

      CALL QENTER('QMMMCRO_M2')

      KTRES = 1
      KLAST = KTRES + N2ORBX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMCRO_M2',-KLAST,LWRK2)

C     Wake up slaves
      IPRTYP = 33
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

C     Send info
C     Integers for arrays/matrixes
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)

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

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)

      CALL MPIXBCAST(IATNOW,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)

      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)

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
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>
#include <gnrinf.h>
#include <infdim.h>

      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)

      CALL QENTER('QMMMCRO_S2')

      QM3 = .TRUE.

      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NATOMS,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2ASHX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NASHT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)

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
      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)

      CALL MPIXBCAST(IATNOW,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*NATOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,NATOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)

      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,NATOMS,'INTEGER',MASTER)
      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
C    
      KUCMO   = 1
      KTRMO   = KUCMO   + NORBT*NBAST
      KUTR    = KTRMO   + NNORBX
      KMAT    = KUTR    + N2ORBX
      KTRES   = KMAT    + 3*NNBASX
      KTLMA   = KTRES   + N2ORBX
      KTLMB   = KTLMA   + N2ORBX
      KZYM2   = KTLMB   + N2ORBX
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

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL DZERO(WRK(KTRES),N2ORBX)

 101  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(J, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(LRI, 1, 'INTEGER', MASTER, MPTAG2)

      IF (J.GT.0) THEN

         DIST2 = (MMCORD(1,J)-QMCOM(1))**2 +
     *            (MMCORD(2,J)-QMCOM(2))**2 +
     *            (MMCORD(3,J)-QMCOM(3))**2
         DIST = SQRT(DIST2)

         IF (DIST .GT. RCUTMM) THEN
            GOTO 101
         ENDIF

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
               DO 124 M=1,NUCIND
                  DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &                    (DIPORG(2)-CORD(2,M))**2 +
     &                    (DIPORG(3)-CORD(3,M))**2
                  IF (DISTC .LE. DIST) THEN
                     DIST = DISTC
                     MHIT = M
                  ENDIF
 124           CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
               DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                (DIPORG(2)-QMCOM(2))**2 +
     &                (DIPORG(3)-QMCOM(3))**2
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

C        x-component
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)
          
         CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         IF (ISYMT.EQ.ISYMV1) THEN 
            F1 = -(1.0D0/2.0D0)*WRK(KIND1+LRI-1+0)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV3,WRK(KZYM3),WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV2))
            CALL DAXPY(N2ORBX,F1,WRK(KTLMB),1,WRK(KTRES),1)
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            F2 = -(1.0D0/2.0D0)*WRK(KIND2+LRI-1+0)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,WRK(KZYM3),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL DAXPY(N2ORBX,F2,WRK(KTLMA),1,WRK(KTRES),1)
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            F3 = -(1.0D0/6.0D0)*WRK(KIND3+LRI-1+0)
            CALL DAXPY(N2ORBX,F3,WRK(KUTR),1,WRK(KTRES),1)
         ENDIF

C        y-component
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         IF (ISYMT.EQ.ISYMV1) THEN
            F1 = -(1.0D0/2.0D0)*WRK(KIND1+LRI-1+1)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV3,WRK(KZYM3),WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV2))
            CALL DAXPY(N2ORBX,F1,WRK(KTLMB),1,WRK(KTRES),1)
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            F2 = -(1.0D0/2.0D0)*WRK(KIND2+LRI-1+1)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,WRK(KZYM3),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL DAXPY(N2ORBX,F2,WRK(KTLMA),1,WRK(KTRES),1)
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            F3 = -(1.0D0/6.0D0)*WRK(KIND3+LRI-1+1)
            CALL DAXPY(N2ORBX,F3,WRK(KUTR),1,WRK(KTRES),1)
         ENDIF

C        z-component
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &              NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         IF (ISYMT.EQ.ISYMV1) THEN
            F1 = -(1.0D0/2.0D0)*WRK(KIND1+LRI-1+2)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV3,WRK(KZYM3),WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV2))
            CALL DAXPY(N2ORBX,F1,WRK(KTLMB),1,WRK(KTRES),1)
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            F2 = -(1.0D0/2.0D0)*WRK(KIND2+LRI-1+2)
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,WRK(KZYM3),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL DAXPY(N2ORBX,F2,WRK(KTLMA),1,WRK(KTRES),1)
         ENDIF

         IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            F3 = -(1.0D0/6.0D0)*WRK(KIND3+LRI-1+2)
            CALL DAXPY(N2ORBX,F3,WRK(KUTR),1,WRK(KTRES),1)
         ENDIF

         GO TO 101
      ENDIF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KTRES),WRK(KTRES),N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ


      CALL QEXIT('QMMMCRO_S2')
      RETURN
      END
