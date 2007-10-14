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
C-----------------------------------------------------------------------
C  /* Deck qm3lno */
      SUBROUTINE QM3LNO(NOSIM,BOVECS,CREF,CMO,XINDX,UDV,DV,EVECS,
     &               WRK,LWRK)
C
C  CBN+JK, Nov. 05
C
C  Purpose:  Calculate SCF/DFT E^2 contribution from a surrounding
C            polarizable MM medium to an orbital trial vector.
C  
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <qm3.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>

      DIMENSION BOVECS(*), CMO(*), XINDX(*), UDV(*), EVECS(KZYVAR,*)
      DIMENSION WRK(*), DV(*), CREF(*)
      LOGICAL LOPEN 

      LOPEN = .FALSE.

      CALL QENTER('QM3LNO')
      IF ( (.NOT. (OLDTG)) .AND. (LOSPC) ) THEN
         CALL QEXIT('QM3LNO')
         RETURN
      END IF

      IF (TRPLET) CALL QUIT('TRIPLET NOT IMPLEMENTED FOR QM/MM YET')

      KUCMO   = 1
      KUBO    = KUCMO   + NORBT*NBAST
      KFSOLMO = KUBO    + NOSIM*N2ORBX
      KUFSOL  = KFSOLMO + NNORBX
      KUTGX   = KUFSOL  + N2ORBX
      KRXYO   = KUTGX   + N2ORBX
      KRRAOx  = KRXYO   + NOSIM*N2ORBX
      KTRMO   = KRRAOx  + NNBASX
      KUTR    = KTRMO   + NNORBX
      KUTRX   = KUTR    + N2ORBX
      KOVLP   = KUTRX   + N2ORBX
      KWRK    = KOVLP   + NOSIM
      LWRK1   = LWRK    - KWRK      

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3LNO',-KWRK,LWRK)

      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KUBO),NOSIM*N2ORBX)
      CALL DZERO(WRK(KFSOLMO),NNORBX)
      CALL DZERO(WRK(KUFSOL),N2ORBX)
      CALL DZERO(WRK(KUTGX),N2ORBX)
      CALL DZERO(WRK(KRXYO),NOSIM*N2ORBX)
      CALL DZERO(WRK(KRRAOx),NNBASX)
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      CALL DZERO(WRK(KUTRX),N2ORBX)
      CALL DZERO(WRK(KOVLP),NOSIM)

      IF (IPRRSP .GE. 40) THEN
         WRITE (LUPRI,'(//A)') ' >>> TEST OUTPUT FROM QM3LNO >>>'
      END IF

      IF (IPRRSP .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> QM3LNO - svec(orb,1) on entry'
         WRITE (LUPRI,'(/A)') ' Z - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Y - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZVAR+KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,
     *               1,1,LUPRI)
      END IF


      CALL UPKCMO(CMO,WRK(KUCMO))

      IF (NOSIM.GT.0) THEN
         CALL RSPZYM(NOSIM,BOVECS,WRK(KUBO))
         CALL DSCAL(NOSIM*N2ORBX,-1.0D0,WRK(KUBO),1)
      END IF

C--------------------------------------------------------
C 1. Construct one-index transformed T^g. 
C    First read qm/mm mo fock matrix from interface file.
C    This includes N_s contribution if OLDTG = .TRUE.. 
C    If OLDTG = .FALSE. N_s is included in one-electron 
C    part of vacuum Hamiltonian
C--------------------------------------------------------

      IF (LUSIFC .LE. 0) THEN
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         LOPEN = .TRUE.
      END IF
     
      REWIND(LUSIFC)

      MMORBT = MAX(4,NNORBT)
      IF (MMORBT .NE. NNORBX) THEN
        CALL QUIT('ERROR IN DIMENSION OF QM/MM FOCK CONTR. IN QM3LNO')
      END IF

      CALL MOLLAB('SOLVTMAT',LUSIFC,LUPRI)
      CALL READT (LUSIFC,NNORBX,WRK(KFSOLMO))

      IF (LOPEN) THEN
         CALL GPCLOSE(LUSIFC,'KEEP')
      END IF

      CALL DSPTSI(NORBT,WRK(KFSOLMO),WRK(KUFSOL))

      DO 200 IOSIM = 1,NOSIM

         IF (IREFSY .EQ. KSYMST) THEN
            WRK(KOVLP-1+IOSIM) = DDOT(KZCONF,CREF,1,EVECS(1,IOSIM),1)
         ELSE
            WRK(KOVLP-1+IOSIM) = 0.0D0
         END IF

         JRXYO = KRXYO + (IOSIM-1)*N2ORBX
         JUBO  = KUBO  + (IOSIM-1)*N2ORBX

         CALL DZERO(WRK(KUTGX),N2ORBX)

         CALL ONEXH1(WRK(JUBO),WRK(KUFSOL),WRK(KUTGX))

         FAC = 1.0D0
         CALL DAXPY(N2ORBX,FAC,WRK(KUTGX),1,WRK(JRXYO),1)

C-------------------------------------------------------
C 2. Construct -alpha \sum_a <0|Q_a|0>Rra
C-------------------------------------------------------
         LUCCEF = -1
         CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND(LUCCEF)

C
C        Readin Rra integrals (ao) and transform to mo
C        ---------------------------------------------
C
         LM = 0

         DO 520 I = 1, ISYTP
           IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
             DO 521 J = NSYSBG(I), NSYSED(I)
               DO 522 K = 1, NUALIS(I)
                 LM = LM + 1

                 FAC = -ALPIMM(I,K)
C x-component
                 CALL DZERO(WRK(KRRAOx),NNBASX)
                 CALL DZERO(WRK(KTRMO),NNORBX)
                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)


                 CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))

                 CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                     NBAST,NORBT)

                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

                 CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))

                 FACx = FAC * QM3QLM(WRK(KUTRX))

                 CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,WRK(JRXYO),1)

C y-component
                 CALL DZERO(WRK(KRRAOx),NNBASX)
                 CALL DZERO(WRK(KTRMO),NNORBX)
                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)

                 CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                 CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                     NBAST,NORBT)
                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                 CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))

                 FACy = FAC * QM3QLM(WRK(KUTRX))

                 CALL DAXPY(N2ORBX,FACy,WRK(KUTR),1,WRK(JRXYO),1)
C z-component
                 CALL DZERO(WRK(KRRAOx),NNBASX)
                 CALL DZERO(WRK(KTRMO),NNORBX)
                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)

                 CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                 CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                     NBAST,NORBT)
                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                 CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))

                 FACz = FAC * QM3QLM(WRK(KUTRX))

                 CALL DAXPY(N2ORBX,FACz,WRK(KUTR),1,WRK(JRXYO),1)

  522       CONTINUE
  521     CONTINUE
        END IF
  520 CONTINUE
C
      CALL GPCLOSE(LUCCEF,'KEEP')

  200 CONTINUE

C--------------------------------------------
C 3.  Add effective operators to response
C--------------------------------------------

      CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYO))

      IF (IPRRSP .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> QM3LNO - svec(orb,1) on exit'
         WRITE (LUPRI,'(/A)') ' Z - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Y - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZVAR+KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,
     *               1,1,LUPRI)
      END IF

      DO 900 IOSIM = 1,NOSIM
         TEST = DDOT(KZCONF,CREF,1,EVECS(1,IOSIM),1)
     *        - WRK(KOVLP-1+IOSIM)
         IF (ABS(TEST) .GT. 1.D-8) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(/A,I5,/A,1P,D12.4)')
     *         ' >>> SLVLNO WARNING, for IOSIM =',IOSIM,
     *         ' <CREF | SVEC_solvent(iosim) > =',TEST
         END IF
  900    CONTINUE

      CALL QEXIT('QM3LNO')
      RETURN
      END 
C
C-----------------------------------------------------------------------
C  /* Deck qm3qro */
      SUBROUTINE QM3QRO(VEC1,VEC2,ETRS,XINDX,ZYM1,ZYM2,
     &                 UDV,WRK,LWRK,KZYVR,KZYV1,KZYV2,
     &                 IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP)
C
C  CBN+JK, Dec. 05
C
C  Purpose:  Calculate SCF/DFT E^3 contribution from a surrounding
C            polarizable MM medium to an orbital trial vector.
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
#include <qm3.h>
#include <priunit.h>
#include <ccinftap.h>

      PARAMETER ( D0=0.0D0, D1=1.0D0 )
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(LWRK),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2)
      DIMENSION MJWOP(2,MAXWOP,8)
      LOGICAL LCON, LORB, LREF, LOPEN

      CALL QENTER('QM3QRO')
      IF ( (.NOT. (OLDTG)) .AND. (LOSPC) ) THEN
         CALL QEXIT('QM3QRO')
         RETURN
      END IF

      LOPEN = .FALSE.
  
      IF (TRPLET) CALL QUIT('TRIPLET NOT IMPLEMENTED FOR QM/MM YET')

      KCREF   = 1
      KTRES   = KCREF   + NCREF
      KUCMO   = KTRES   + N2ORBX
      KFSOLMO = KUCMO   + NORBT*NBAST
      KUFSOL  = KFSOLMO + NNORBX
      KTLMA   = KUFSOL  + N2ORBX
      KTLMB   = KTLMA   + N2ORBX
      KRRAOx  = KTLMB   + N2ORBX
      KTRMO   = KRRAOx  + NNBASX
      KUTR    = KTRMO   + NNORBX
      KWRK    = KUTR    + N2ORBX
      LWRK1   = LWRK    - KWRK 

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3QRO',-KWRK,LWRK1)

      CALL DZERO(WRK(KCREF),NCREF)
      CALL DZERO(WRK(KTRES),N2ORBX)
      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KFSOLMO),NNORBX)
      CALL DZERO(WRK(KUFSOL),N2ORBX)
      CALL DZERO(WRK(KTLMA),N2ORBX)
      CALL DZERO(WRK(KTLMB),N2ORBX)
      CALL DZERO(WRK(KRRAOx),NNBASX)
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      NSIM = 1

C NB NB NB We assume no symmetry in the DFT/MM calculations

      ISYMT = 1

C Get the reference state

      CALL GETREF(WRK(KCREF),MZCONF(1))


C Unpack the response vectors

      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)


C Unpack symmetry blocked CMO

      CALL UPKCMO(CMO,WRK(KUCMO))

C--------------------------------------------------------
C 1. Construct two-index transformed T^g.
C    First read qm/mm mo fock matrix from interface file.
C    This includes N_s contribution if OLDTG = .TRUE..
C    If OLDTG = .FALSE. N_s is included in one-electron
C    part of vacuum Hamiltonian
C--------------------------------------------------------

      IF (LUSIFC .LE. 0) THEN
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         LOPEN = .TRUE.
      END IF

      MMORBT = MAX(4,NNORBT)
      IF (MMORBT .NE. NNORBX) THEN
        CALL QUIT('ERROR IN DIMENSION OF QM/MM FOCK CONTR. IN QM3QRO')
      END IF
C
      REWIND(LUSIFC)

      CALL MOLLAB('SOLVTMAT',LUSIFC,LUPRI)
      CALL READT (LUSIFC,NNORBX,WRK(KFSOLMO))

      IF (LOPEN) THEN
         CALL GPCLOSE(LUSIFC,'KEEP')
      END IF

      CALL DSPTSI(NORBT,WRK(KFSOLMO),WRK(KUFSOL))

C     Create the effective operator:
C
C     TRES = Tg(k1,k2) + A1(k2) + A12

      IF (ISYMT.EQ.ISYMV1 .AND. ISYMT.EQ.ISYMV2
     &    .AND. ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
          CALL DZERO(WRK(KTLMA),N2ORBX)
          CALL DZERO(WRK(KTLMB),N2ORBX)
          CALL OITH1(ISYMV1,ZYM1,WRK(KUFSOL),WRK(KTLMA),ISYMT)
          CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &              MULD2H(ISYMT,ISYMV1))
          CALL DAXPY(N2ORBX,0.5D0,WRK(KTLMB),1,WRK(KTRES),1)
      END IF

      CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUCCEF)

C     Readin Rra integrals (ao) and transform to mo
C     ---------------------------------------------

      LM = 0

      DO 520 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 521 J = NSYSBG(I), NSYSED(I)
            DO 522 K = 1, NUALIS(I)
              LM = LM + 1

              FAC = -ALPIMM(I,K)

C             A1(k2) = F1*Rra(k2)
C             A12 = F2*Rra
C
C             F1 = 2*<0| Rra(k1) |0>
C             F2 = <0| Rra(k1,k2) |0>

C x-component
              F1=D0
              F2=D0
              CALL DZERO(WRK(KRRAOx),NNBASX)
              CALL DZERO(WRK(KTRMO),NNORBX)
              CALL DZERO(WRK(KUTR),N2ORBX)
C
              CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
              CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                     NBAST,NORBT)
              CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
C
              IF (ISYMT.EQ.ISYMV1) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL MELONE(WRK(KTLMA),1,UDV,D1,F1,200,'QM3QRO')
                 F1 = 1.0D0*F1*FAC

                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)

                 CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)
              END IF

              IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL DZERO(WRK(KTLMB),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
                 CALL MELONE(WRK(KTLMB),1,UDV,D1,F2,200,'C3SOL')
                 F2 = 0.5D0*F2*FAC

                 CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)
              END IF

C y-component
              F1=D0
              F2=D0
              CALL DZERO(WRK(KRRAOx),NNBASX)
              CALL DZERO(WRK(KTRMO),NNORBX)
              CALL DZERO(WRK(KUTR),N2ORBX)

              CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
              CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                  NBAST,NORBT)
              CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

              IF (ISYMT.EQ.ISYMV1) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL MELONE(WRK(KTLMA),1,UDV,D1,F1,200,'QM3QRO')
                 F1 = 1.0D0*F1*FAC

                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)

                 CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)
              END IF

              IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL DZERO(WRK(KTLMB),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
                 CALL MELONE(WRK(KTLMB),1,UDV,D1,F2,200,'QM3QRO')
                 F2 = 0.5D0*F2*FAC

                 CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)
              END IF

C z-component
              F1=D0
              F2=D0
              CALL DZERO(WRK(KRRAOx),NNBASX)
              CALL DZERO(WRK(KTRMO),NNORBX)
              CALL DZERO(WRK(KUTR),N2ORBX)

              CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
              CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                  NBAST,NORBT)
              CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

              IF (ISYMT.EQ.ISYMV1) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL MELONE(WRK(KTLMA),1,UDV,D1,F1,200,'C3SOL')
                 F1 = 1.0D0*F1*FAC

                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)

                 CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)   
              END IF

              IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL DZERO(WRK(KTLMB),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
                 CALL MELONE(WRK(KTLMB),1,UDV,D1,F2,200,'C3SOL')
                 F2 = 0.5D0*F2*FAC

                 CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)
              END IF

  522       CONTINUE
  521     CONTINUE
        END IF
  520 CONTINUE

      CALL GPCLOSE(LUCCEF,'KEEP')

C
C       Make the gradient
C
C     / <0| [qj ,TRES] |0> \
C     |          0         |
C     | <0| [qj+,TRES] |0> |
C      \         0         /
C
      ISYMDN = 1
      OVLAP  = D1
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
      CALL QEXIT('QM3QRO')
      RETURN
      END
C*****************************************************************************************
C  /* Deck qm3qlm */
      FUNCTION QM3QLM(RLM)
C
C QUADRATIC PACKED DENSITY MATRIX AND ACTIVE INTEGRALS
C
#include <implicit.h>
      DIMENSION RLM(NORBT,*)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
C
C Used from common blocks:
C  INFORB: NISHT,NASHT
C  INFIND: IROW(*), ISX(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <priunit.h>

      TELM  = D0

      DO 300 IW = 1,NISHT
         I = ISX(IW)
         TELM = TELM + D2*RLM(I,I)
  300 CONTINUE

      QM3QLM = TELM
      RETURN
      END
C*********************************************************************************************
