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
C      CBN+JK
C*************************************************************************
C  /* Deck qm3fck */
      SUBROUTINE QM3FCK(DCAO,DVAO,FSOL,EQM3T,ENSEL,EPOL,EELEL,
     *                  WRK,LWRK,IPRINT)

#include <implicit.h>
#include <dummy.h>
#include <inftap.h>
#include <priunit.h>
#include <qm3.h>
#include <mxcent.h>
#include <thrzer.h>
#include <iratdef.h>
#include <codata.h>
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>

      DIMENSION DCAO(*), DVAO(*)
      DIMENSION FSOL(*), WRK(LWRK)
      LOGICAL HSFOCK
      DIMENSION DFTNS(MXQM3)


      CALL QENTER('QM3FCK')
      KDTAO = 1
      KTAO = KDTAO + NNBASX
      KEND = KTAO + NNBASX
      LWRK1 = LWRK - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3FCK',-KEND,LWRK)

      HSFOCK=((DODFT .AND. NASHT .GE. 1) .OR. HSROHF)

C     Get total density      
      IF (NASHT .EQ. 0) THEN
            CALL PKSYM1(WRK(KDTAO),DCAO,NBAS,NSYM,-1)
      ELSE
            DO I = 1,NNBAST
               IF (HSFOCK) THEN
                  WRK(KTAO-1+I) = DCAO(I)
               ELSE
                  WRK(KTAO-1+I) = DCAO(I) + DVAO(I)
               END IF
            END DO
            CALL PKSYM1(WRK(KDTAO),WRK(KTAO),NBAS,NSYM,-1)
      END IF

C     From COMMON
      ENSQM3 = 0.0D0
      EPOQM3 = 0.0D0

C     Calculate RRa, Obar and ENSA and write these to file
      CALL QM3_RRA(WRK(KDTAO),WRK(KEND),LWRK1,IPRINT) 

C     Calculate Ns and keep in memory (DFTNS). 
      CALL QM3_NS(WRK(KDTAO),DFTNS,WRK(KEND),LWRK1,IPRINT)

C     Modify the fock operator. Modification returned in FSOL. 
      CALL QM3_FMO(FSOL,WRK(KEND),LWRK1,IPRINT)

C     Calculate the QM3 contribution to the energy (returned in EQM3T)
      CALL QM3_ENERGY(DFTNS,ENSEL,EPOL,EELEL,EQM3T,WRK(KEND),LWRK1,
     &                IPRINT)

      CALL QEXIT('QM3FCK')
      RETURN
      END
C ***********************************************************************
C  /* Deck qm3_rra */
      SUBROUTINE QM3_RRA(DCAO,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>

      DIMENSION WRK(LWRK)
      DIMENSION EELEC(3,MXQM3)
      DIMENSION FFJ(3)
      LOGICAL   LOINDM

      IF (LOSPC) RETURN

      CALL QENTER('QM3_RRA')

      KRRAOx = 1
      KRRAOy = KRRAOx + NNBASX
      KRRAOz = KRRAOy + NNBASX
      KRAx   = KRRAOz + NNBASX
      KRAy   = KRAx   + NCOMS
      KRAz   = KRAy   + NCOMS
      KWRK1  = KRAz   + NCOMS
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3_RRA',-KWRK1,LWRK)

      CALL DZERO(WRK(KRRAOx),3*NNBASX)
      CALL DZERO(WRK(KRAx),3*NCOMS)

      DO 879 I = 1, MXQM3
        DO 880 J = 1, 3
          EELEC(J,I) = 0.0D0
  880   CONTINUE
  879 CONTINUE

      IF (IQM3PR .GT. 15) THEN
        WRITE(LUPRI,*) 'QM3RAINT: Read in integrals'
      ENDIF

      CALL GPOPEN(LUQM3E,'ELFDMM','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUQM3E)

      LM = 0

      DO 520 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 521 J = NSYSBG(I), NSYSED(I)
            DO 522 K = 1, NUALIS(I)
              LM = LM + 1

              CALL READT(LUQM3E,NNBASX,WRK(KRRAOx))
              WRK(KRAx+LM-1) = DDOT(NNBASX,WRK(KRRAOx),1,DCAO,1)

              CALL READT(LUQM3E,NNBASX,WRK(KRRAOy))
              WRK(KRAy+LM-1) = DDOT(NNBASX,WRK(KRRAOy),1,DCAO,1)

              CALL READT(LUQM3E,NNBASX,WRK(KRRAOz))
              WRK(KRAz+LM-1) = DDOT(NNBASX,WRK(KRRAOz),1,DCAO,1)

  522       CONTINUE
  521     CONTINUE
        END IF
  520 CONTINUE

      IF (IQM3PR .GE. 5) THEN
        WRITE(LUPRI,'(/A/A/A)')
     *  ' +==========================================================+',
     *  ' | COM| <Rra>_x         | <Rra>_y         | <Rra>_z         |',
     *  ' +==========================================================+'

         NUM = 0

         DO 215 I = 1, ISYTP
           IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
             DO 216 J = NSYSBG(I), NSYSED(I)
               DO 217 K=1,NUALIS(I)

                 NUM = NUM + 1

                 WRITE(LUPRI,'(A,I3,A,3(F16.10,A)/A)')
     *  ' | ', J,'|', WRK(KRAx + NUM - 1),' |',
     *           WRK(KRAy + NUM - 1),' |', WRK(KRAz + NUM - 1),' |',
     *  ' +----------------------------------------------------------+'
  217          CONTINUE
  216        CONTINUE
           END IF
  215    CONTINUE
         WRITE(LUPRI,'(//,A)')
      END IF

      IF (LM .NE. NCOMS) THEN
        CALL QUIT('Error in no. of center of masses in QM3RAINT')
      END IF

      DO 534 LM = 1,NCOMS
        EELEC(1,LM) = WRK(KRAx+LM-1)
        EELEC(2,LM) = WRK(KRAy+LM-1)
        EELEC(3,LM) = WRK(KRAz+LM-1)
  534 CONTINUE

C     If RELMOM is true we want to include the external field(s) in
C     the determination of the induced dipole moments

        FFJ(1) = 0.0D0
        FFJ(2) = 0.0D0
        FFJ(3) = 0.0D0

      IF (RELMOM) THEN
        DO 330 IF =1, NFIELD
          IF (LFIELD(IF) .EQ. 'XDIPLEN') FFJ(1) = FFJ(1) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'YDIPLEN') FFJ(2) = FFJ(2) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'ZDIPLEN') FFJ(3) = FFJ(3) + EFIELD(IF)
  330   CONTINUE
      END IF

      IF (FIXMOM) THEN
        WRITE(LUPRI,'(/A)')'FIXMOM: NO ITER. DETERM. OF MYIND'
      ELSE IF (LOSPC) THEN
        WRITE(LUPRI,'(/A)')'All MM models are SPC, INDMOM not called'
      ELSE
        LOINDM = .FALSE.
        CALL INDMOM(EELEC,LOINDM,FFJ)
      END IF

      CALL QM3_OBAR(FFJ)

      CALL CC_PUT31('CC_RA',NCOMS,WRK(KRAx),WRK(KRAy),WRK(KRAz))

      CALL GPCLOSE(LUQM3E,'KEEP')

      CALL QEXIT('QM3_RRA')
      RETURN
      END
C*********************************************************************
C  /* Deck qm3_ns */
      SUBROUTINE QM3_NS(DCAO,DFTNS,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>

      DIMENSION WRK(LWRK) 
      DIMENSION DFTNS(*)
      
      CALL QENTER('QM3_NS')

      KNSAO  = 1
      KWRK1  = KNSAO + NNBASX
      LWRK1  = LWRK  - KWRK1

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3_NS',-KWRK1,LWRK)

      CALL DZERO(WRK(KNSAO),NNBASX)

      ENSEL = 0.0D0
       
      CALL GPOPEN(LUQM3P,'POTMM','UNKNOWN',' ',
     &          'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND (LUQM3P)

      FAC1 = -1.0D0

      L = 0

      DO 525 I = 1, ISYTP
        IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
          DO 526 J = NSYSBG(I), NSYSED(I)
            DO 527 K = 1,NSISY(I)

              L = L +1

              CALL READT(LUQM3P,NNBASX,WRK(KNSAO))

              DFTNS(L) = DDOT(NNBASX,DCAO,1,WRK(KNSAO),1)
              ENSEL     = ENSEL - DFTNS(L) 

              IF (IQM3PR .GE. 21) THEN
                WRITE (LUPRI,'(/A,I3,A)')
     *             ' N(',L,')_ao matrix: '
                CALL OUTPAK(WRK(KNSAO),NBAST,1,LUPRI)
              END IF

  527       CONTINUE
  526     CONTINUE
        END IF
  525 CONTINUE       

      CALL GPCLOSE(LUQM3P,'KEEP')

C-------------------
C     Print section:
C-------------------
C
      IF (IQM3PR .GT. 3) THEN
        WRITE(LUPRI,'(/A/A/A)')
     *        ' +======================+',
     *        ' |Site| <N_s>           |',
     *        ' +======================+'

        LS = 0

        DO 215 I = 1, ISYTP
          IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
            DO 216 J = NSYSBG(I), NSYSED(I)
              DO 217 K = 1, NSISY(I)

                LS = LS + 1

                WRITE(LUPRI,'(A,I3,A,F16.10,A/A)')
     *                      ' | ', LS,'|', DFTNS(LS),' |',
     *                      ' +----------------------+'
  217         CONTINUE
  216       CONTINUE
          END IF
  215   CONTINUE
        WRITE(LUPRI,'()')
      END IF

      CALL QEXIT('QM3_NS')
      RETURN
      END
C*********************************************************************
C  /* Deck qm3_fmo */
      SUBROUTINE QM3_FMO(FSOL,WRK,LWRK,IPRINT)   

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>

      DIMENSION WRK(LWRK)
      DIMENSION FSOL(*)

      IF (LOSPC .AND. (.NOT. OLDTG)) RETURN

      CALL QENTER('QM3_FMO')

      KRRAOx = 1
      KRRAOy = KRRAOx + NNBASX
      KRRAOz = KRRAOy + NNBASX
      KRAx   = KRRAOz + NNBASX
      KRAy   = KRAx   + NCOMS
      KRAz   = KRAy   + NCOMS
      KENSAx = KRAz   + NCOMS
      KENSAy = KENSAx + NCOMS
      KENSAz = KENSAy + NCOMS
      KTAO   = KENSAz + NCOMS
      KWRK1  = KTAO   + NNBASX
      LWRK1  = LWRK   - KWRK1

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3_FMO',-KWRK1,LWRK)

      CALL DZERO(WRK(KRRAOx),3*NNBASX)
      CALL DZERO(WRK(KRAx),6*NCOMS)
      CALL DZERO(WRK(KTAO),NNBASX)

      IF (.NOT. LOSPC) THEN  

        CALL CC_GET31('CC_RA',NCOMS,
     *                 WRK(KRAx),WRK(KRAy),WRK(KRAz))

        CALL CC_GET31('ENSAFILE',NCOMS,
     *                 WRK(KENSAx),WRK(KENSAy),WRK(KENSAz))

        IF (IQM3PR .GT. 5) THEN
          WRITE(LUPRI,'(/A/A/A)')
     *  ' +==========================================================+',
     *  ' | COM| <Rra>_x         | <Rra>_y         | <Rra>_z         |',
     *  ' +==========================================================+'

          NUM = 0

          DO 215 I = 1, ISYTP
            IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
              DO 216 J = NSYSBG(I), NSYSED(I)
                DO 217 K=1,NUALIS(I)

                  NUM = NUM + 1

                  WRITE(LUPRI,'(A,I3,A,F16.10,A,F16.10,A,F16.10,A/A)')
     *  ' | ', J,'|', WRK(KRAx + NUM - 1),' |',
     *              WRK(KRAy + NUM - 1),' |', WRK(KRAz + NUM - 1),' |',
     *  '+----------------------------------------------------------+'
  217           CONTINUE
  216         CONTINUE
            END IF
  215     CONTINUE

          WRITE(LUPRI,'(/A/A/A)')
     *  ' +==========================================================+',
     *  ' | COM| ENSA_x          | ENSA_y          | ENSA_z          |',
     *  ' +==========================================================+'

          NUM = 0

          DO 415 I = 1, ISYTP
            IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN

              NUM = NUM + 1

              DO 416 J = NSYSBG(I), NSYSED(I)
                DO 417 K=1,NUALIS(I)
                  WRITE(LUPRI,'(A,I3,A,F16.10,A,
     &                          F16.10,A,F16.10,A/A)')
     *  ' | ', J,'|', WRK(KENSAx + NUM - 1),' |',
     *           WRK(KENSAy + NUM - 1),' |', WRK(KENSAz + NUM - 1),' |',
     *  ' +----------------------------------------------------------+'
  417           CONTINUE
  416         CONTINUE
            END IF
  415     CONTINUE
        END IF

        CALL GPOPEN(LUQM3E,'ELFDMM','OLD',' ',
     &              'UNFORMATTED',IDUMMY,.FALSE.)
        REWIND(LUQM3E)

        LM = 0

        DO 520 I = 1, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 521 J = NSYSBG(I), NSYSED(I)
              DO 522 K = 1, NUALIS(I)
                LM = LM + 1

                CALL READT(LUQM3E,NNBASX,WRK(KRRAOx))

                IF (IQM3PR .GE. 21) THEN
                  WRITE (LUPRI,'(/A)') ' Rra_ao x matrix:'
                  CALL OUTPAK(WRK(KRRAOx),NBAST,1,LUPRI)
                END IF
  
                CALL READT(LUQM3E,NNBASX,WRK(KRRAOy))
  
                IF (IQM3PR .GE. 21) THEN
                  WRITE (LUPRI,'(/A)') ' Rra_ao y matrix:'
                  CALL OUTPAK(WRK(KRRAOy),NBAST,1,LUPRI)
                END IF

                CALL READT(LUQM3E,NNBASX,WRK(KRRAOz))
  
                IF (IQM3PR .GE. 21) THEN
                  WRITE (LUPRI,'(/A)') ' Rra_ao z matrix:'
                  CALL OUTPAK(WRK(KRRAOz),NBAST,1,LUPRI)
                END IF

                IF (MDLWRD(I) .EQ. 'SPC_E01') THEN
  
                  FACx =  - ALPIMM(I,K) * (WRK(KRAx + LM - 1)
     &                  + 0.5D0 * WRK(KENSAx + LM - 1))
                  FACy =  - ALPIMM(I,K) * (WRK(KRAy + LM - 1)
     &                  + 0.5D0 * WRK(KENSAy + LM - 1))
                  FACz =  - ALPIMM(I,K) * (WRK(KRAz + LM - 1)
     &                  + 0.5D0 * WRK(KENSAz + LM - 1))

                  CALL DAXPY(NNBASX,FACx,WRK(KRRAOx),1,WRK(KTAO),1)
  
                  CALL DAXPY(NNBASX,FACy,WRK(KRRAOy),1,WRK(KTAO),1)
 
                  CALL DAXPY(NNBASX,FACz,WRK(KRRAOz),1,WRK(KTAO),1)
                END IF
  522         CONTINUE
  521       CONTINUE
          END IF
  520   CONTINUE

        CALL GPCLOSE(LUQM3E,'KEEP')

        IF (IQM3PR .GE. 15) THEN
            WRITE (LUPRI,'(/A)') ' QM/MM contribution:'
            CALL OUTPAK(WRK(KTAO),NBAST,1,LUPRI)
        END IF

      END IF ! LOSPC

      IF (OLDTG) THEN

C       We use RRAOx in this case for the pot.energy integrals
        CALL DZERO(WRK(KRRAOx),NNBASX)

        CALL GPOPEN(LUQM3P,'POTMM','UNKNOWN',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
        REWIND (LUQM3P)

        FAC1 = -1.0D0

        L = 0

        DO 525 I = 1, ISYTP
          IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
            DO 526 J = NSYSBG(I), NSYSED(I)
              DO 527 K = 1,NSISY(I)

                L = L +1

                CALL READT(LUQM3P,NNBASX,WRK(KRRAOx))

                IF (IQM3PR .GE. 21) THEN
                  WRITE (LUPRI,'(/A,I3,A)')
     *               ' N(',L,')_ao matrix: '
                  CALL OUTPAK(WRK(KRRAOx),NBAST,1,LUPRI)
                END IF

                CALL DAXPY(NNBASX,FAC1,WRK(KRRAOx),
     *                     1,WRK(KTAO),1)
  527         CONTINUE
  526       CONTINUE
          END IF
  525   CONTINUE

        CALL GPCLOSE(LUQM3P,'KEEP')

      END IF

      CALL PKSYM1(WRK(KTAO),FSOL,NBAS,NSYM,+1)

      CALL QEXIT('QM3_FMO')
      RETURN
      END
C*********************************************************************
C  /* Deck qm3_energy */
      SUBROUTINE QM3_ENERGY(DFTNS,EEL,EPOL,EELEL,EQM3T,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>

      DIMENSION WRK(LWRK)
      DIMENSION FFJ(3), DFTNS(*)

      CALL QENTER('QM3_ENERGY')

      IF ( .NOT. (LOSPC) ) THEN

        KOMMSx = 1
        KOMMSy = KOMMSx + NCOMS
        KOMMSz = KOMMSy + NCOMS
        KRAx   = KOMMSz + NCOMS
        KRAy   = KRAx   + NCOMS
        KRAz   = KRAy   + NCOMS
        KWRK1  = KRAz   + NCOMS
        LWRK1  = LWRK   - KWRK1

        IF (LWRK1.LT.0) CALL QUIT( 'Too little work in QM3_ENERGY')

        CALL DZERO(WRK(KOMMSx),6*NCOMS)

      END IF

C     First we see if any static fields are to be added

        FFJ(1) = 0.0D0
        FFJ(2) = 0.0D0
        FFJ(3) = 0.0D0

      IF (RELMOM) THEN
        DO 330 IF =1, NFIELD
          IF (LFIELD(IF) .EQ. 'XDIPLEN') FFJ(1) = FFJ(1) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'YDIPLEN') FFJ(2) = FFJ(2) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'ZDIPLEN') FFJ(3) = FFJ(3) + EFIELD(IF)
  330   CONTINUE
      END IF

C-------------------------------------
C     Add up the energy contributions:
C     1) Epol:
C-------------------------------------

      EQM3T = 0.0D0

      IF (.NOT. LOSPC) THEN
        CALL CC_GET31('OBARFILE',NCOMS,
     *                WRK(KOMMSx),WRK(KOMMSy),WRK(KOMMSz))

        CALL CC_GET31('CC_RA',NCOMS,
     *                WRK(KRAx),WRK(KRAy),WRK(KRAz))

        KS = 0

        DO 110 I = 1, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 120 J = NSYSBG(I), NSYSED(I)
              DO 121 K = 1, NUALIS(I)
                KS = KS+1
                RAT = 0.0D0
                RAT =  0.5D0 * ALPIMM(I,K) * WRK(KRAx + KS - 1) *
     &                (WRK(KRAx + KS - 1) + WRK(KOMMSx + KS - 1))
     &               + 0.5D0 * ALPIMM(I,K) * WRK(KRAy + KS - 1) *
     &                (WRK(KRAy + KS - 1) + WRK(KOMMSy + KS - 1))
     &               + 0.5D0 * ALPIMM(I,K) * WRK(KRAz + KS - 1) *
     &                (WRK(KRAz + KS - 1) + WRK(KOMMSz + KS - 1))
                EQM3T = EQM3T - RAT
  121         CONTINUE
  120       CONTINUE
          END IF
  110   CONTINUE

        TEMP1 = EQM3T

        CALL QM3_OTILDE(OTILDE,FFJ)

      ELSE

        OTILDE = 0.0D0
        TEMP1 = 0.0D0

      END IF

      EQM3T = EQM3T + OTILDE
      EPOL = EQM3T

C--------
C 2) Eel:
C--------
C
      ETEMP = 0.0D0
      L = 0

      DO 130 I = 1, ISYTP
        IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
          DO 140 J = NSYSBG(I), NSYSED(I)
            DO 150 K = 1, NSISY(I)
              L = L + 1
              ETEMP = ETEMP - DFTNS(L)
  150       CONTINUE
  140     CONTINUE
        END IF
  130 CONTINUE

      EELEL = ETEMP

      EQM3T = EQM3T + ENUQM3 + EELEL

      EEL = ENUQM3 + EELEL

C-------------
C 3) E(QM/MM):
C-------------

      EQM3T = EQM3T + ECLVDW

C-----------------------------------------------
C   Testing the energy contributions one by one:
C-----------------------------------------------

      IF (IQM3PR .GT. 1) THEN
        WRITE(LUPRI,'(//A/A/A)')
     *       ' +====================================='
     *      //'======================================'
     *      //'==================================+',

     *       ' |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
     *      //'The different energy contributions outlined'
     *      //' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|',

     *       ' +====================================='
     *      //'======================================'
     *      //'==================================+'
        WRITE(LUPRI,'(A)')
     *       ' | Evdw                | E-nuclear           |'
     *      //' -Sum_s <N_s>        |'
     *      //' -alpha*Sum_a[    >  | O-tilde             |',

     *       ' | contribution        | contribution        |'
     *      //'                     |'
     *      //'<Rra>*{<Rra>+OmmS}]  | contribution        |',

     *       ' +-------------------------------------'
     *      //'--------------------------------------'
     *      //'----------------------------------+'
        WRITE(LUPRI,'(A,F20.16,A,F20.16,A,F20.16,A,
     &                F20.16,A,F20.16,A)')
     *       ' |',ECLVDW,' |',ENUQM3,' |',EELEL,' |',TEMP1,' |',
     *      OTILDE,' |'
        WRITE(LUPRI,'(A/)')
     *       ' +====================================='
     *      //'======================================'
     *      //'==================================+'
      END IF

      CALL QEXIT('QM3_ENERGY')
      RETURN
      END
C*************************************************************************
C  /* Deck qm3grad */
      SUBROUTINE QM3GRAD(CREF,CMO,INDXCI,G,EQM3T,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <inftap.h>
#include <infpri.h>
#include <inflin.h>

      DIMENSION CREF(*),CMO(*),INDXCI(*)
      DIMENSION G(*),WRK(LWRK)
      LOGICAL FIRST, FNDLAB
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0, DCVAL = D2 )
      CHARACTER*8 STAR8, SOLVDI, EODATA
      SAVE        FIRST
      DATA        FIRST/.TRUE./, STAR8/'********'/
      DATA        SOLVDI/'SOLVDIAG'/, EODATA/'EODATA  '/

      CALL QENTER('QM3GRAD')

      KDIASH  = 1
      KGRDLM  = KDIASH  + NVAR
      KUCMO   = KGRDLM  + NVARH
      KFSOLAO = KUCMO   + NORBT*NBAST
      KFSOLMO = KFSOLAO + NNBASX
      KDIALM  = KFSOLMO + NNORBX
      KDV     = KDIALM  + NVAR
      KDENS   = KDV     + N2BASX
      KDVS    = KDENS   + NNBASX
      KWRK    = KDVS    + NNBASX
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3GRAD',-KWRK,LWRK)

C     Construct the ao density matrix from the molecular orbial
C     coefficients.

      CALL FCKDEN((NISHT.GT.0),.FALSE.,WRK(KDV),
     *            DUMMY,CMO,DUMMY,WRK(KWRK),LWRK1)

      CALL DGEFSP(NBAST,WRK(KDV),WRK(KDVS))
      CALL PKSYM1(WRK(KDVS),WRK(KDENS),NBAS,NSYM,1)
    

      CALL DZERO(WRK(KDIASH),NVAR)
      CALL DZERO(WRK(KGRDLM),NVARH)

      EQM3T = 0.0D0
      ENSEL = 0.0D0
      EPOL = 0.0D0
      EELEL = 0.0D0

C     Make effective solvent operator and transform to mo basis
      CALL QM3FCK(WRK(KDENS),WRK(KDENS),WRK(KFSOLAO),EQM3T,ENSEL,
     &            EPOL,EELEL,WRK(KWRK),LWRK1,IPRINT)
      CALL UPKCMO(CMO,WRK(KUCMO))
      CALL UTHU(WRK(KFSOLAO),WRK(KFSOLMO),WRK(KUCMO),WRK(KWRK),
     &             NBAST,NORBT)
 
      IF (.NOT. OLDTG) EQM3T = EQM3T - EELEL 

C     Calculate gradients
      IF (NWOPT .GT. 0) THEN
         CALL SOLGO(D2,WRK(KDENS),WRK(KFSOLMO),WRK(KGRDLM+NCONF))
      ENDIF

      IF (NWOPPT .GT. 0) THEN
         CALL OR1DIA(WRK(KDENS),WRK(KFSOLMO),WRK(KDIALM+NCONST),
     &               IPRINT)
      END IF

      DO 420 I = 0,(NVAR-1)
         WRK(KDIASH+I) = WRK(KDIASH+I)
     &                 + WRK(KDIALM+I)
     &                 + WRK(KGRDLM+I) * WRK(KGRDLM+I)
  420    CONTINUE

      TEST = DDOT(NCONF,CREF,1,WRK(KGRDLM),1)
      IF (ABS(TEST) .GT. 1.D-8) THEN
          NWARN = NWARN + 1
          WRITE (LUPRI,'(/A,I5,/A,1P,D12.4)')
     &      ' >>> QM3GRD WARNING <CREF | GRAD > =',TEST
         END IF

      FAC=1.0D0
      CALL DAXPY(NVARH,FAC,WRK(KGRDLM),1,G,1)

      IF (LUIT2 .GT. 0) THEN
         NW4 = MAX(NWOPT,4)
         REWIND LUIT2
         IF (FNDLAB(EODATA,LUIT2)) BACKSPACE LUIT2
         WRITE (LUIT2) STAR8,STAR8,STAR8,SOLVDI
         IF (NWOPT .GT. 0) CALL WRITT(LUIT2,NW4,WRK(KDIASH+NCONF))
         WRITE (LUIT2) STAR8,STAR8,STAR8,EODATA
      END IF

      FIRST = .FALSE.
      CALL QEXIT('QM3GRAD')
      RETURN
      END
C
C*************************************************************************
C  /* Deck qm3lin */
      SUBROUTINE QM3LIN(NOSIM,BOVECS,CREF,CMO,INDXCI,
     &                  SOVECS,ORBLIN,WRK,LWRK)
C
#include <implicit.h>
#include <inflin.h>
#include <inforb.h>
#include <dummy.h>
#include <priunit.h>
#include <infpri.h>

      DIMENSION BOVECS(*),CREF(*),CMO(*),INDXCI(*)
      DIMENSION SOVECS(*),WRK(LWRK)
      LOGICAL   ORBLIN

      CALL QENTER('QM3LIN')

      KDV   = 1
      KDVS  = KDV   + N2BASX
      KDENS = KDVS  + NNBASX
      KWRK  = KDENS + NNBASX
      LWRK1 = LWRK  - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3LIN',-KWRK,LWRK)

C     Construct the ao density matrix from the molecular orbial
C     coefficients.
      CALL FCKDEN((NISHT.GT.0),.FALSE.,WRK(KDV),
     *            DUMMY,CMO,DUMMY,WRK(KWRK),LWRK1)

      CALL DGEFSP(NBAST,WRK(KDV),WRK(KDVS))
      CALL PKSYM1(WRK(KDVS),WRK(KDENS),NBAS,NSYM,1)

C
      IF ( NOSIM .GT.0 ) THEN
         IF ( .NOT.ORBLIN ) THEN
            NSVAR  = NVARPT
         ELSE
            NSVAR  = NWOPPT
         END IF


      CALL QM3HESS(NOSIM,BOVECS,CREF,CMO,INDXCI,
     &             WRK(KDENS),SOVECS,NSVAR,WRK(KWRK),LWRK1)

      END IF
      CALL QEXIT('QM3LIN')
      RETURN
      END
C
C*************************************************************************
C  /* Deck qm3hess */
      SUBROUTINE QM3HESS(NOSIM,BOVECS,CREF,CMO,INDXCI,
     &                   DV,SVEC,NSVEC,WRK,LWRK)
C
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <maxash.h>
#include <infind.h>
#include <qm3.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
C
      DIMENSION BOVECS(NWOPPT,*),CREF(*),CMO(*),INDXCI(*)
      DIMENSION DV(*),SVEC(NSVEC,*),WRK(LWRK)
      LOGICAL   FULHES
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
C
      CALL QENTER('QM3HESS')
C
      KRXYO   = 1
      KUCMO   = KRXYO   + NOSIM*NNORBX
      KUBO    = KUCMO   + NORBT*NBAST
      KFSOLAO = KUBO    + NOSIM*N2ORBX
      KFSOLMO = KFSOLAO + NNBASX
      KUFSOL  = KFSOLMO + NNORBX
      KTGX    = KUFSOL  + N2ORBX
      KRRAOx  = KTGX    + NNORBX
      KRRAOy  = KRRAOx  + NNBASX
      KRRAOz  = KRRAOy  + NNBASX
      KTRAO   = KRRAOz  + NNBASX
      KTRMO   = KTRAO   + NNBASX
      KUTRX   = KTRMO   + NNORBX
      KUTR    = KUTRX   + N2ORBX
      KTRX    = KUTR    + N2ORBX
      KWRK    = KTRX    + NNORBX
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3HESS',-KWRK,LWRK)

      IF (JWOPSY .NE. LSYMPT .OR. NWOPPT .NE. NWOPT) THEN
         WRITE (LUPRI,*) 'QM3HESS ERROR: JWOPSY .ne. LSYMPT .or.'//
     *      ' NWOPPT .ne. NWOPT'
         WRITE (LUPRI,*) 'LSYMPT,JWOPSY:',LSYMPT,JWOPSY
         WRITE (LUPRI,*) 'NWOPPT,NWOPT :',NWOPPT,NWOPT
         CALL QUIT('QM3HESS ERROR,lsympt.ne.jwopsy or nwoppt.ne.nwopt')
      END IF

C     Determine if full Hessian or only orbital Hessian

      FULHES = (NSVEC .EQ. NVARPT)
      IF (FULHES) THEN
         JSOVEC = 1 + NCONST
      ELSE
         JSOVEC = 1
      END IF

      IF (IQM3PR .GT. 15) THEN
         WRITE (LUPRI,'(/A)') ' >>> SOLLNO - svec(orb) on entry'
         CALL OUTPUT(SVEC(JSOVEC,1),1,NWOPPT,1,NOSIM,
     *               NSVEC,NOSIM,1,LUPRI)
      END IF

C     1. Construct effective operators Tgxo and RRayo

      CALL DZERO(WRK(KRXYO),NOSIM*NNORBX )
      CALL DZERO(WRK(KUBO),NOSIM*N2ORBX)

C     1.a. Unpack symmetry blocked CMO
      CALL UPKCMO(CMO,WRK(KUCMO))
C     1.b. Calculate unpacked orbital trial vectors in UBO

      IF (NOSIM.GT.0) THEN
         DO 210 IOSIM = 1,NOSIM
            JUBO = KUBO + (IOSIM-1)*N2ORBX
            CALL UPKWOP(NWOPPT,JWOP,BOVECS(1,IOSIM),WRK(JUBO))
  210    CONTINUE
      END IF


      CALL DZERO(WRK(KFSOLMO),NNORBX)
      CALL DZERO(WRK(KUFSOL),N2ORBX)

C     1.c. Construct the Tg operator in ao and transform to mo
      CALL QM3FCKMO(CMO,WRK(KFSOLMO),WRK(KWRK),LWRK1,IPRINT)

      CALL DSPTSI(NORBT,WRK(KFSOLMO),WRK(KUFSOL))

      DO 200 IOSIM = 1,NOSIM
         JRXYO = KRXYO + (IOSIM-1)*NNORBX
         JUBO = KUBO + (IOSIM-1)*N2ORBX

C        1.d. Calculate one-index transformed Tg integrals

         CALL DZERO(WRK(KUTRX),N2ORBX)
         CALL TR1UH1(WRK(JUBO),WRK(KUFSOL),WRK(KUTRX),1)
         CALL DGETSP(NORBT,WRK(KUTRX),WRK(KTGX))

         FAC = 1.0D0
         CALL DAXPY(NNORBX,FAC,WRK(KTGX),1,WRK(JRXYO),1)


         WRITE (LUPRI,'(/A,I3,A,I3)')
     *      ' >>> SOLLNO - Rxo matrix no.',IOSIM,' of',NOSIM
         CALL OUTPAK(WRK(JRXYO),NORBT,1,LUPRI)


         CALL GPOPEN(LUQM3E,'ELFDMM','OLD',' ',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND(LUQM3E)

C        1.e. Readin Rra integrals and transform to mo

         LM = 0

         FAC = 1.0D0

         DO 520 I = 1, ISYTP
           IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
             DO 521 J = NSYSBG(I), NSYSED(I)
               DO 522 K = 1, NUALIS(I)
                 LM = LM + 1

                 FAC = -ALPIMM(I,K)

                 CALL DZERO(WRK(KRRAOx),NNBASX)
                 CALL READT(LUQM3E,NNBASX,WRK(KRRAOx))

C                Calculate one-index transformed Rra.x integrals

                 CALL DZERO(WRK(KTRMO),NNORBX)
                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)
                 CALL DZERO(WRK(KTRX),NNORBX)
                 CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                     NBAST,NORBT)
                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                 CALL TR1UH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX),1)
                 CALL DGETSP(NORBT,WRK(KUTRX),WRK(KTRX))

                 FACx = 0.0D0
                 DO 523 IW = 1,NISHT
                   IG = ISX(IW)
                   FACx = FACx + 2.0D0*WRK(KTRX -1 + IROW(IG+1))
  523            CONTINUE
                 FACx = FACx * FAC
  
                 CALL DAXPY(NNORBX,FACx,WRK(KTRMO),1,WRK(JRXYO),1)

                 CALL DZERO(WRK(KRRAOy),NNBASX)
                 CALL READT(LUQM3E,NNBASX,WRK(KRRAOy))

C                Calculate one-index transformed Rra.y integrals

                 CALL DZERO(WRK(KTRMO),NNORBX)
                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)
                 CALL DZERO(WRK(KTRX),NNORBX)
                 CALL UTHU(WRK(KRRAOy),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                     NBAST,NORBT)
                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                 CALL TR1UH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX),1)
                 CALL DGETSP(NORBT,WRK(KUTRX),WRK(KTRX))

                 FACy = 0.0D0
                 DO 524 IW = 1,NISHT
                   IG = ISX(IW)
                   FACy = FACy + 2.0D0*WRK(KTRX -1 + IROW(IG+1))
  524            CONTINUE
                 FACy = FACy * FAC

                 CALL DAXPY(NNORBX,FACy,WRK(KTRMO),1,WRK(JRXYO),1)

                 CALL DZERO(WRK(KRRAOz),NNBASX)
                 CALL READT(LUQM3E,NNBASX,WRK(KRRAOz))

C                Calculate one-index transformed Rra.z integrals

                 CALL DZERO(WRK(KTRMO),NNORBX)
                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)
                 CALL DZERO(WRK(KTRX),NNORBX)
                 CALL UTHU(WRK(KRRAOz),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                     NBAST,NORBT)
                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                 CALL TR1UH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX),1)
                 CALL DGETSP(NORBT,WRK(KUTRX),WRK(KTRX))

                 FACz = 0.0D0
                 DO 525 IW = 1,NISHT
                   IG = ISX(IW)
                   FACz = FACz + 2.0D0*WRK(KTRX -1 + IROW(IG+1))
  525            CONTINUE
                 FACz = FACz * FAC

                 CALL DAXPY(NNORBX,FACz,WRK(KTRMO),1,WRK(JRXYO),1)

  522          CONTINUE
  521        CONTINUE
           END IF
  520    CONTINUE
C
         CALL GPCLOSE(LUQM3E,'KEEP')

  200 CONTINUE

C      IF (IQM3PR .GT. 15) THEN
         DO 600 IOSIM = 1,NOSIM
            JRXYO = KRXYO + (IOSIM-1)*NNORBX
            WRITE (LUPRI,'(/A,I3,A,I3)')
     *         ' >>> SOLLNO - (Rxo + Ryo) matrix no.',IOSIM,' of',NOSIM
            CALL OUTPAK(WRK(JRXYO),NORBT,1,LUPRI)
  600    CONTINUE
C      END IF

C     2. Add effective operators Tgxo and RRayo contribution to SVEC(NSVEC,NOSIM) 
C     Tell SOLGO only to use the NWOPPT first JWOP entries
      DO 800 IOSIM = 1,NOSIM
         JRXYO  = KRXYO  + (IOSIM-1)*NNORBX
         CALL SOLGO(D2,DV,WRK(JRXYO),SVEC(JSOVEC,IOSIM))
  800 CONTINUE
C
      IF (IQM3PR .GT. 15) THEN
         WRITE (LUPRI,'(/A)') ' >>> SOLLNO - svec(orb) on exit'
         CALL OUTPUT(SVEC(JSOVEC,1),1,NWOPPT,1,NOSIM,
     *               NSVEC,NOSIM,1,LUPRI)
      END IF

      CALL QEXIT('QM3HESS')
      RETURN
      END
C*************************************************************************
C  /* Deck qm3fck */
      SUBROUTINE QM3FCKMO(CMO,FSOL,WRK,LWRK,IPRINT)
C
C     Construct the QM/MM contribution to the Fock-matrix in MO basis
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <inforb.h>
#include <infopt.h>
C
      DIMENSION CMO(*), FSOL(*), WRK(LWRK)
C
      CALL QENTER('QM3FCKMO')
C
      KDV     = 1
      KDENS   = KDV     + N2BASX
      KDVS    = KDENS   + NNBASX
      KFSOLAO = KDVS    + NNBASX
      KUCMO   = KFSOLAO + NNBASX
      KWRK    = KUCMO   + NORBT*NBAST
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3FCKMO',-KWRK,LWRK)

C     Construct the density matrix
      CALL FCKDEN((NISHT.GT.0),.FALSE.,WRK(KDV),
     *            DUMMY,CMO,DUMMY,WRK(KWRK),LWRK1)

      CALL DGEFSP(NBAST,WRK(KDV),WRK(KDVS))
      CALL PKSYM1(WRK(KDVS),WRK(KDENS),NBAS,NSYM,1)

C     Construct the QM/MM contribution to the Fock-matrix in AO
      CALL QM3FCK(WRK(KDENS),WRK(KDENS),WRK(KFSOLAO),ESOLT,
     *            ENSQM3,EPOQM3,EELEL,WRK(KWRK),LWRK1,IPRINT)

C     Transform to mo
      CALL UPKCMO(CMO,WRK(KUCMO))
      CALL UTHU(WRK(KFSOLAO),FSOL,WRK(KUCMO),WRK(KWRK),
     &             NBAST,NORBT)
C
      CALL QEXIT('QM3FCKMO')
      RETURN
      END
