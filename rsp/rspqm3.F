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
      SUBROUTINE QM3LNO(NOSIM,BOVECS,CREF,CMO,XINDX,UDV,DV,UDVTR,EVECS,
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
#include <mxcent.h>
#include <pcmdef.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <pcm.h>
#include <pcmlog.h>


      DIMENSION BOVECS(*), CMO(*), XINDX(*), UDV(*), EVECS(KZYVAR,*)
      DIMENSION UDVTR(*)
      DIMENSION WRK(*), DV(*), CREF(*)
      LOGICAL LOPEN


      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      LOGICAL FIRST1, FILE_EXIST
c      LOGICAL FIRST 
c      SAVE    FIRST
c      DATA    FIRST /.TRUE./

      IF (MMPCM) THEN

        KCOUNT = 0
        DO 100 I=1,NSTATES
          KCOUNT = KCOUNT + ICQM3(I)
  100   CONTINUE
        IF (KCOUNT.EQ.NSTATES) THEN
          FIRST1 = .TRUE.
        ELSE 
          FIRST1 = .FALSE.
        ENDIF

        IF (FIRST1) THEN
C         This is a new perturbation and we need to delete old interface
C         files between mm/pcm
C         1) The MM_TO_PCM file
          INQUIRE(FILE='MYFPCM',EXIST=FILE_EXIST)
          IF (FILE_EXIST) THEN
            LFILE = -1
            CALL GPOPEN(LFILE,'MYFPCM','OLD',' ','FORMATTED',IDUMMY,
     &                  .FALSE.)
            CALL GPCLOSE(LFILE,'DELETE')
          ENDIF 
C         2) do the same for the PCM_TO_MM file
          INQUIRE(FILE='QFQM3',EXIST=FILE_EXIST)
          IF (FILE_EXIST) THEN
            KFILE = -1
            CALL GPOPEN(KFILE,'QFQM3','OLD',' ','FORMATTED',IDUMMY,
     &                  .FALSE.)
            CALL GPCLOSE(KFILE,'DELETE')
          ENDIF
        ELSE 
C         open the file for PCM_TO_MM data
          KFILE = -1
          CALL GPOPEN(KFILE,'QFQM3','OLD',' ','FORMATTED',IDUMMY,
     &                .FALSE.)
               REWIND(KFILE)
        ENDIF
      ENDIF

C     See if some vectors have converged
c      DO 101 I=1,NOSIMFIRST
c        IF (ICQM3(I).EQ.-1) THEN
c          write(lupri,*) 'root converged:',I
c        ENDIF
c  101 CONTINUE
c      IF ( (NOSIM.NE.NOSIMOLD) .AND. (.NOT.FIRST1) )THEN
c        write(lupri,*) 'convergence since last iteration'
c      ENDIF

      IF (NASHT .GT. 0) THEN
         CALL QUIT('QM3 only implemented for closed shell')
      ENDIF

      LOPEN = .FALSE.

      CALL QENTER('QM3LNO')
      IF ( (.NOT. (OLDTG)) .AND. (LOSPC) ) THEN
         CALL QEXIT('QM3LNO')
         RETURN
      END IF

      KUCMO   = 1
      KUBO    = KUCMO   + NORBT*NBAST
      KFSOLMO = KUBO    + NOSIM*N2ORBX
      KUFSOL  = KFSOLMO + NNORBX
      KUTGX   = KUFSOL  + N2ORBX
      KRXYO   = KUTGX   + N2ORBX
      IF (TRPLET) THEN
        KRXYOT  = KRXYO   + NOSIM*N2ORBX
        KRRAOx  = KRXYOT  + NOSIM*N2ORBX
      ELSE
        KRRAOx  = KRXYO   + NOSIM*N2ORBX
      ENDIF
      KTRMO   = KRRAOx  + NNBASX
      KUTR    = KTRMO   + NNORBX
      KUTRX   = KUTR    + N2ORBX
      IF (MMPCM) THEN
        KXCOOR  = KUTRX   + N2ORBX
        KYCOOR  = KXCOOR  + NCOMS
        KZCOOR  = KYCOOR  + NCOMS
        KXDIP   = KZCOOR  + NCOMS
        KYDIP   = KXDIP   + NCOMS
        KZDIP   = KYDIP   + NCOMS
        KXPCM   = KZDIP   + NCOMS
        KYPCM   = KXPCM   + NTS*NOSIM
        KZPCM   = KYPCM   + NTS*NOSIM
        KQEX    = KZPCM   + NTS*NOSIM
        KSTATE  = KQEX    + NTS*NOSIM
        KWRK    = KSTATE  + NOSIM
        LWRK1   = LWRK    - KWRK 
      ELSE
        KWRK    = KUTRX   + N2ORBX
        LWRK1   = LWRK    - KWRK
      ENDIF
      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3LNO',-KWRK,LWRK)

      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KUBO),NOSIM*N2ORBX)
      CALL DZERO(WRK(KFSOLMO),NNORBX)
      CALL DZERO(WRK(KUFSOL),N2ORBX)
      CALL DZERO(WRK(KUTGX),N2ORBX)
      CALL DZERO(WRK(KRXYO),NOSIM*N2ORBX)
      IF (TRPLET) CALL DZERO(WRK(KRXYOT),NOSIM*N2ORBX)
      CALL DZERO(WRK(KRRAOx),NNBASX)
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      CALL DZERO(WRK(KUTRX),N2ORBX)
      IF (MMPCM) THEN
        CALL DZERO(WRK(KXCOOR),NCOMS)
        CALL DZERO(WRK(KYCOOR),NCOMS)
        CALL DZERO(WRK(KZCOOR),NCOMS)
        CALL DZERO(WRK(KXDIP),NCOMS)
        CALL DZERO(WRK(KYDIP),NCOMS)
        CALL DZERO(WRK(KZDIP),NCOMS)
        CALL DZERO(WRK(KXPCM),NTS*NOSIM)
        CALL DZERO(WRK(KYPCM),NTS*NOSIM)
        CALL DZERO(WRK(KZPCM),NTS*NOSIM)
        CALL DZERO(WRK(KQEX),NTS*NOSIM)
      ENDIF

      IF (IPRRSP .GE. 140) THEN
         WRITE (LUPRI,'(//A)') ' >>> TEST OUTPUT FROM QM3LNO >>>'
         WRITE (LUPRI,'(/A)') ' >>> QM3LNO - svec(orb,1) on entry'
         WRITE (LUPRI,'(/A)') ' Z - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Y - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZVAR+KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,
     *               1,1,LUPRI)
      END IF
      IF (MMPCM .AND. (.NOT.FIRST1)) THEN
C     Read in the PCM coordinates and charges but only for the states
C     that we are still converging (it might be that some of the
C     states on the PCM_TO_MM interface file have converged after the
C     last PCM correction and before this MM correction) 

         READ(KFILE,*) NOSIMPCM, NTSPCM
         IF (NOSIMPCM.NE.NOSIMOLD) THEN
            CALL QUIT('Error in no. of converged vectors in mm/pcm')
         ENDIF
         IF (NTS.NE.NTSPCM) THEN
            CALL QUIT('Error in no. of tesseras in mm/pcm')
         ENDIF
         KL = 0
         DO 150 IOSIM = 1,NOSIMOLD ! nosimold entries on this file
            READ(KFILE,*) NVECIT
            IF (ICQM3(NVECIT).NE.-1) THEN ! -1 means converged
C           we only enter here NOSIM times even though the loop is over
C           NOSIMOLD   
               WRK(KSTATE+KL) = ICQM3(NVECIT) ! save state number
               DO 151 L = 1,NTSPCM
                  READ(KFILE,'(4(E25.15,2x))') 
     &                 WRK(KXPCM+KL*NTS+L-1),
     &                 WRK(KYPCM+KL*NTS+L-1),
     &                 WRK(KZPCM+KL*NTS+L-1),
     &                 WRK(KQEX+KL*NTS+L-1)
 151           CONTINUE
               KL=KL+1
            ELSE                ! this vector has converged and we do not need this
               DO 152 L = 1,NTSPCM
                  READ(KFILE,'(4(E25.15,2x))')
     &                 WRK(KWRK),
     &                 WRK(KWRK),
     &                 WRK(KWRK),
     &                 WRK(KWRK)
 152           CONTINUE
C     do not update KL in this case
            ENDIF
 150     CONTINUE
      ENDIF

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

       IF (INTDIR) THEN
         OBKPX = DIPORG(1)
         OBKPY = DIPORG(2)
         OBKPZ = DIPORG(3)
       ENDIF

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

      IF (MMPCM) THEN
         INQUIRE(FILE='MYFPCM',EXIST=FILE_EXIST)
         IF (FILE_EXIST) THEN
            LFILE = -1
            CALL GPOPEN(LFILE,'MYFPCM','OLD',' ','FORMATTED',IDUMMY,
     &              .FALSE.)
            CALL GPCLOSE(LFILE,'DELETE')
         ENDIF
         LFILE = -1
         CALL GPOPEN(LFILE,'MYFPCM','NEW',' ','FORMATTED',IDUMMY,
     &              .FALSE.)
         REWIND(LFILE)
         WRITE(LFILE,*) NOSIM, NCOMS
      ENDIF

      DO 200 IOSIM = 1,NOSIM

         JRXYO = KRXYO + (IOSIM-1)*N2ORBX
         IF (TRPLET) JRXYOT = KRXYOT + (IOSIM-1)*N2ORBX
         JUBO  = KUBO  + (IOSIM-1)*N2ORBX

         CALL DZERO(WRK(KUTGX),N2ORBX)

         CALL ONEXH1(WRK(JUBO),WRK(KUFSOL),WRK(KUTGX))

         FAC = 1.0D0
         IF (TRPLET) THEN
           CALL DAXPY(N2ORBX,FAC,WRK(KUTGX),1,WRK(JRXYOT),1)
         ELSE
           CALL DAXPY(N2ORBX,FAC,WRK(KUTGX),1,WRK(JRXYO),1)
         ENDIF
 
C        Since only hf/dft closed shell no contribution can arrise from 
C        the one-index transformed electric field if triplet.
C        Also, if LGSPOL neglect the time-dependece of the expectation value 
C        contained in the interaction operator, i.e. keep the ground state polarization only. 
C        
  
         IF (TRPLET. OR. LGSPOL) GOTO 300

C-------------------------------------------------------
C 2. Construct -alpha \sum_a <0|Q_a|0>Rra
C-------------------------------------------------------

         IF (.NOT. INTDIR) THEN
           LUCCEF=-1
           CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &                 'UNFORMATTED',IDUMMY,.FALSE.)
           REWIND(LUCCEF)
         ENDIF
C
C        Readin Rra integrals (ao) and transform to mo
C        ---------------------------------------------
C
         LM = 0

         IF (INTDIR) L = NUSITE + NUALIS(0)

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

                 IF (INTDIR) THEN
                   KMAT = KWRK
                   KLAST = KMAT + 3*NNBASX
                   LWRK2 = LWRK - KLAST + 1
                   IATNOW = NUCIND + L + LM

                   KPATOM = 0
                   NOSIMI = 3
                   TOFILE = .FALSE.
                   TRIMAT = .TRUE.
                   EXP1VL = .FALSE.
                   DIPORG(1) = CORD(1,IATNOW)
                   DIPORG(2) = CORD(2,IATNOW)
                   DIPORG(3) = CORD(3,IATNOW)

                   FORQM3 = .TRUE.
                   CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIMI,WRK(KLAST),
     &                         LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                         KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
                   FORQM3 = .FALSE.

                   CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &                       NBAST,NORBT)
                   CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                   CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
                   FACx = FAC * QM3QLM(WRK(KUTRX))
                   IF ( (MMPCM) .AND. (.NOT.FIRST1) ) THEN
                     TEMP = 0.0D0
                     DO 523 IPCM = 1, NTS
                       XPCMNOW = WRK(KXPCM+(IOSIM-1)*NTS+IPCM-1)
                       YPCMNOW = WRK(KYPCM+(IOSIM-1)*NTS+IPCM-1)
                       ZPCMNOW = WRK(KZPCM+(IOSIM-1)*NTS+IPCM-1)
                       QPCMNOW = WRK(KQEX+(IOSIM-1)*NTS+IPCM-1)
                       DIST2 = (CORD(1,IATNOW)-XPCMNOW)**2
     &                       + (CORD(2,IATNOW)-YPCMNOW)**2
     &                       + (CORD(3,IATNOW)-ZPCMNOW)**2
                       DIST1 = SQRT(DIST2)
                       DISTM3 = 1.0/(DIST1**3)
                       TEMP = TEMP 
     &                      + QPCMNOW*(CORD(1,IATNOW)-XPCMNOW)*DISTM3
  523                CONTINUE
                     FACx = FACx + FAC*TEMP
                   ENDIF
                   CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,WRK(JRXYO),1)
                   IF (MMPCM) THEN
                      WRK(KXCOOR+LM-1) = CORD(1,IATNOW)
                      WRK(KYCOOR+LM-1) = CORD(2,IATNOW)
                      WRK(KZCOOR+LM-1) = CORD(3,IATNOW)
                      WRK(KXDIP+LM-1)  = -FACx
                   ENDIF 
                 ELSE
                   CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                   CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KWRK),NBAST,NORBT)
                   CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                   CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
                   FACx = FAC * QM3QLM(WRK(KUTRX))
                   CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,WRK(JRXYO),1)
                 ENDIF

C y-component
                 CALL DZERO(WRK(KRRAOx),NNBASX)
                 CALL DZERO(WRK(KTRMO),NNORBX)
                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)

                 IF (INTDIR) THEN
                   CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KLAST),NBAST,NORBT)
                 ELSE
                    CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                    CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),
     &                        WRK(KWRK),NBAST,NORBT)
                 ENDIF
                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                 CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
                 FACy = FAC * QM3QLM(WRK(KUTRX))
                 IF ( (MMPCM) .AND. (.NOT.FIRST1) ) THEN
                   TEMP = 0.0D0
                   DO 524 IPCM = 1, NTS
                     XPCMNOW = WRK(KXPCM+(IOSIM-1)*NTS+IPCM-1)
                     YPCMNOW = WRK(KYPCM+(IOSIM-1)*NTS+IPCM-1)
                     ZPCMNOW = WRK(KZPCM+(IOSIM-1)*NTS+IPCM-1)
                     QPCMNOW = WRK(KQEX+(IOSIM-1)*NTS+IPCM-1)
                     DIST2 = (CORD(1,IATNOW)-XPCMNOW)**2
     &                     + (CORD(2,IATNOW)-YPCMNOW)**2
     &                     + (CORD(3,IATNOW)-ZPCMNOW)**2
                     DIST1 = SQRT(DIST2)
                     DISTM3 = 1.0/(DIST1**3)
                     TEMP = TEMP 
     &                    + QPCMNOW*(CORD(2,IATNOW)-YPCMNOW)*DISTM3
  524              CONTINUE
                   FACy = FACy + FAC*TEMP
                 ENDIF
                 CALL DAXPY(N2ORBX,FACy,WRK(KUTR),1,WRK(JRXYO),1)
                 IF (MMPCM) WRK(KYDIP+LM-1)  = -FACy

C z-component
                 CALL DZERO(WRK(KRRAOx),NNBASX)
                 CALL DZERO(WRK(KTRMO),NNORBX)
                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)

                 IF (INTDIR) THEN
                   CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KLAST),NBAST,NORBT)
                 ELSE
                   CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                   CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KWRK),NBAST,NORBT)
                 ENDIF
                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                 CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
                 FACz = FAC * QM3QLM(WRK(KUTRX))
                 IF ( (MMPCM) .AND. (.NOT.FIRST1) ) THEN
                   TEMP = 0.0D0
                   DO 525 IPCM = 1, NTS
                     XPCMNOW = WRK(KXPCM+(IOSIM-1)*NTS+IPCM-1)
                     YPCMNOW = WRK(KYPCM+(IOSIM-1)*NTS+IPCM-1)
                     ZPCMNOW = WRK(KZPCM+(IOSIM-1)*NTS+IPCM-1)
                     QPCMNOW = WRK(KQEX+(IOSIM-1)*NTS+IPCM-1)
                     DIST2 = (CORD(1,IATNOW)-XPCMNOW)**2
     &                     + (CORD(2,IATNOW)-YPCMNOW)**2
     &                     + (CORD(3,IATNOW)-ZPCMNOW)**2
                     DIST1 = SQRT(DIST2)
                     DISTM3 = 1.0/(DIST1**3)
                     TEMP = TEMP
     &                    + QPCMNOW*(CORD(3,IATNOW)-ZPCMNOW)*DISTM3
  525              CONTINUE
                   FACz = FACz + FAC*TEMP
                 ENDIF
                 CALL DAXPY(N2ORBX,FACz,WRK(KUTR),1,WRK(JRXYO),1)
                 IF (MMPCM) WRK(KZDIP+LM-1)  = -FACz

  522       CONTINUE
  521     CONTINUE
        END IF
  520 CONTINUE
C
      IF (.NOT. INTDIR) CALL GPCLOSE(LUCCEF,'KEEP')

  300 CONTINUE 

      IF (MMPCM) THEN
         IF (FIRST1) THEN
            WRITE(LFILE,*) IOSIM
         ELSE
            WRITE(LFILE,*) WRK(KSTATE+IOSIM-1)
         ENDIF 

         DO 820 L = 1,NCOMS
            WRITE(LFILE,'(6(E25.15,2x))') WRK(KXCOOR+L-1),
     &           WRK(KYCOOR+L-1),
     &           WRK(KZCOOR+L-1),WRK(KXDIP+L-1),WRK(KYDIP+L-1),
     &           WRK(KZDIP+L-1)
 820     CONTINUE
      ENDIF

  200 CONTINUE

      IF ( (MMPCM) .AND. (.NOT.FIRST1) ) CALL GPCLOSE(KFILE,'KEEP')
      IF (MMPCM) CALL GPCLOSE(LFILE,'KEEP')

      IF (INTDIR) THEN
        DIPORG(1) = OBKPX
        DIPORG(2) = OBKPY
        DIPORG(3) = OBKPZ
      ENDIF

C--------------------------------------------
C 3.  Add effective operators to response
C--------------------------------------------


      IF (TRPLET) THEN
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYOT))
      ELSE
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYO))
      ENDIF

      IF (IPRRSP .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> QM3LNO - svec(orb,1) on exit'
         WRITE (LUPRI,'(/A)') ' Z - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Y - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZVAR+KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,
     *               1,1,LUPRI)
      END IF

C     Save nosim for use with mm/pcm

      IF (FIRST1) NOSIMFIRST = NOSIM
      NOSIMOLD = NOSIM

c      IF (FIRST) THEN
c         FIRST = .FALSE.
c      END IF
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
#include <mxcent.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <ccinftap.h>
#include "nuclei.h"

      PARAMETER ( D0=0.0D0, D1=1.0D0 )
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(LWRK),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2)
      DIMENSION MJWOP(2,MAXWOP,8)
      LOGICAL LCON, LORB, LREF, LOPEN
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QM3QRO')
      IF ( (.NOT. (OLDTG)) .AND. (LOSPC) ) THEN
         CALL QEXIT('QM3QRO')
         RETURN
      END IF

      LOPEN = .FALSE.
  
      IF (TRPLET) CALL QUIT('TRIPLET NOT IMPLEMENTED FOR QR-QM/MM YET')

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

C We assume no symmetry in the DFT/MM calculations although we have kept 
C some symmetry options below ...

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

C     If LGSPOL neglect the time-dependece of the expectation value 
C     contained in the interaction operator, i.e. keep the ground state polarization only. 

      IF (LGSPOL) GOTO 600


      IF (.NOT. INTDIR) THEN
        CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &                 'UNFORMATTED',IDUMMY,.FALSE.)
        REWIND(LUCCEF)
      ENDIF

C     Readin Rra integrals (ao) and transform to mo
C     ---------------------------------------------

      LM = 0

      IF (INTDIR) THEN
        L = NUSITE + NUALIS(0)
        OBKPX = DIPORG(1)
        OBKPY = DIPORG(2)
        OBKPZ = DIPORG(3)
      ENDIF 


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

              IF (INTDIR) THEN
                KMAT = KWRK
                KLAST = KMAT + 3*NNBASX
                LWRK2 = LWRK - KLAST + 1
                IATNOW = NUCIND + L + LM

                KPATOM = 0
                NOSIM = 3
                TOFILE = .FALSE.
                TRIMAT = .TRUE.
                EXP1VL = .FALSE.
                DIPORG(1) = CORD(1,IATNOW)
                DIPORG(2) = CORD(2,IATNOW)
                DIPORG(3) = CORD(3,IATNOW)

                FORQM3 = .TRUE.
                CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIM,WRK(KLAST),
     &                      LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                      KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
                FORQM3 = .FALSE.
                CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &                    NBAST,NORBT)
              ELSE
                CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                    NBAST,NORBT)
              ENDIF 

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

              IF (INTDIR) THEN
                CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                    WRK(KLAST),NBAST,NORBT)
              ELSE
                CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                    NBAST,NORBT)
              ENDIF

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

              IF (INTDIR) THEN
                CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                    WRK(KLAST),NBAST,NORBT)
              ELSE
                CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                    NBAST,NORBT)
              ENDIF 

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

      IF (INTDIR) THEN
        DIPORG(1) = OBKPX
        DIPORG(2) = OBKPY
        DIPORG(3) = OBKPZ
      ENDIF

      IF (.NOT. INTDIR) CALL GPCLOSE(LUCCEF,'KEEP')

  600 CONTINUE ! if lgspol

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
