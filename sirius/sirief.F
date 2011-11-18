C/* Deck iefcmmv */
      SUBROUTINE IEFCMMV(DMATM1,CPAK,SI,SPAK,DI,DPAK,
     $     WORK,IPVT,CORE,LWRK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <pi.h>
#include <infpri.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <symmet.h>
      PARAMETER (TPI = 2.0D0*PI, FPI = 4.0D0*PI)
C
      DIMENSION DMATM1(NTS,*),WORK(NTS),IPVT(NTS)
      DIMENSION SI(NTS,*),DI(NTS,*)
      DIMENSION CPAK(NTSIRR,*),DPAK(NTSIRR,*),SPAK(NTSIRR,*)
      DIMENSION DET(2)
C
C     CORE and CPAK point both to the same memory location. CORE is used
C     first as a scratch vector for the symmetrization e CPAK contains
C     the packed form of the A matrix (and A inverted after inversion).
C
      DIMENSION CORE(*)
#include <pcm.h>
#include <pcmlog.h>
C
      DATA ZERO,ONE,TWO/0.0D+00,1.0D+00,2.0D+00/
C
C        compute the C and then C-inverse matrix
C
ckr         WRITE(LUPRI,*) ' -------------------'
ckr         WRITE(LUPRI,*) ' -- MATRIX C^(-1) --'
ckr         WRITE(LUPRI,*) ' -------------------'
C
C   compute the C matrix, stored at -dmatm1-
C   isotropic dielectrics
C
      IF (NONEQ) THEN
         EPS0=EPS
         EPS=EPSINF
      ENDIF
      FACT=(EPS+1.D0)/(EPS-1.D0)
C      
      DO 140 I=1,NTS
C
C -- Calculation of the diagonal terms
C
        LI=ISPHE(I)
        SSI=1.07D0*DSQRT(AS(I)/FPI)
        DDI=-SSI*AS(I)/(2.0D0*RE(LI))
        SI(I,I)=SSI
        DI(I,I)=DDI
        XI=XTSCOR(I)
        YI=YTSCOR(I)
        ZI=ZTSCOR(I)
C
C -- Calculation of the off-diagonal terms
C
        DO 150 J=1,NTS
C
          IF (J.EQ.I) GOTO 150                  
          LJ=ISPHE(J) 
C
          XJ=XTSCOR(J)
          YJ=YTSCOR(J)
          ZJ=ZTSCOR(J)
          CXJ=(XJ-XE(LJ))/RE(LJ)
          CYJ=(YJ-YE(LJ))/RE(LJ)
          CZJ=(ZJ-ZE(LJ))/RE(LJ)
          DRIJ=DSQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
          DRIJ3=DRIJ**3
          SSI=AS(J)/(FPI*DRIJ)
          DDI=AS(I)*AS(J)
     *        *((XI-XJ)*CXJ+(YI-YJ)*CYJ+(ZI-ZJ)*CZJ)/(FPI*DRIJ3)
          SI(I,J)=SSI
          DI(I,J)=DDI
  150   CONTINUE
  140 CONTINUE
C
Clf Check if the off-diagonal terms of D and S are too big
C   If they are not, then are properly stored
C   if they are, then are set to zero
C
      IF(NPCMMT.GE.1) THEN
         DO I=1,NTS
            DO J=1,NTS
               IF(I.NE.J) THEN
                  SSI=SI(I,J)
                  DDI=DI(I,J)
                  AVGS=(ABS(SI(I,I))+ABS(SI(J,J)))/2
                  AVGD=(FACT-ABS(DI(I,I))-ABS(DI(J,J)))/2
                  IF (ABS(SSI).LE.AVGS) THEN
                     SI(I,J)=SSI
                  ELSE
                     WRITE(LUPRI,*)'Warning, element ',
     $                    I,J, 'of SI too big: set to zero'
                     SI(I,J)=ZERO
                  ENDIF
                  IF (ABS(DDI).LE.AVGD) THEN
                  ELSE
                     WRITE(LUPRI,*)'Warning, element ',
     $                    I,J, 'of DI too big: set to zero'
                     DI(I,J)=ZERO
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      END IF
C
C Symmetrization of SI and DI matrices
C
      CALL SYMPCM(DI,NTS,CORE,LWRK,NTSIRR)
      CALL PAKPCM(DI,NTS,DPAK,NTSIRR,MAXREP)
      CALL SYMPCM(SI,NTS,CORE,LWRK,NTSIRR)
      CALL PAKPCM(SI,NTS,SPAK,NTSIRR,MAXREP)
C
C Calculation of the C matrix in a symmetrical fashion
C
      DO LSYM = 0, MAXREP
         LSTART = NTSIRR * LSYM + 1
         CALL MATPCM(DPAK(1,LSTART),SPAK(1,LSTART),AS,CPAK(1,LSTART),
     $        NTSIRR,FACT,1)
      ENDDO               
C
Clf Check if the off-diagonal terms of DMATM1 are too big
C   If they are not, then are properly stored
C   if they are, then are set to zero
C     
      IF (NPCMMT.GE.2) THEN
         write(lupri,*) 'npcmmt grater than 1'
         DO LSYM =0, MAXREP
            LSTART = NTSIRR * LSYM
            DO JTS = 1,NTSIRR
               J = JTS + LSTART
               DO ITS=1,NTSIRR
                  I = ITS + LSTART
                  IF (ITS.NE.JTS) THEN
                     AVGD1=(ABS(CPAK(ITS,I))+ABS(CPAK(JTS,J)))/2
                     IF (ABS(CPAK(ITS,J)).GT.AVGD1) THEN
                        WRITE(LUPRI,*)'Warning, element ',
     $                       I,J, 'of DMATM1 too big: set to zero'
                        DMATM1(I,J)=ZERO
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END IF
C
      IF(IPRPCM.GE.5) THEN
         WRITE(LUPRI,*) ' IDMM: THE A MATRIX IS'
         CALL OUTPUT(DMATM1,1,NTS,1,NTS,NTS,NTS,1,LUPRI)
      END IF
      INFO=0
      DO ISYM = 0, MAXREP
         ISTART = NTSIRR * ISYM + 1
         CALL DGEFA(CPAK(1,ISTART),NTSIRR,NTSIRR,IPVT,INFO)
         IF(INFO.NE.0) THEN
            WRITE(LUPRI,*) ' THE A MATRIX IS SINGULAR'
            CALL QUIT('Singular A matrix in IEFCMMV')
         END IF
         CALL DGEDI(CPAK(1,ISTART),NTSIRR,NTSIRR,IPVT,DET,WORK,01)
      ENDDO
      IF(IPRPCM .GE. 2) THEN
         WRITE(LUPRI,*) ' IDMM: THE A INVERSE MATRIX IS'
         CALL OUTPUT(CPAK,1,NTSIRR,1,NTS,NTSIRR,NTS,1,LUPRI)
      END IF
C
Clf symmetry-blocked multiplication
C
      DO LSYM = 0, MAXREP
         LSTART = NTSIRR * LSYM + 1
         CALL MATPCM(DPAK(1,LSTART),SPAK(1,LSTART),AS,CPAK(1,LSTART),
     $        NTSIRR,FACT,2)
      ENDDO               
      IF(IPRPCM .GE. 5) THEN
         WRITE(LUPRI,*) ' IDMM: THE A INVERSE MATRIX IS'
         CALL OUTPUT(CPAK,1,NTSIRR,1,NTS,NTSIRR,NTS,1,LUPRI)
      END IF
      IF (NONEQ) EPS=EPS0
C
      RETURN
      END
C/* Deck icvev */
      SUBROUTINE ICVEV(DMATM1,SI,SPAK,BEMPOT,VEC,VPOT,VCAM,
     *                 NTSX,L2,SOME,WORK,LWORK)
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <codata.h>
      PARAMETER (FPI = 4.0D0*PI)
#include <infpri.h>
#include <nuclei.h>
#include <inforb.h>
#include <pcmnuclei.h>
#include <pcm.h>
#include <pcmlog.h>
#include <orgcom.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <symmet.h>
C
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION DMATM1(NTSIRR,*),SPAK(NTSIRR,*),
     *          SI(NTS,*),VEC(*),
     *          VPOT(*),VCAM(*),BEMPOT(L2)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WORK(LWORK)
C
      LOGICAL SOME
C
      PARAMETER (ANTOAU=1.0D+00/XTANG)
      LOGICAL TOFILE, TRIMAT, VCHECK
C
      DATA ZERO, ONE, TWO  /0.0D+00, 1.0D+00, 2.0D+00/
C
Clf Do not insert other definitions below the following include
C
#include <ibtfun.h>
ClfC
ClfC     We use as a quick way of transfering tessera coordinates to HERMIT
ClfC     the dipole origin. Need to be restored.
ClfC
Clf      XI = DIPORG(1)
Clf      YI = DIPORG(2)
Clf      ZI = DIPORG(3)
C
C  2) Calculation of apparent charges generated by the solute's nuclei.
C
      
      DO  ITS = 1, NTSIRR
Clf         DIPORG(1) = XTSCOR(ITS)
Clf         DIPORG(2) = YTSCOR(ITS)
Clf         DIPORG(3) = ZTSCOR(ITS)
Clf         NCOMP = 1
Clf         NTESP = 1
Clf         TOFILE = .TRUE.
Clf         TRIMAT = .TRUE.
Clf         KPATOM = 0
ClfC
ClfC        Calculates nuclear potential energy integrals (in AO basis) for
ClfC        the given tessera
ClfC
Clf         L=1
Clf         CALL GET1IN(DUMMY,'NPETES ',0,WORK,LWORK,LABINT,INTREP,
Clf     &                     INTADR,L,.TRUE.,KPATOM,.TRUE.,IPRPCM)
         L=ISPHE(ITS)
         CX=(XTSCOR(ITS)-XE(L))/RE(L)
         CY=(YTSCOR(ITS)-YE(L))/RE(L)
         CZ=(ZTSCOR(ITS)-ZE(L))/RE(L)
         VCAM(ITS) = ZERO
         VPOT(ITS) = ZERO
         DO  JATOM = 1, NUCDEP  
            XL=PCMCORD(1,JATOM)    
            YL=PCMCORD(2,JATOM)   
            ZL=PCMCORD(3,JATOM)  
            RIL=DSQRT((XTSCOR(ITS)-XL)**2+(YTSCOR(ITS)-YL)**2
     &                +(ZTSCOR(ITS)-ZL)**2) 
            XM=((XTSCOR(ITS)-XL)*CX+(YTSCOR(ITS)-YL)*CY
     &         +(ZTSCOR(ITS)-ZL)*CZ)/RIL**3 
            VPOT(ITS) = VPOT(ITS) - PCMCHG(JATOM)/RIL
            VCAM(ITS) = VCAM(ITS) + XM*PCMCHG(JATOM) 
         ENDDO
      ENDDO 
ClfC
ClfC     We now collect all potential integrals in one file (AOPROPER), sorted by
ClfC     symmetry transformation properties
ClfC
Clf      CALL PCMSRT(WORK,NNBASX)
ClfC
Clf      DIPORG(1) = XI
Clf      DIPORG(2) = YI
Clf      DIPORG(3) = ZI
      NOP = MAXREP + 1
      SQTNOP = DBLE(NOP)
C
C Potentials and field components are reproduced for the equivalent tesserae
C
      DO ISYM = 1, MAXREP
         ISTART = ISYM * NTSIRR
         DO ITS = 1, NTSIRR
            JTS = ITS + ISTART
            VPOT(JTS) = VPOT(ITS)
            VCAM(JTS) = VCAM(ITS)
         ENDDO
      ENDDO
C
Clf   
C     Symmetrization of nuclear potentials on tesserae: the actual
C     copmutational procedure is not the shortest possible. The shortest
C     is the computation of the potential for the tesserae from 1 to
C     NTSIRR and a scaling of them by SQTNOP to account for
C     normalization. This is because the totally symmetric representation
C     is the only one that has non zero elements. I want to keep the
C     current structure because it is formally cleaner. If a faster code
C     will be needed I will change it.
Clf
C
      CALL DZERO(WORK,2*NTS)
      DO ISYM = 0, MAXREP
         ISTART = ISYM * NTSIRR
         DO JSYM = 0, MAXREP
            JSTART = JSYM * NTSIRR
            UMATIJ = PT(IBTAND(ISYM,JSYM)) / SQTNOP 
            DO K = 1, NTSIRR
               ITS = ISTART + K
               JTS = JSTART + K
               WORK(ITS) = WORK(ITS) + UMATIJ * VPOT(JTS)
               WORK(ITS + NTS) = WORK(ITS + NTS) + UMATIJ * VCAM(JTS)
            ENDDO
         ENDDO
      ENDDO
      VCHECK = .FALSE.
      DO I = 1, NTSIRR
         VPOT(I) = WORK(I)
         VCAM(I) = WORK(I + NTS)
      ENDDO
      DO I = NTSIRR + 1, NTS
         IF (WORK(I) .GT. 1.0D-8) THEN
            WRITE(lupri,*) 'NUC POT TOO BIG', I, WORK(I)
            VCHECK = .TRUE.
         ENDIF
         IF (WORK(I+NTS) .GT. 1.0D-8) THEN
            WRITE(lupri,*) 'NUC FIELD TOO BIG', I, WORK(I + NTS)
            VCHECK = .TRUE.
         ENDIF
         VCAM(I)=ZERO
         VPOT(I)=ZERO
      ENDDO
      IF (VCHECK) CALL QUIT ('ICVEV: SYMMETRIZATION ERROR')
C
C     If excited state calculations, we need to include corrections
C     from the "slow" charges. We do this by creating a new effective
C     nuclear charge
C
Clf NON EQUILIBRIUM PART TO BE CHECKED FOR SYMMETRIZATION
C
      IF (NONEQ) THEN
         FACT=(EPS-EPSINF)/(EPS-1.0D0)
C
C Computation and symmetrization of SI
C
         DO I = 1, NTS
            XI = XTSCOR(I)
            YI = YTSCOR(I)
            ZI = ZTSCOR(I)
            SI(I,I) = 1.07D0 * DSQRT (FPI / AS(I))
            DO J = 1, I - 1
               XJ = XTSCOR(J)
               YJ = YTSCOR(J)
               ZJ = ZTSCOR(J)
               DRIJ = DSQRT((XI-XJ)**2 + (YI-YJ)**2 + (ZI-ZJ)**2)
               SI(I,J) = 1/DRIJ 
               SI(J,I) = SI(I,J)
            ENDDO
         ENDDO
         CALL SYMPCM(SI,NTS,WORK,LWORK,NTSIRR)
         CALL PAKPCM(SI,NTS,SPAK,NTSIRR,MAXREP)
C
Clf I suppose we store symmetrized charges (probably zero beyond NTSIRR)
C
         CALL PCMINR(WORK,WORK(NTS + 1))
Clf
         DO ISYM = 0, MAXREP
            ISTART = ISYM * NTSIRR + 1
C
C Perform the operation Vnuc = Vnuc + FACT * SI * Qinr
C
            CALL DGEMV('N',NTSIRR,NTSIRR,FACT,SI(1,ISTART),NTSIRR,WORK
     $           ,1,ONE,VPOT(ISTART),1)
         ENDDO
      END IF
C
      QNTOT = ZERO
      QNUC2 = ZERO
C
Clf Not symmetrized version of nuclear charges calculation
C
      DO ISYM = 0, MAXREP
         ISTART = NTSIRR * ISYM + 1
C
C Perform the operation Qnuc = - DMATM1 * Vnuc
C
            CALL DGEMV('N',NTSIRR,NTSIRR,-1.0D0,DMATM1(1,ISTART),NTSIRR,
     $        VPOT(ISTART),1,ZERO,QSN(ISTART),1)
      ENDDO
C
      DO ITS = 1, NTSIRR
C
C For some bizarre reason VPOT is computed with the wrong sign!!
C
         POTCAVNUC(ITS) = -(MAXREP+1)*VPOT(ITS)
         QSN(ITS) = QSN(ITS) * AS(ITS)
         QNTOT = QNTOT + QSN(ITS)
         QNUC2 = QNUC2 + VCAM(ITS)*AS(ITS)
         IF(NONEQ) THEN
            QSN(ITS)=QSN(ITS)-FACT*WORK(ITS)
         ENDIF
      ENDDO
      QNUC2 = QNUC2/FPI
      IF(SOME) WRITE(LUPRI,1200) QNUC2 * SQTNOP
C
C  4) Renormalization of `nuclear' apparent charges.
C     If ICOMPCM=1 the correction is proportional to the tessera's area.
C     If ICOMPCM=2 or ICOMPCM=3 the correction is equal for each tessera.
C
      CHG = ZERO
      DO JATOM = 1, NUCDEP
         CHG = CHG + PCMCHG(JATOM)
      ENDDO
C
      TCH = - CHG * (EPS - ONE) / EPS
      QNTOTN = ZERO
C
C     Option 1:
C
      FN = ONE
Clf dirty trick
Clf      ICOMPCM = 0
CLF
      IF(ICOMPCM.EQ.1) THEN
         SUPTOT = STOT * ANTOAU * ANTOAU
         DIFF = TCH - QNTOT
         DO ITS = 1, NTS
            QSN(ITS) = QSN(ITS) + DIFF * AS(ITS) / SUPTOT
            QNTOTN = QNTOTN + QSN(ITS)
         ENDDO
C
C     Option 2:
C
      ELSE IF(ICOMPCM.EQ.2) THEN
         FN = TCH / QNTOT
         DO ITS = 1, NTS
            QSN(ITS) = QSN(ITS) * FN
            QNTOTN = QNTOTN + QSN(ITS)
         ENDDO
      END IF
C
C     Total apparent charge before and after renormalization.
C
      IF(SOME) THEN
         IF(ICOMPCM.NE.0) WRITE(LUPRI,1000) SQTNOP * QNTOT, TCH, QNTOTN
         IF(ICOMPCM.EQ.0) WRITE(LUPRI,1100) SQTNOP * QNTOT, TCH
      END IF
      RETURN
C
 1000 FORMAT(' NUCLEAR APPARENT CHARGE',F10.5/
     *       ' THEORETICAL',F10.5,' RENORMALIZED',F10.5)
Clf 1100 FORMAT(1X,'NUCLEAR APPARENT CHARGE',F10.5,' THEORETICAL',
Clf     *          F10.5,' NOT RENORMALIZED',F10.5)
 1100 FORMAT(' NUCLEAR APPARENT CHARGE',F10.5,/
     *       ' THEORETICAL',F10.5,' NOT RENORMALIZED')
 1200 FORMAT(/' ESTIMATE OF NUCLEAR CHARGE',F15.5)
      END
C/* Deck pcminr */
      SUBROUTINE PCMINR(QSEGR,POTGR)
C
C     Reads in ground-state charges
C
#include <implicit.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <priunit.h>
#include <pcm.h>
#include <pcmlog.h>
      DIMENSION QSEGR(NTS), POTGR(NTS)
#include <inftap.h>
C
      CALL QENTER('PCMINR')
C
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,FNSIFC,'OLD',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND (LUSIFC)
      CALL MOLLAB('PCMINFO ',LUSIFC,lupri)
      READ (LUSIFC) MTS
      IF (MTS .NE. NTS) THEN
         WRITE (LUPRI,'(/A)') 
     &        ' Inconsistent number of tesseras on LUSIFC'
         WRITE (LUPRI,'(A,I6,/A,I6)') ' Tesseras on LUSIFC : ',MTS,
     &                                ' Tesseras in the run: ',NTS
         CALL QUIT('Inconsistent number of tesseras found on LUSIFC')
      END IF
C     Here we read the total surface charges (nuc + el) of the reference
C     (the ground state in most of the cases) state and (i suppose) the
C     molecular potential (V(\rho^{el}+\gamma^{nuc})) of the same state
C     There are some points which are not already clear at the moment
C
C 1) The individual charges are saved with the right sign
C 2) They are NOT scaled for any factor (apart from renormalization if requested)
C 3) In the old version they were saved as difference because one of the two
C    had wrong sign. Now we save them as a sum with the "-" sign simply
C    because in this way the program works well :)
      READ (LUSIFC) 
      READ (LUSIFC)
      READ (LUSIFC) QSEGR 
      READ (LUSIFC) POTGR
Clf write ground state charges and potentials on output
      IF (IPRPCM .GE. 20) THEN
         WRITE(LUPRI,*) 'QSEGR AND POTGR'
         DO I=1,NTS
            WRITE(LUPRI,*) I,QSEGR(I),POTGR(I)
         ENDDO
      END IF
Clf
C
      CALL QEXIT('PCMINR')
      RETURN
      END
C/* Deck sympcm */
      SUBROUTINE SYMPCM(DMAT,NTS,WORK,LWORK,NTSIRR)
C
C Authors: Luca Frediani and Kenneth Ruud
C
C Date: 2001 November 9th
C
C Purpose: Put the IEF-PCM matrices in a block-diagonal form
C 
C     This matrices are originally of dimension NTS * NTS where NTS is
C     the number of tesserae of the cavity. After the
C     block-diagonalization we get MAXREP+1 blocks of dimension
C     NTS/(MAXREP+1) each. This is done computing a matrix which is the
C     character table of the point group (ordered in the Dalton way
C     instead as it is found on books).  
C     The transformation is carried out with a unitary matrix whose elements
C     are the UIJ derived from the character table.
C     For the time being we leave it
C     unpacked in the same memory location of the input matrix. We
C     do not think it is convenient to pack it.
C    
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
#include <symmet.h>
      DIMENSION DMAT(NTS,*), WORK(NTS,*), UMAT(0:7,0:7)
      LOGICAL NOWSTOP
#include <ibtfun.h>
C
C      
      NOP = MAXREP + 1
      IF (LWORK .LT. NTS*NTS) THEN
         WRITE(LUPRI,300) LWORK,NTS*NTS
         CALL QUIT('Not enough memory to symmetrize PCM matrix')
      ENDIF
C
C Computation of UMAT
C
      CALL DZERO(WORK(1,1),NTS*NTS)
      DO ISYM = 0, MAXREP
         DO JSYM = 0, MAXREP
            UMAT(ISYM, JSYM) = ONE * PT(IBTAND(ISYM,JSYM))
         ENDDO
      ENDDO

      

C
C Computation of WORK = DMAT * UMAT
C      
ckr
ckr      ntsirr = nts
ckr      
      IF (NTSIRR * NOP .NE. NTS) THEN
         WRITE(LUPRI,200) NTS, NTSIRR, NOP
         CALL QUIT(
     $        'PCM: Wrong no. of tesserae for a symmmetry calculation')
      ENDIF
      DO ISYM = 0, MAXREP
         ISTART = NTSIRR * ISYM
         DO KSYM = 0, MAXREP
            KSTART = NTSIRR * KSYM
            DO JSYM = 0, MAXREP
               JSTART = NTSIRR * JSYM
C
C Instead of scaling twice by NOP ** 1/2 we scale only once by NOP
C
               UJK = UMAT(JSYM,KSYM) / NOP
               DO I=1, NTSIRR
                  ITS = ISTART + I
                  DO II = 1, NTSIRR
                     JTS = JSTART + II
                     KTS = KSTART + II
                     WORK(ITS,KTS) = WORK(ITS,KTS) + 
     $                    DMAT(ITS,JTS) * UJK
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C Computation of DMAT = UMAT * WORK
C      
      CALL DZERO(DMAT(1,1),NTS*NTS)
      DO ISYM = 0, MAXREP
         ISTART = NTSIRR * ISYM
         DO KSYM = 0, MAXREP
            KSTART = NTSIRR * KSYM
            DO JSYM = 0, MAXREP
               JSTART = NTSIRR * JSYM
               UIJ = UMAT(ISYM,JSYM)
               DO I=1, NTSIRR
                  ITS = ISTART + I
                  JTS = JSTART + I
                  DO II = 1, NTSIRR
                     KTS = KSTART + II
                     DMAT(ITS,KTS) = DMAT(ITS,KTS) + 
     $                               UIJ * WORK(JTS,KTS)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C         
C     Check that we have done the block-diagonalization in the proper way 
C (this test will be skipped after debugging)  
C
      NOWSTOP = .FALSE.
      DO ISYM = 0, MAXREP
         ISTART = NTSIRR * ISYM
         DO JSYM = 0, MAXREP
            JSTART = NTSIRR * JSYM
            IF (ISYM .NE. JSYM) THEN
               DO I=1, NTSIRR
                  ITS = ISTART + I
                  DO II = 1, NTSIRR
                     JTS = JSTART + II
                     IF (DMAT(ITS,JTS) .GT. 1.0D-10) THEN
                        WRITE(LUPRI,*) 'SIMMETRY BREAKING ELEMENT',
     $                       ITS,JTS,DMAT(ITS,JTS)
                        NOWSTOP=.TRUE.
                     ELSE
                        DMAT(ITS,JTS) = ZERO
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
ckr      IF (NOWSTOP) CALL QUIT(
ckr     $     'Not able to get a block-diagonal PCM matrix')
 100  FORMAT(8F7.3)
 200  FORMAT('Wrong no. of tesserae for a symmetry calculation',3I5)
 300  FORMAT('Not enough memory in sympcm. Free: ',I8,'Required: ',I8)
      RETURN
      END
C
C/* Deck pakpcm */
      SUBROUTINE PAKPCM(DUNPAK,NTS,DPAK,NTSIRR,MAXREP)
C
C Authors: Luca Frediani and Kenneth Ruud
C
C Date: 2001 November 11th
C
C Purpose: Put the IEF-PCM matrices in a packed form
C 
C     The input of this routine is the unpacked block-diagonal matrix of
C     the PCM calculation in case we have symmetry. The output is a
C     matrix of dimension NTSIRR * NTS (where NTSIRRR = NTS/(MAXREP + 1) )
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION DUNPAK(NTS,*), DPAK(NTSIRR,*)
      
      IF (NTSIRR*(MAXREP+1) .NE. NTS) THEN
         WRITE(LUPRI,200) NTS, NTSIRR, MAXREP
         CALL QUIT(
     $        'PCM: Wrong no. of tesserae for a symmmetry calculation')
      ENDIF      
C      CALL DZERO(DPAK(1,1),NTSIRR*NTS)
      DO ISYM=0, MAXREP
         ISTART = NTSIRR * ISYM  + 1
         CALL MCOPY(NTSIRR, NTSIRR, DUNPAK(ISTART,ISTART), NTS,
     $        DPAK(1,ISTART), NTSIRR)
C             SUBROUTINE MCOPY(NROWA,NCOLA,A,NRDIMA,B,NRDIMB)
      ENDDO
C
 200  FORMAT('Wrong no. of tesserae for a symmetry calculation',3I5)
      RETURN
      END
C
C/* Deck matpcm */
      SUBROUTINE MATPCM(D,S,A,C,N,FACT,IFLAG)
C
C Authors: Luca Frediani and Kenneth Ruud
C
C Date: 2001 November 11th
C
C Purpose: MARICES MANIPULATION IN THE IEF CALCULATION
C 
C     The input of this routine is one block of the packed
C     block-diagonal matrix of the PCM calculation (in case we have a
C     certain symmetry) and the complete matrix if symmetry is absent.
C     
C     Operation 1
C
C     C = ((FACT*A/2)-D)(A^(-1)*S)
C
C     FACT = (eps + 1) / (eps - 1) A = diagonal matrix of the areas of
C     the tesserae (passed as the diagonal vector 
C     D = matrix DI in iefcmmv 
C     S = matrix SI in iefcmmv
C
#include <implicit.h>
#include <pi.h>
C
      PARAMETER (FPI = 4.0D0 * PI, ZERO = 0.0D0, TWO = 2.0D0)
      DIMENSION D(N,*), S(N,*), A(N), C(N,*)
C
      IF (IFLAG .EQ. 1) THEN
         CALL DZERO(C(1,1),N*N)
         DO I = 1, N
            DO K = 1, N
               DO J = 1, N
                  DELTAIJ = ZERO
                  IF (I.EQ.J) DELTAIJ = A(I) * FACT / TWO
Clf*                  C(I,K) = C(I,K)
Clf*     $                 + FPI * (DELTAIJ - D(I,J)) * S(J,K) / A(J)
                  C(I,K) = C(I,K)
     $                 + FPI * (DELTAIJ - D(I,J)) * S(J,K)
               END DO
            END DO
         END DO
      ELSE IF (IFLAG .EQ. 2) THEN
         DO I=1,N
            DO K=1,N
               DMAT=0.0D0
               DO J=1,N
                  DELTAJK = ZERO
                  IF (K .EQ. J) DELTAJK = A(J) / TWO
                  DMAT = DMAT - C(I,J) * (DELTAJK - D(J,K))
               ENDDO
               S(I,K)=DMAT
            ENDDO
         ENDDO
         DO J=1,N
            DO I=1,N
               C(I,J)= S(I,J)
            ENDDO
         ENDDO
      ELSE
         CALL QUIT('MATPCM: FLAG VALUE NOT ALLOWED')
      END IF
      RETURN
      END
CC/* Deck V2EQMV1 */
C      SUBROUTINE V2EQMV1(V2,V1,AMAT,N)
C#include <implicit.h>
C      DIMENSION V1(*),V2(*),AMAT(N,*)
C      CALL DZERO (V2,N)
C      DO J = 1, N
C         DO I = 1, N
C            V2(I) = V2(I) + AMAT(I,J) * V1(J)
C         ENDDO
C      ENDDO
C      RETURN
C      END
C
C/* Deck J4INT */
C
      SUBROUTINE J4INT(WLOCFL,WORK,LWORK)
C
#include <implicit.h>
#include <pi.h>
      PARAMETER (D0=0.0D0,DP5=0.50D0,FPI=4.0D0*PI)
#include <dummy.h>
#include <priunit.h>
#include <inftap.h>
#include <mxcent.h>
#include <inforb.h>
#include <chrnos.h>
#include <pcmdef.h>
      CHARACTER*8 RTNLBL(2), LFLABEL, LABINT(9*MXCENT)
      LOGICAL FNDLAB, EXP1VL, TOFILE, TRIMAT, FIRST
#include <pcm.h>
#include <pcmlog.h>
#include "orgcom.h"
#include <maxorb.h>
#include <maxaqn.h>
#include <symmet.h>
C
      DIMENSION WLOCFL(*), WORK(*), INTREP(9*MXCENT), INTADR(9*MXCENT)
      LPCMBK = 0
      LPRPBK = 0
      LPCMBK = LUPCMD
      IF (LUPCMD .LT. 0) CALL GPOPEN(LUPCMD,'PCMDATA','UNKNOWN',
     &     'SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND (LUPCMD)
      KFREE = 1
      KFRSAV = KFREE
      LFRSAV = LWORK
      NTNIRR = NTS * NTSIRR
      IF(.NOT.OUTFLD) 
     $     CALL MEMGET('REAL',LDMATM,NTNIRR,WORK,KFREE,LWORK)
      CALL MEMGET('REAL',LDNVEC,NTSIRR,WORK,KFREE,LWORK)
      CALL MEMGET('REAL',KQFLD,NTSIRR,WORK,KFREE,LWORK)
      IF(.NOT.OUTFLD) THEN
         IF (FNDLAB('D-LOCINV',LUPCMD)) THEN
            CALL READT(LUPCMD,NTNIRR,WORK(LDMATM))
         ELSE
            WRITE (LUPRI,'(/A)') ' Matrix label D-LOCINV not '//
     &           'found on file PCMDATA'
            CALL QUIT('Integral label not found in J4INT')
         END IF
      END IF
      IF(IPRPCM .GE. 10) THEN
         WRITE(LUPRI,*) '========== LF MATRIX INV PACK =========='
         CALL OUTPUT(WORK(LDMATM),1,NTSIRR,1,NTS,NTSIRR,NTS,1,LUPRI)
      END IF
C     
      DO IC = 1,3
         ISYM = ISYMAX(IC,1)+1

         IF (LUPROP .LT. 0) then
            CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',
     &           'SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
         END IF
         DO ITS=1,NTSIRR
            L = ISPHE(ITS)
            XP = XTSCOR(ITS)
            YP = YTSCOR(ITS)
            ZP = ZTSCOR(ITS)
            IF (IC.EQ.1) THEN
               WORK(LDNVEC+ITS-1) = (XE(L) - XP)/RE(L)
            ELSE IF (IC.EQ.2) THEN
               WORK(LDNVEC+ITS-1) = (YE(L) - YP)/RE(L)
            ELSE IF (IC.EQ.3) THEN
               WORK(LDNVEC+ITS-1) = (ZE(L) - ZP)/RE(L)
            END IF
         ENDDO
         
         LVEC2 = 1 + NNBASX * (IC - 1)
         CALL DZERO(WLOCFL(LVEC2),NNBASX)
C     
         REWIND (LUPROP)
         XI = DIPORG(1)
         YI = DIPORG(2)
         ZI = DIPORG(3)
         KDMATM = LDMATM + ISYMAX(IC,1) * NTSIRR * NTSIRR
         DO ITS = 1, NTSIRR
            IF(.NOT.OUTFLD) THEN
               DQ = D0
               DO JTS = 1,NTSIRR
                  INDEX = (JTS-1) * NTSIRR + ITS - 1 + KDMATM
                  DQ = DQ + WORK(INDEX)*WORK(LDNVEC+JTS-1)
               END DO
               WORK(KQFLD + ITS - 1) = DQ*AS(ITS)
Clf
               QLOC(ITS,IC) = WORK(KQFLD + ITS - 1)
            ELSE
               WORK(KQFLD + ITS - 1) = WORK(LDNVEC+ITS-1)*AS(ITS)
Clf
               QLOC(ITS,IC) = WORK(KQFLD + ITS - 1)
            END IF
         ENDDO
         CALL J1INT(WORK(KQFLD),.FALSE.,WLOCFL(LVEC2),1,.FALSE.,
     $        'NPETES ',ISYM,WORK(KFREE),LWORK)

      END DO
      DIPORG(1) = XI
      DIPORG(2) = YI
      DIPORG(3) = ZI
      CALL MEMREL('J4INT',WORK,KFRSAV,KFRSAV,KFREE,LWORK)
      IF (LPCMBK .LT. 0) CALL GPCLOSE(LUPCMD,'KEEP')
      RETURN
      END
C  /* Deck diplfn */
      SUBROUTINE DIPLFN(IPRINT)
C
C     Calculates nuclear part of the local-field induced
C     dipole moment (based on DIPNUC)
C
C     LF feb 2005
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0,D100=100.0D0)
#include <orgcom.h>
#include <abainf.h>
#include <nuclei.h>
#include <dipole.h>
#include <symmet.h>
#include <dorps.h>
#include <pcmdef.h>
#include <pcm.h>
#include <ibtfun.h>
C
C     Dipole moment
C
      CALL DZERO(DLFN,3)
      DO ICOOR = 1, 3
         DO ITS = 1, NTSIRR
            VNUC = D0
            XI = XTSCOR(ITS)
            YI = YTSCOR(ITS)
            ZI = ZTSCOR(ITS)
            DO IATOM = 1, NUCIND
               CHA = CHARGE(IATOM)
               IF (ABS(CHA) .LT. D100) THEN
                  FAC = FMULT(ISTBNU(IATOM))*CHA
                  IF (ISYMAX(ICOOR,1) .EQ. 0) THEN
                     XN=CORD(1,IATOM)
                     YN=CORD(2,IATOM)
                     ZN=CORD(3,IATOM)
                     DIST=DSQRT((XI-XN)**2+(YI-YN)**2+(ZI-ZN)**2)
                     VNUC = VNUC + FAC/DIST
                  END IF
               END IF
            END DO
            DLFN(ICOOR) = DLFN(ICOOR) + QLOC(ITS,ICOOR) * VNUC
         END DO
      END DO
      RETURN
      END
