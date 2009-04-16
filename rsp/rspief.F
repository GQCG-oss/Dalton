C  /* Deck pcmltr */
      SUBROUTINE PCMLTR(NCSIM,NOSIM,BCVECS,BOVECS,CREF,CMO,INDXCI,
     &                  UDV,DV,UDVTR,DVTR,DTV,DTVTR,SCVECS,SOVECS,
     &                  WRK,LWRK)
C
C LF BM RC 24 oct 01 based on PCMLIN
C
C Common driver for PCMLNC and PCMLNO
C
#include <implicit.h>
#include <dummy.h>

      DIMENSION BCVECS(*),BOVECS(*),CREF(*),CMO(*),INDXCI(*)
      DIMENSION UDV(*),DV(*),DTV(*),SCVECS(*),SOVECS(*),WRK(LWRK)
      DIMENSION UDVTR(*),DVTR(*),DTVTR(*)
C
C
#include <maxorb.h>
#include <priunit.h>
#include <infpri.h>
#include <infrsp.h>
#include <inftap.h>
#include <wrkrsp.h>
C
      CALL QENTER('PCMLTR')
C
      CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ',
     *        'UNFORMATTED',IDUMMY,.FALSE.)
      IF (NCSIM .GT. 0) THEN
Clara         IF (IPRRSP.GT.101) THEN
         IF (.true.) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFIRMED CONFIGURATION VECTOR'
            WRITE(LUPRI,*)' **** BEFORE IEFLNC **** '
            CALL OUTPUT(SCVECS,1,KZYVAR,1,NCSIM,KZYVAR,NCSIM,1,LUPRI)
         END IF
C         write (lupri,*) 'lara bef. IEFLNC'
         CALL IEFLNC(NCSIM,BCVECS,CREF,CMO,INDXCI,
     *               UDV,DV,UDVTR,DVTR,DTV,DTVTR,SCVECS,WRK,LWRK)
Clara         IF (IPRRSP.GT.101) THEN
         IF (.true.) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFIRMED CONFIGURATION VECTOR'
            WRITE(LUPRI,*)' **** AFTER IEFLNC **** '
            CALL OUTPUT(SCVECS,1,KZYVAR,1,NCSIM,KZYVAR,NCSIM,1,LUPRI)
         END IF
      END IF
      IF ( NOSIM .GT.0 ) THEN
         IF (IPRRSP.GT.101) THEN
C         IF (.true.) THEN
             WRITE(LUPRI,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
             WRITE(LUPRI,*)' **** BEFORE IEFLNO **** '
             CALL OUTPUT(SOVECS,1,KZYVAR,1,NOSIM,KZYVAR,NOSIM,1,LUPRI)
          END IF
C         PRINT *, 'lara bef. IEFLNO'
         CALL IEFLNO(NCSIM,NOSIM,BOVECS,CREF,CMO,INDXCI,
     *               UDV,DV,UDVTR,DVTR,SOVECS,WRK,LWRK)
C         PRINT *, 'lara aft. IEFLNO'
         IF (IPRRSP.GT.101) THEN
C         IF (.true.) THEN
             WRITE(LUPRI,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
             WRITE(LUPRI,*)' **** AFTER IEFLNO **** '
             CALL OUTPUT(SOVECS,1,KZYVAR,1,NOSIM,KZYVAR,NOSIM,1,LUPRI)
          END IF
      END IF
      IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      CALL QEXIT('PCMLTR')
      RETURN
      END
C  /* Deck ieflnc */
      SUBROUTINE IEFLNC(NCSIM,BCVEC,CREF,CMO,INDXCI,
     *                  UDV,DV,UDVTR,DVTR,DTV,DTVTR,SVEC,
     *                  WRK,LFREE)
C
C 24 OCT 01
C
C  Purpose:  Calculate MCSCF E2  contribution from a
C            pcm-ief surrounding medium to a csf trial vector.
C
#include <implicit.h>
      DIMENSION BCVEC(*),  CREF(*), CMO(*)
      DIMENSION INDXCI(*), UDV(*), DV(*),   DTV(N2ASHX,*)
      DIMENSION SVEC(KZYVAR,*),     WRK(*)
      DIMENSION UDVTR(*),DVTR(*),DTVTR(N2ASHX,*)
C
#include <iratdef.h>
C
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0 , THRZER = 1.0D-14 )
      LOGICAL FNDLAB
C
C
C  Used from common blocks:
C    INFINP : NLMSOL, LSOLMX, INERSI
C    INFORB : NNASHX, NNORBX, NNBASX, etc.
C    INFTAP : LUSOL, LBSYMB
C    INFPRI : LUERR
C
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <priunit.h>
#include <orgcom.h>
#include <infinp.h>
#include <inforb.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <inftap.h>
#include <infpri.h>
#include <pcm.h>
#include <pcmlog.h>
#include "sdpre.h"

C
Clf this is wrong!!! Must be fixed!!!!!11
      DIMENSION VTEX(NTS,NCSIM)
C
      CALL QENTER('IEFLNC')
C
C     Core allocation
C
      KUCMO  = 1
Clara define KJPXAO 
      KJPXAO = KUCMO  + NORBT*NBAST
      KJPX   = KJPXAO + NNBASX
Clara
C      KJPX   = KUCMO  + NORBT*NBAST
      KJPXAC = KJPX   + NNORBX
      KPCMXB = KJPXAC + NNASHX
      KXBAC  = KPCMXB + NCSIM*NNORBX
      KEXBAC = KXBAC  + NCSIM*NNASHX
      KW10   = KEXBAC + NCSIM
C
      KCHRG   = KW10
      KW20    = KCHRG + NTS
C     2.1 read rlmao in ao basis and transform to rlm in mo basis
      KDV     = KW20
Clara add definition of  KDW 
      KDW     = KDV  + NCSIM*NNBASX 
      KW30    = KDW  + NCSIM*N2BASX
C
      KJPCMAO = KW10
      KW50    = KJPCMAO + NNBASX
      LW30    = LFREE   + 1 - KW30
C
C     4.0 rspsor
      KURXC  = KW10   + NNORBX
      KURYC  = KURXC  + N2ORBX
      KW40   = KURYC  + N2ORBX
      KNEED  = MAX(KW20,KW30)
      KNEED  = MAX(KNEED,KW40)
      IF (KNEED .GT. LFREE) CALL ERRWRK('IEFLNC',KNEED,LFREE)
C
C
C
      CALL DZERO(WRK(KJPX),NNORBX)
Clara clear space for charges, jpxao
      CALL DZERO(WRK(KCHRG),NTS)
      CALL DZERO(WRK(KJPXAO),NNBASX)
C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
CLF Read Potential integrals on tesserae
C
      CALL DCOPY(NTS,QSN,1,WRK(KCHRG),1)
      CALL DAXPY(NTS,1.0D0,QSE,1,WRK(KCHRG),1)
Clara multiply by -1 the tot charges
      CALL DSCAL(NTS,-1.0D0,WRK(KCHRG),1)
C         write (lupri,*),'qn'
C         call output(qsn,1,nts,1,1,nts,1,1,lupri)
C         write (lupri,*),'qe'
C         call output(qse,1,nts,1,1,nts,1,1,lupri)
      CALL J1INT(WRK(KCHRG),.FALSE.,WRK(KJPXAO),1,.FALSE.,'NPETES ',
     &           1,WRK(KW30),LW30)
      CALL UTHU(WRK(KJPXAO),WRK(KJPX),WRK(KUCMO),WRK(KW30),NBAST,NORBT)
C
      IF (NASHT .GT. 0) THEN
         CALL GETAC2(WRK(KJPX),WRK(KJPXAC))
      END IF
      IF (IPRRSP .GE. 15) THEN
         WRITE (LUPRI,'(/A)') ' J+X_mo matrix:'
         CALL OUTPAK(WRK(KJPX),  NORBT,1,LUPRI)
         IF (NASHT .GT. 0) THEN
            WRITE (LUPRI,'(/A)') ' J+X_ac matrix:'
            CALL OUTPAK(WRK(KJPXAC),NASHT,1,LUPRI)
         END IF
      END IF
C
C     Expectation value of J+X
C
      EXPJPX = SOLELM(DV,WRK(KJPXAC),WRK(KJPX),ACJPX)
      IF (IPRRSP .GE. 6) THEN
         WRITE (LUPRI,'(A,F17.8)')
     *        ' --- J+X expectation value :',EXPJPX
         WRITE (LUPRI,'(A,F17.8)')
     *        ' --- active part of J1X    :',ACJPX
      END IF
C
      DO ICSIM = 1, NCSIM
Clara first clear for KDV, KDW and charges
         CALL DZERO(WRK(KCHRG),NTS)
         CALL DZERO(WRK(KDW),N2BASX)
         CALL DZERO(WRK(KDV),NNBASX)
Clara
#ifdef PCM_DEBUG
      write (lupri,*) 'lara dtv='
      call output(dtv(1,icsim),1,nasht,1,nasht,nasht,nasht,1,lupri)
#endif
clara
Clara change call to fckden with fckden2 and KDV with KDW
         CALL FCKDEN2(.FALSE.,.TRUE.,DUMMY,WRK(KDW),CMO,DTV(1,ICSIM),
     &               WRK(KW30),LW30)
Clara add call to dgfsp
Clara
#ifdef PCM_DEBUG
      write (lupri,*) 'lara kdw='
      call output(wrk(kdw),1,nbast,1,nbast,nbast,nbast,1,lupri)
#endif
clara
         CALL DGEFSP(NBAST,WRK(KDW),WRK(KDV))
clara 
#ifdef PCM_DEBUG
      write (lupri,*) 'lara kdv='
      call outpak (wrk(kdv),nbast,1,lupri)
#endif
clara 
C      if (abs(wrk(kdv+ 3)-0.00000000d0)<1d-7.and.
C     &    abs(wrk(kdv+10)+0.12249707d0)<1d-7) wrk(kdv:kdv+21) =
C     &(/ 0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 1.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0/)
C      if (abs(wrk(kdv+ 3)-0.12249707d0)<1d-7.and.
C     &    abs(wrk(kdv+10)+0.12249707d0)<1d-7) wrk(kdv:kdv+21) =
C     &(/ 0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0,
C     &   1.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 1.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0,
C     &   0.00000000000000d0, 0.00000000000000d0, 0.00000000000000d0/)
clara 
#ifdef PCM_DEBUG
      write (lupri,*) 'lara kdv after='
      call outpak (wrk(kdv),nbast,1,lupri)
#endif
clara 
         CALL J1INT(WRK(KCHRG),.TRUE.,WRK(KDV),1,.FALSE.,'NPETES ',
     &              KSYMOP,WRK(KW30),LW30)
         CALL DSCAL(NTS,-1.0D0,WRK(KCHRG),1)
Clara
#ifdef PCM_DEBUG
      write (lupri,*) 'lara kchrg='
      call output (wrk(kchrg),1,nts,1,1,nts,1,1,lupri)
#endif
clara
         CALL V2Q(WRK(KW30),WRK(KCHRG),VTEX(1,ICSIM),QET,NEQRSP)
         IF (KSYMOP .NE. 1) THEN
            CALL DCOPY(NTSIRR,VTEX((KSYMOP - 1)*NTSIRR + 1,ICSIM),1,
     &                        VTEX(1,ICSIM),1)
            CALL DZERO(VTEX(NTSIRR+1,ICSIM),NTS-NTSIRR)
         END IF
      END DO
      CALL GPCLOSE(LUPCMD,'KEEP')
C
      CALL DZERO(WRK(KPCMXB),NCSIM*NNORBX)
C
      DO ICSIM = 1, NCSIM
         CALL DZERO(WRK(KJPCMAO),NNBASX)
Clara
#ifdef PCM_DEBUG
      write (lupri,*) 'lara vtex='
      call output(vtex,1,nts,1,ncsim,nts,ncsim,1,lupri)
C      write (lupri,'(<1+(nts-1)/4>(4f18.12/))') vtex(:,icsim)
#endif
clara
         CALL J1INT(VTEX(1,ICSIM),.FALSE.,WRK(KJPCMAO),1,.FALSE.,
     &              'NPETES ',KSYMOP,WRK(KW30),LW30)
         JPCMXB = KPCMXB + (ICSIM - 1)*NNORBX
Clara
#ifdef PCM_DEBUG
      write (lupri,*) 'lara icsim=', icsim,' nnbasx=', nnbasx, 'jpcmao='
      call outpak(wrk(kjpcmao),nbast,1,lupri)
#endif
clara
         CALL UTHU(WRK(KJPCMAO),WRK(JPCMXB),WRK(KUCMO),WRK(KW30),
     &             NBAST,NORBT)
      END DO
C
Clara clear for kxbac
      CALL DZERO(WRK(KXBAC),NNASHX*NCSIM)
      DO ICSIM = 1, NCSIM
         JPCMXB = KPCMXB + (ICSIM - 1)*NNORBX
         JXBAC  = KXBAC  + (ICSIM - 1)*NNASHX
         CALL GETAC2(WRK(JPCMXB),WRK(JXBAC))
Clara
#ifdef PCM_DEBUG
      write (lupri,*) 'lara icsim=',icsim,' nnashx=',nnashx, 'jxbac='
      call outpak(wrk(jxbac),nbast,1,lupri)
#endif
clara
         IF (TRPLET) THEN
            TEMP = SOLELM(DVTR,WRK(JXBAC),WRK(JPCMXB),XBAC)
            TEMP = XBAC
         ELSE
            TEMP = SOLELM(DV,WRK(JXBAC),WRK(JPCMXB),XBAC)
         ENDIF
         WRK(KEXBAC - 1 + ICSIM) = XBAC
      END DO
C
      IF (.true.) THEN 
clara      IF (IPRRSP.GT.101) THEN
         WRITE(LUPRI,*)' LINEAR TRANSFORMED CONFIGURATION VECTOR'
         WRITE(LUPRI,*)' **** BEFORE SLVSC in IEFLNC **** '
         CALL OUTPUT(SVEC,1,KZYVAR,1,NCSIM,KZYVAR,NCSIM,1,LUPRI)
      END IF
C
      CALL SLVSC(NCSIM,0,NNASHX,BCVEC,CREF,SVEC,WRK(KXBAC),WRK(KJPXAC),
     *           WRK(KEXBAC),ACJPX,INDXCI,WRK(KW30),LW30)
C     CALL SLVSC(NCSIM,NOSIM,NNASHX,BCVECS,CREF,SVECS,
C    *           RXAC,RYAC,TRXAC,TRYAC,INDXCI,WRK,LWRK)
C
      IF (.true.) THEN 
Clara      IF (IPRRSP.GT.101) THEN
         WRITE(LUPRI,*)' LINEAR TRANSFORMED CONFIGURATION VECTOR'
         WRITE(LUPRI,*)' **** AFTER SLVSC in IEFLNC **** '
         CALL OUTPUT(SVEC,1,KZYVAR,1,NCSIM,KZYVAR,NCSIM,1,LUPRI)
      END IF
C
C     ... orbital part of sigma vector(s)
C
      IF (KZWOPT .GT. 0) THEN
         JPCMXB   = KPCMXB
C
         CALL DSPTSI(NORBT,WRK(KJPX),WRK(KURYC))
         DO 800 ICSIM = 1,NCSIM
C
            CALL DSPTSI(NORBT,WRK(JPCMXB),WRK(KURXC))
            IF (TRPLET) THEN
               CALL SLVSOR(.TRUE.,.FALSE.,1,UDVTR,
     *                     SVEC(1,ICSIM),WRK(KURXC))
            ELSE
               CALL SLVSOR(.TRUE.,.TRUE.,1,UDV,
     *                     SVEC(1,ICSIM),WRK(KURXC))
            END IF
            IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)' **** AFTER SLVSOR  in IEFLNC **** '
               WRITE(LUPRI,*)
     *         ' orbital part of LINEAR TRANSFORMED CONF VEC No',ICSIM
               CALL OUTPUT(SVEC(1,ICSIM),1,KZYVAR,1,1,KZYVAR,1,1,LUPRI)
            END IF
            IF (TRPLET) THEN
               CALL SLVSOR(.FALSE.,.FALSE.,1,DTVTR(1,ICSIM),
     *                      SVEC(1,ICSIM),WRK(KURYC))
            ELSE
               CALL SLVSOR(.FALSE.,.FALSE.,1,DTV(1,ICSIM),
     *                     SVEC(1,ICSIM),WRK(KURYC))
            END IF
            IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)
     *         ' orbital part of LINEAR TRANSFORMED CONF VEC No',ICSIM
               CALL OUTPUT(SVEC(1,ICSIM),1,KZYVAR,1,1,KZYVAR,1,1,LUPRI)
            END IF

            JPCMXB   = JPCMXB   + NNORBX
  800    CONTINUE
         IF (IPRRSP.GT.101) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFIRMED CONFIGURATION VECTOR'
            WRITE(LUPRI,*)' **** AFTER SLVSOR  in IEFLNC **** '
            CALL OUTPUT(SVEC,1,KZYVAR,1,NCSIM,KZYVAR,NCSIM,1,LUPRI)
         END IF
      END IF
C
      CALL QEXIT('IEFLNC')
      RETURN
C     end of IEFLNC.
      END
C  /* Deck ieflno */
      SUBROUTINE IEFLNO(NCSIM,NOSIM,BOVECS,CREF,CMO,INDXCI,
     *                  UDV,DV,UDVTR,DVTR,SVEC,
     *                  WRK,LFREE)
C
C
C  Purpose:  Calculate MCSCF E2 contribution from a
C            surrounding ief-pcm medium to an orbital trial vector.
C
C
#include <implicit.h>
      DIMENSION BOVECS(*), CREF(*), CMO(*)
      DIMENSION INDXCI(*),   UDV(*),   DV(*)
      DIMENSION SVEC(KZYVAR,*),    WRK(*)
      DIMENSION UDVTR(*),DVTR(*)
C
#include <iratdef.h>
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0, D2 = 2.0D0, DP5 = 0.50D0)
      LOGICAL FNDLAB, TOFILE, EXP1VL, TRIMAT
#include <dummy.h>
C
C
C  Used from common blocks:
C    INFINP : NLMSOL, LSOLMX, INERSI
C    INFORB : NNASHX, NNORBX, NNBASX, etc.
C    INFVAR : JWOP
C    INFRSP :
C    WRKRSP :
C    INFTAP : LUSOL,  LBSYMB
C    INFPRI : LUERR
C
#include <pcmdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL FILE_EXIST
#include "orgcom.h"
#include <priunit.h>
#include <infinp.h>
#include <pcm.h>
#include <pcmlog.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <inftap.h>
#include <infpri.h>
#include "infpar.h"
#include "qm3.h"
C
C
C
      CALL QENTER('IEFLNO')
C
C     Determine if full Hessian or only orbital Hessian
C
C
      IF (IPRRSP .GE. 40) THEN
         WRITE (LUPRI,'(//A)') ' >>> TEST OUTPUT FROM IEFLNO >>>'
      END IF
Clf      IF (IPRRSP .GE. 140) THEN
      IF (.true.) THEN
         WRITE (LUPRI,'(/A)') ' >>> IEFLNO - svec(1,nosim) on entry'
         CALL OUTPUT(SVEC,1,KZYVAR,1,NOSIM,KZYVAR,NOSIM,1,LUPRI)
      END IF
C
C     Arnfinn Nov. -08 ->
C     Read the value of NOSIM and NCOMS from the QM3 calculations 
C     (from the file MYFPCM).
      IF (MMPCM) THEN
         INQUIRE(FILE='MYFPCM',EXIST=FILE_EXIST)
         IF (FILE_EXIST) THEN
            LQM3PCM = -1
            CALL GPOPEN(LQM3PCM,'MYFPCM','OLD',' ','FORMATTED ',IDUMMY,
     &                  .FALSE.)
            REWIND(LQM3PCM)
            READ(LQM3PCM,*) NOSIMMM,NCOMS
Cahs         ELSE
Cahs            CALL QUIT('There is no input from QM/MM response 
Cahs     &                to PCM response')
         ENDIF
      ENDIF

      IF (NOSIMMM .NE. NOSIM) THEN
        CALL QUIT('NOSIM error in MM/PM interface IEF')
      ENDIF
C <- Arnfinn

C
C     Core allocation
C
      KVTEX  = 1
      KQTEX  = KVTEX + NTS*NOSIM
      KUBO   = KQTEX + NTS*NOSIM
      KW10   = KUBO  + NOSIM*N2ORBX
C
      KUCMO  = KW10
      KPCMX  = KUCMO  + NORBT*NBAST

      IF (TRPLET) THEN
         KPCMXT  = KPCMX  + NOSIM*N2ORBX
         KPCMXA  = KPCMXT + NOSIM*N2ORBX
         KPCMXAT = KPCMXA + NOSIM*N2ASHX
         KJPXAC  = KPCMXAT + NOSIM*N2ASHX
      ELSE
         KPCMXA = KPCMX  + NOSIM*N2ORBX
         KJPXAC = KPCMXA + NOSIM*N2ASHX
      ENDIF
      IF (TRPLET) THEN
         KJPXACT = KJPXAC + NOSIM
         KW20    = KJPXACT + NOSIM
      ELSE
         KW20   = KJPXAC + NOSIM
      ENDIF
C
      KJ1AO  = KW20
      KJ1SQ  = KJ1AO  + NNBASX*NSYM
      KJ1    = KJ1SQ  + N2ORBX
      KJ1XSQ = KJ1    + NNORBX
      KJ1XAC = KJ1XSQ + N2ORBX*NOSIM
      IF (TRPLET) THEN
         KJ1XACT = KJ1XAC   + NOSIM * N2ASHX
         KW29    = KJ1XACT  + NOSIM * N2ASHX
      ELSE
         KW29    = KJ1XAC   + NOSIM * N2ASHX
      ENDIF
C Arnfinn Nov. -08 ->
C Allocate memory and read coordinats and dipole moments (\mu^{[-1]}) from MM
      IF(MMPCM) THEN

         KMMCOR = KW29
         KMMDIP = KMMCOR + NOSIMMM * 3 * NCOMS
         KSTATE = KMMDIP + NOSIMMM * 3 * NCOMS
         KW30   = KSTATE + NOSIMMM

         CALL DZERO(WRK(KMMCOR),NOSIMMM*3*NCOMS)
         CALL DZERO(WRK(KMMDIP),NOSIMMM*3*NCOMS)

         IF (FILE_EXIST) THEN
            L = 0
            DO I = 1,NOSIMMM
               DO J = 1, NCOMS
                  L = L + 1
                  READ(LQM3PCM,'(6(E25.17,2x))') 
     &                 WRK(KMMCOR + 3 * ( L - 1 )     ),
     &                 WRK(KMMCOR + 3 * ( L - 1 ) + 1 ),
     &                 WRK(KMMCOR + 3 * ( L - 1 ) + 2 ),
     &                 WRK(KMMDIP + 3 * ( L - 1 )     ),
     &                 WRK(KMMDIP + 3 * ( L - 1 ) + 1 ),
     &                 WRK(KMMDIP + 3 * ( L - 1 ) + 2 )
Clf                  Write(lupri,*) 'transformed dipoles in rspief',
Clf     &                 l, j, i, ncoms, nosimm
Clf                  WRITE(LUPRI,'(6(E25.15,2x))') 
Clf     &                 WRK(KMMCOR + 3 * ( L - 1 )     ),
Clf     &                 WRK(KMMCOR + 3 * ( L - 1 ) + 1 ),
Clf     &                 WRK(KMMCOR + 3 * ( L - 1 ) + 2 ),
Clf     &                 WRK(KMMDIP + 3 * ( L - 1 )     ),
Clf     &                 WRK(KMMDIP + 3 * ( L - 1 ) + 1 ),
Clf     &                 WRK(KMMDIP + 3 * ( L - 1 ) + 2 )
               ENDDO
            ENDDO
            CALL GPCLOSE(LQM3PCM,'KEEP')
         ELSE
            CALL DZERO(WRK(KMMCOR),NOSIMMM*3*NCOMS)
            CALL DZERO(WRK(KMMDIP),NOSIMMM*3*NCOMS)
         ENDIF
      ENDIF

C <- Arnfinn
C
C     3.0 SOLSC
      KOVLP = KW30
      KW49   = KOVLP  + NOSIM
      LW49   = LFREE  + 1 - KW49
C
      KNEED = MAX(KW30,KW49)
      IF (KNEED .GT. LFREE) CALL ERRWRK('IEFLNO',KNEED,LFREE)
C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
C     Calculate unpacked orbital trial vectors in UBO
C
      IF (NOSIM.GT.0) THEN
         CALL RSPZYM(NOSIM,BOVECS,WRK(KUBO))
         CALL DSCAL(NOSIM*N2ORBX,-1.0D0,WRK(KUBO),1)
Clf         IF (IPRRSP .GE. 55) THEN
         IF (.true.) THEN
            DO 210 IOSIM = 1,NOSIM
               JUBO = KUBO + (IOSIM-1)*N2ORBX
               WRITE (LUPRI,2110) IOSIM,NOSIM
               CALL OUTPUT(WRK(JUBO),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                     LUPRI)
 210        CONTINUE
         END IF
      END IF
 2110 FORMAT (/,' Orbital trial vector unpacked to matrix form (no.',
     *        I3,' of',I3,')')
C
C     Contributions from all tesserae are included.
C
C
      CALL DZERO(WRK(KPCMX),NOSIM*N2ORBX)
      IF (TRPLET) CALL DZERO(WRK(KPCMXT),NOSIM*N2ORBX)
      REWIND LUPROP
      XI = DIPORG(1)
      YI = DIPORG(2)
      ZI = DIPORG(3)
#if defined (VAR_MPI)
      IF (NODTOT .GE. 1) THEN
         CALL J1XP(NOSIM,WRK(KVTEX),WRK(KUCMO),WRK(KUBO),UDV,UDVTR,
     &             WRK(KPCMX),WRK(KPCMXT),.FALSE.,TRPLET,JWOPSY,
     &             WRK(KW49),LW49)
      ELSE
#endif
C Arnfinn Nov. -08 ->
C Convert MM dipole moments to potentials
      IF (MMPCM) THEN
         KMMPOT = KW49
         KQNMM  = KMMPOT + NTSIRR
         KQMAT  = KQNMM + NTSIRR * NOSIMMM
         KW50   = KQMAT + NTS * NTSIRR
         LW50   = LFREE + 1 - KW50
         KNEED  = MAX(KW49,KW50)
         IF (KNEED .GT. LFREE) CALL ERRWRK('IEFLNO',KNEED,LFREE)

         CALL DZERO(WRK(KMMPOT),NTSIRR)
         CALL DZERO(WRK(KQNMM),NTSIRR*NOSIMMM)
         CALL DZERO(WRK(KQMAT),NTS*NTSIRR)
C
C false false true makes potential from dipole moments only
         CALL COMP_NUC_POT_CAV_MM(WRK(KMMPOT),.FALSE.,.FALSE.,
     &         .TRUE.,WRK(KMMCOR),WRK(KMMDIP),NOSIMMM)
C Convert potentials to charges, q
         CALL V2Q(WRK(KW50),WRK(KMMPOT),WRK(KQNMM),QNMMT,.FALSE.)
      ENDIF
C <- Arnfinn
      DO I = 1, NTSIRR
C Read AO potential integrals on tesserae
         L = 1
         NCOMP     = NSYM
         DIPORG(1) = XTSCOR(I)
         DIPORG(2) = YTSCOR(I)
         DIPORG(3) = ZTSCOR(I)
         EXP1VL    = .FALSE.
         TOFILE    = .FALSE.
         KPATOM    = 0
         TRIMAT    = .TRUE.
         CALL GET1IN(WRK(KJ1AO),'NPETES ',NCOMP,WRK(KW50),LW50,
     &               LABINT,INTREP,INTADR,L,TOFILE,KPATOM,TRIMAT,
     &               DUMMY,EXP1VL,DUMMY,IPRRSP)
         JJ1AO = KJ1AO
C
C     The totally symmetric case is always performed
C
Clara
C         do j=1, nsym 
C            WRITE (LUPRI,'(/A)') ' JJ1AO before UTHU'
C            CALL OUTPAK(WRK(JJ1AO+(j-1)*NNBASX),NBAST,1,LUPRI)
C         end do
C     Transform AO pot. int. into MO basis  V(AO) --> V(MO)
C
         CALL UTHU(WRK(JJ1AO),WRK(KJ1),WRK(KUCMO),WRK(KW30),
     *             NBAST,NORBT)
C Transform V(MO) from triangular to square format
         CALL DSPTSI(NORBT,WRK(KJ1),WRK(KJ1SQ))
         CALL DZERO(WRK(KJ1XSQ),N2ORBX*NOSIM)
         DO IOSIM = 1, NOSIM
            JUBO = KUBO + (IOSIM - 1) * N2ORBX
            JJ1X = KJ1XSQ + (IOSIM - 1) * N2ORBX
            JJ1XAC = KJ1XAC + (IOSIM - 1) * N2ASHX
            CALL ONEXH1(WRK(JUBO),WRK(KJ1SQ),WRK(JJ1X))
C     
            IF (NASHT .GT. 0) THEN
               CALL GETACQ(WRK(JJ1X),WRK(JJ1XAC))
            END IF
            IF (IPRRSP .GE. 15) THEN
               WRITE (LUPRI,'(/A)') ' J1X_mo matrix:'
               CALL OUTPUT(WRK(JJ1X),1,NORBT,1,NORBT,
     &              NORBT,NORBT,1,LUPRI)
               IF (NASHT .GT. 0) THEN
                  WRITE (LUPRI,'(/A)') ' J1X_ac matrix:'
                  CALL OUTPUT(WRK(JJ1XAC),1,NASHT,1,NASHT,
     &                 NASHT,NASHT,1,LUPRI)
               END IF
            END IF
C
C     Expectation value of transformed potential on tesserae:
C                     <0|\tilde{V}|0> 
C
            IF (KSYMOP .EQ. 1) THEN
               IF (TRPLET) THEN
                  FACTOR = SLVTLM(UDVTR,WRK(JJ1XAC),WRK(JJ1X),TJ1XAC)
               ELSE
                  FACTOR = SLVQLM(UDV,WRK(JJ1XAC),WRK(JJ1X),TJ1XAC)
               END IF
               WRK(KVTEX + (I-1) + (IOSIM-1)*NTS) = FACTOR
               IF (IPRRSP .GE. 6) THEN
                  WRITE (LUPRI,'(A,F17.8)')
     *                 ' --- J1X expectation value :',
     &                 WRK(KVTEX + (I-1) + (IOSIM-1)*NTS)
                  WRITE (LUPRI,'(A,F17.8)')
     *                 ' --- active part of J1X    :',TJ1XAC
               END IF
            END IF
         ENDDO
C Arnfinn Nov. -08 ->
C Added QSM 
         QFACTOR = QSN(I)+QSE(I)
         IF (MMPCM) QFACTOR = QSN(I)+QSE(I)+QSMM(I)
         IF (TRPLET) THEN
            CALL DAXPY(NOSIM*N2ORBX,QFACTOR,WRK(KJ1XSQ),1,WRK(KPCMXT),1)
         ELSE
            CALL DAXPY(NOSIM*N2ORBX,QFACTOR,WRK(KJ1XSQ),1,WRK(KPCMX),1)
         ENDIF
C KPCMX: \tilde{J} + \tilde{X}
C
C     For non-totally symmetric perturbation operators
C
         IF (KSYMOP .NE. 1) THEN
            ITS = (KSYMOP - 1)*NTSIRR + I
C     Transform AO pot. int. into MO basis  V(AO) --> V(MO)
            JJ1AO = KJ1AO + (KSYMOP - 1)*NNBASX
            CALL UTHU(WRK(JJ1AO),WRK(KJ1),WRK(KUCMO),WRK(KW50),
     *                 NBAST,NORBT)
C Transform V(MO) from triangular to square format
            CALL DSPTSI(NORBT,WRK(KJ1),WRK(KJ1SQ))
            CALL DZERO(WRK(KJ1XSQ),N2ORBX*NOSIM)
            DO IOSIM = 1, NOSIM
               JUBO = KUBO + (IOSIM - 1) * N2ORBX
               JJ1X = KJ1XSQ + (IOSIM - 1) * N2ORBX
               JJ1XAC = KJ1XAC + (IOSIM - 1) * N2ASHX
               CALL ONEXH1(WRK(JUBO),WRK(KJ1SQ),WRK(JJ1X))
C     
               IF (NASHT .GT. 0) THEN
                  CALL GETACQ(WRK(JJ1X),WRK(JJ1XAC))
               END IF
               IF (IPRRSP .GE. 15) THEN
                  WRITE (LUPRI,'(/A)') ' J1X_mo matrix :'
                  CALL OUTPUT(WRK(JJ1X),1,NORBT,1,NORBT,
     &                 NORBT,NORBT,1,LUPRI)
                  IF (NASHT .GT. 0) THEN
                     WRITE (LUPRI,'(/A)') ' J1X_ac matrix :'
                     CALL OUTPUT(WRK(JJ1XAC),1,NASHT,1,NASHT,
     &                    NASHT,NASHT,1,LUPRI)
                  END IF
               END IF
C
C     Expectation value of transformed potential on tesserae:
C                     <0|\tilde{V}|0> 
C
               IF (TRPLET) THEN
                  FACTOR = SLVTLM(UDVTR,WRK(JJ1XAC),WRK(JJ1X),TJ1XAC)
               ELSE
                  FACTOR = SLVQLM(UDV,WRK(JJ1XAC),WRK(JJ1X),TJ1XAC)
               END IF
               WRK(KVTEX + (ITS-1) + (IOSIM-1)*NTS) = FACTOR
               IF (IPRRSP .GE. 6) THEN
                  WRITE (LUPRI,'(A,F17.8)')
     *                 ' --- J1X expectation value :',
     &                 WRK(KVTEX + (ITS-1) + (IOSIM-1)*NTS)
                  WRITE (LUPRI,'(A,F17.8)')
     *                 ' --- active part of J1X    :',TJ1XAC
               END IF
            ENDDO
C KPCMX: \tilde{J} + \tilde{X}
         END IF
      ENDDO
#if defined (VAR_MPI)
      END IF
#endif
Clf
         write(lupri,*) 'MM surface charges RESPONSE'
         do i=1,nts
            write(lupri,*) i, qsmm(i)
         enddo
C
Clf      IF (IPRRSP .GE. 50) THEN
      IF (.true.) THEN
         DO IOSIM = 1,NOSIM
            JPCMX = KPCMX + (IOSIM-1)*N2ORBX
            WRITE (LUPRI,'(/A,I3,A,I3)')
     *           ' >>> IEFLNO - (JPCMX) half.',IOSIM,' of',NOSIM
            CALL OUTPUT(WRK(JPCMX),1,NORBT,1,NORBT,NORBT,NORBT,
     *                  1,LUPRI)
         END DO
      END IF
      DO IOSIM = 1, NOSIM
         CALL V2Q(WRK(KW50),WRK(KVTEX + NTS*(IOSIM-1)),
     &        WRK(KQTEX + NTS*(IOSIM-1)),
     &        QTEXS,NEQRSP)
         IF (IPRRSP .GE. 6) THEN
            DO I = 1, NTS
               WRITE (LUPRI,'(A,F17.8)')
     *              ' --- J2X expectation value :',
     &          WRK(KQTEX + NTS*(IOSIM-1) + I - 1)
            END DO
         END IF
      ENDDO
      CALL GPCLOSE(LUPCMD,'KEEP')
C
C Arnfinn Nov. -08 ->
C Make the file QFQM3 with charges for use by QM3 response calculations
      IF (MMPCM) THEN
         INQUIRE(FILE='QFQM3',EXIST=FILE_EXIST)
         IF (FILE_EXIST) THEN
            LPCMQM3 = -1
            CALL GPOPEN(LPCMQM3,'QFQM3','OLD',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
            CALL GPCLOSE(LPCMQM3,'DELETE')
         ENDIF
         LPCMQM3 = -1
         CALL GPOPEN(LPCMQM3,'QFQM3','NEW',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
         WRITE (LPCMQM3,*) NOSIMMM, NTS
         DO I = 1, NOSIMMM
            NSTATE = INT(WRK(KSTATE+I-1))
            WRITE (LPCMQM3,*) NSTATE
            DO ITS = 1, NTS
C luca this will not work with more than one perturbation 
C (see above for how WRK(KQNMM) is created)
               INDEX = NTS * (I-1) + ITS - 1
               WRK(KQTEX + INDEX) = 
     &             WRK(KQTEX + INDEX) - 2.0D0 * WRK(KQNMM + INDEX)
               WRITE (LPCMQM3,'(4(E25.15,2x))') XTSCOR(ITS),YTSCOR(ITS),
     &              ZTSCOR(ITS), WRK(KQTEX + INDEX) 
            ENDDO
         ENDDO
         CALL GPCLOSE(LPCMQM3,'KEEP')
      ENDIF
C <- Arnfinn
C
C     TRANSFORMED CHARGES MULTIPIED TO POTENTIALS.
C     There are only one-index-transformed charges for the totally
C     symmetric irrep.
C
      CALL DZERO(WRK(KJ1XSQ),NOSIM*NNBASX)
      IF (KSYMOP .GT. 1) THEN
         DO IOSIM = 1, NOSIM
            KTSSIM = (IOSIM-1) * NTS
            KTSSYMOP = (KSYMOP-1) * NTSIRR
            CALL DCOPY(NTSIRR,WRK(KQTEX + KTSSIM + KTSSYMOP),1,
     &                 WRK(KQTEX + KTSSIM))
            CALL DZERO(WRK(KQTEX + KTSSIM + KTSSYMOP), NTSIRR)
         END DO
      END IF
      CALL J1INT(WRK(KQTEX),.FALSE.,WRK(KJ1XSQ),NOSIM,.FALSE.,
     &           'NPETES ',KSYMOP,WRK(KW50),LW50)
      DO IOSIM = 1, NOSIM
         JAOP  = KJ1XSQ + (IOSIM - 1) * NNBASX
         JPCMX = KPCMX  + (IOSIM - 1) * N2ORBX
         CALL UTHU(WRK(JAOP),WRK(KJ1),WRK(KUCMO),WRK(KW50),
     *             NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KJ1),WRK(KJ1SQ))
         CALL DAXPY(N2ORBX,-D1,WRK(KJ1SQ),1,
     &              WRK(JPCMX),1)
      END DO
C
C     Restore dipole origin
C
      DIPORG(1) = XI
      DIPORG(2) = YI
      DIPORG(3) = ZI
C
Clf      IF (IPRRSP .GE. 50) THEN
      IF (.true.) THEN
         DO IOSIM = 1,NOSIM
            JPCMX = KPCMX + (IOSIM-1)*N2ORBX
            WRITE (LUPRI,'(/A,I3,A,I3)')
     *           ' >>> IEFLNO - (JPCMX) matrix no.',IOSIM,' of',NOSIM
            CALL OUTPUT(WRK(JPCMX),1,NORBT,1,NORBT,NORBT,NORBT,
     *                  1,LUPRI)
         END DO
      END IF
C
Clara add the case of KZCONF.GT.0 
      IF ((.NOT.TDHF).AND.(KZCONF.GT.0)) THEN
C
C        ... CSF part of sigma vectors
C
         DO 700 IOSIM = 1,NOSIM
            JPCMX   = KPCMX   + (IOSIM-1)*N2ORBX
            JPCMXA  = KPCMXA  + (IOSIM-1)*N2ASHX
            IF (TRPLET) THEN
               JPCMXT  = KPCMXT  + (IOSIM-1)*N2ORBX
               JPCMXAT = KPCMXAT + (IOSIM-1)*N2ASHX
            ENDIF
            IF (IREFSY .EQ. KSYMST) THEN
               WRK(KOVLP-1+IOSIM) = DDOT(KZCONF,CREF,1,SVEC(1,IOSIM),1)
            ELSE
               WRK(KOVLP-1+IOSIM) = D0
            END IF
            CALL GETACQ(WRK(JPCMX),WRK(JPCMXA))
C <0|\tilde{Z}|0>  (Z=J+X+V*Q just to remind)
            TJPX   = SLVQLM(UDV,WRK(JPCMXA),WRK(JPCMX),TJPXAC)
            WRK(KJPXAC-1+IOSIM) = TJPXAC
            IF (TRPLET) THEN
               CALL GETACQ(WRK(JPCMXT),WRK(JPCMXAT))
               TJPXT  = SLVQLM(UDVTR,WRK(JPCMXAT),WRK(JPCMXT),TJPXACT)
               WRK(KJPXACT-1+IOSIM) = TJPXACT
            ELSE
               IF (IPRRSP .GE. 40) WRITE (LUPRI,'(/A,I5,3F15.10)')
     *              ' IOSIM, C_OVLP, TJPXAC, TJPX :',
     *              IOSIM,WRK(KOVLP-1+IOSIM),TJPXAC,TJPX
            ENDIF
 700     CONTINUE
C  \sigma_{co} = <i|\tilde{Z}|0>
C KW30???
         CALL SLVSC(0,NOSIM,N2ASHX,DUMMY,CREF,SVEC,WRK(KPCMXA),
     *        WRK(KPCMXAT),
     *              WRK(KJPXAC),DUMMY,INDXCI,WRK(KW50),LW50)
C        CALL SLVSC(NCSIM,NOSIM,BCVECS,CREF,SVECS,
C    *              RXAC,RYAC,TRXAC,TRYAC,INDXCI,WRK,LWRK)
Clf         IF (IPRRSP.GT.101) THEN
         IF (.true.) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
            WRITE(LUPRI,*)' **** AFTER SLVSC in IEFLNO **** '
            CALL OUTPUT(SVEC,1,KZYVAR,1,NOSIM,KZYVAR,NOSIM,1,LUPRI)
         END IF
      END IF
C \sigma_{oo} = <0|[q_j,\tilde{Z}]|0> 
      IF (.NOT. LADDMM) THEN
        CALL DZERO(WRK(KPCMX),NOSIM*N2ORBX)
        IF (TRPLET) CALL DZERO(WRK(KPCMXT),NOSIM*N2ORBX)
      ENDIF

      IF(TRPLET) THEN
         CALL SLVSOR(.TRUE.,.FALSE.,NOSIM,UDVTR,SVEC(1,1),WRK(KPCMX))
         CALL SLVSOR(.TRUE.,.TRUE., NOSIM,UDV,  SVEC(1,1),WRK(KPCMXT))
      ELSE
         CALL SLVSOR(.TRUE.,.TRUE., NOSIM,UDV,  SVEC(1,1),WRK(KPCMX))
      ENDIF
      IF (KZCONF.GT.0 .AND. IREFSY.EQ.KSYMST) THEN
         IF (NCREF .NE. KZCONF) CALL QUIT('IEFLNO: NCREF .ne. KZCONF')
C
C        ... test orthogonality
C
         DO 900 IOSIM = 1,NOSIM
            TEST = DDOT(KZCONF,CREF,1,SVEC(1,IOSIM),1)
     *           - WRK(KOVLP-1+IOSIM)
            IF (ABS(TEST) .GT. 1.D-8) THEN
               NWARN = NWARN + 1
               WRITE (LUPRI,'(/A,I5,/A,1P,D12.4)')
     *              ' >>> IEFLNO WARNING, for IOSIM =',IOSIM,
     *              ' <CREF | SVEC_solvent(iosim) > =',TEST
            END IF
 900     CONTINUE
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        ... test print
C
Clf      IF (IPRRSP .GE. 140) THEN
      IF (.true.) THEN
         WRITE (LUPRI,'(/A)') ' >>> IEFLNO - svec(ci,1) on exit'
         DO 930 I = 1,KZCONF
            IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *         ' conf #',I,SVEC(I,1)
  930    CONTINUE
      END IF
Clf      IF (IPRRSP .GE. 140) THEN
      IF (.true.) THEN
         WRITE (LUPRI,'(/A)') ' >>> IEFLNO - svec(orb,1) on exit'
         WRITE (LUPRI,'(/A)') ' Z - PART OR VECTOR '
         CALL OUTPUT(SVEC(KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Y - PART OR VECTOR '
         CALL OUTPUT(SVEC(KZVAR+KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,
     *               1,1,LUPRI)
      END IF
C
      CALL QEXIT('IEFLNO')
      RETURN
C     ... end of IEFLNO.
      END

C  /* Deck pcm_comp_pot */
      SUBROUTINE PCM_COMP_POT(MOCOEF, DENMAT, DENTRP, UBO, XJSQ, XJSQT, 
     &                        NOSIM, JWOPSY, KSYMP, WORK, LWORK)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <thrzer.h>
#include <iratdef.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
#include <qm3.h>
#include <orgcom.h>
#include <dftcom.h>
#include <inforb.h>
#include <infrsp.h>

      IF (VSNFLG .EQ. -1) THEN
         CALL DZERO(VSN, NTS)
         CALL PCM_COMP_NUC_POT(VSN)
         VSNFLG = 1
      END IF
      IF (VSEFLG .EQ. -1) THEN 
         CALL DZERO(VSE, NTS)
         CALL PCM_COMP_EL_POT(VSE, DENMAT, NOSIM, NBAS, NSYM, NNBASX, 
     &                        KSYMP, WORK, LWORK)
         VSEFLG = 1
      END IF
      IF (VSE1FLG .EQ. -1) THEN 
         CALL DZERO(VSE1, NTS)
c      SUBROUTINE PCM_COMP_EL_POT1(VSE1, NOSIM, DENMAT, DENTRP, JWOPSY, 
c     &                            UCMO, BOVECS, XJSQ, XJSQT, TRPLET, 
c     &                            WORK, LWORK)
         CALL PCM_COMP_EL_POT1(VSE1, NOSIM, DENMAT, DENTRP, JWOPSY, 
     &                         UCMO, BOVECS, XJSQ, XJSQT, TRPLET,
     &                         WORK, LWORK)
         VSE1FLG = 1
      END IF
      IF (VSMMFLG .EQ. -1) THEN
         CALL DZERO(VSMM, NTS)
         CALL PCM_COMP_MM_POT(VSMM)
         VSMMFLG = 1
      END IF
      IF (VSMMFLG .EQ. -1) THEN
         CALL DZERO(VSMM1, NTS)
         CALL PCM_COMP_MM_POT(VSMM1)
         VSMMFLG = 1
      END IF
      RETURN
      END

C  /* Deck pcm_comp_nuc_pot */
      SUBROUTINE PCM_COMP_NUC_POT(VSN)
#include <implicit.h>
      CALL COMP_NUC_POT_CAV(VSN, .TRUE., .FALSE.,.FALSE.)
      RETURN
      END

C /* Deck pcm_comp_el_pot */
      SUBROUTINE PCM_COMP_EL_POT(VSE, DENMAT, NOSIM, NBAS, NSYM, NNBASX,
     &                           KSYMP, WORK, LWORK)
#include <implicit.h>
      DIMENSION WORK(*), DENMAT(*), VSE(*)
      KDEN = 1
      KFREE = KDEN + NNBASX
      LFREE = LWORK - KFREE + 1
      IF (LFREE .LT. 0) THEN
         CALL QUIT('Not enough memory in pcm_comp_el_pot')
      ENDIF
      
      print *,'luca',nbas,nsym,nnbasx,nosim,ksymp
      CALL PKSYM1(WORK(KDEN),DENMAT,NBAS,NSYM,-1)

      CALL J1INT(VSE, .TRUE., WORK(KDEN), NOSIM, .FALSE., 'NPETES ',
     &                 KSYMP, WORK(KFREE), LFREE)
      RETURN
      END

C /* Deck pcm_comp_el_pot1 */
      SUBROUTINE PCM_COMP_EL_POT1(VSE1, NOSIM, DENMAT, DENTRP, JWOPSY, 
     &                            UCMO, BOVECS, XJSQ, XJSQT, TRPLET, 
     &                            WORK, LWORK)
#include <implicit.h>
      LOGICAL TRPLET
      DIMENSION BOVECS(*)

      KUBO = 1
      KFREE = KUBO + NOSIM * N2ORBX
      LFREE = LWORK - KFREE + 1
      IF (LFREE .LT. 0) THEN 
         CALL QUIT('PCM_COMP_EL_POT1: insufficient memory')
      END IF
      
      IF (NOSIM.GT.0) THEN
         CALL RSPZYM(NOSIM,BOVECS,WORK(KUBO))
         CALL DSCAL(NOSIM*N2ORBX,-1.0D0,WORK(KUBO),1)
      END IF

#if defined (VAR_MPI)
      CALL J1XP(NOSIM,VSE1,UCMO,WORK(KUBO),DENMAT,DENTRP,XJ1SQ,XJ1SQT,
     &     .FALSE.,TRPLET,JWOPSY,WORK,LWORK)
#else
      CALL J1X(NOSIM,VSE1,UCMO,WORK(KUBO),DENMAT,UDVTR,XJ1SQ,XJ1SQT,
     &     .FALSE.,TRPLET,JWOPSY,WORK,LWORK)
#endif
      SUBROUTINE J1X(NOSIM,VTEX,UCMO,UBO,UDV,UDVTR,XJ1SQ,XJ1SQT,
     &                TOFILE,TRPLET,JWOPSY,WORK,LWORK)
      RETURN
      END

C /* Deck pcm_comp_mm_pot */
      SUBROUTINE PCM_COMP_MM_POT(VCAV, NTS, XTS, YTS, ZTS, XMM, YMM, 
     &                           ZMM, SOURCE, NSOURCES, TYPE)
#include <implicit.h>
      CHARACTER*3 TYPE
      IF (TYPE .EQ. 'CHG') THEN
         CALL PCM_COMP_POT_CHG(VCAV, NTS, XTS, YTS, ZTS, XMM, YMM,
     &                         ZMM, SOURCE, NSOURCES)
      ELSE IF (TYPE .EQ. 'DIP') THEN
         CALL PCM_COMP_POT_DIP(VCAV, NTS, XTS, YTS, ZTS, XMM, YMM,
     &                         ZMM, SOURCE, NSOURCES)
Clf   quadrupole not yet imlemented
Clf      ELSE IF (TYPE .EQ. 'QDR') THEN
Clf         PCM_COMP_POT_QDR(VCAV, NTS, COORD, SOURCE, NSOURCES)
Clf   quadrupole components order xx, xy, xz, yy, yz, zz
Clf   quadrupole potential: sum_{a,b} 1/2*(3*Ra*Rb-R^2*delta_ab)/R^5 Q_ab
Clf
      ELSE
         CALL QUIT('Wrong type in PCM_COMP_MM_POT')
      END IF
      RETURN
      END

C /* Deck pcm_comp_mm_pot */
      SUBROUTINE PCM_COMP_POT_CHG(VCAV, NTS, XTS, YTS, ZTS, XMM, YMM,
     &                            ZMM, SOURCE, NSITES)
#include <implicit.h>
      DIMENSION VCAV(*), XTS(*), YTS(*), ZTS(*),XMM(*),YMM(*),ZMM(*),
     &          SOURCE(*)
      DO IATOM = 1, NSITES
         XL = XMM(IATOM)    
         YL = YMM(IATOM)    
         ZL = ZMM(IATOM)
         CHGMM = SOURCE(IATOM)
         DO ITS = 1, NTS
            RIL=DSQRT((XTS(ITS)-XL)**2+(YTS(ITS)-YL)**2
     &           +(ZTS(ITS)-ZL)**2) 
            VCAV(ITS) = VCAV(ITS) + CHGMM/RIL
         ENDDO
      ENDDO
      RETURN
      END
      
C /* Deck pcm_comp_pot_dip */
      SUBROUTINE PCM_COMP_POT_DIP(VCAV, NTS, XTS, YTS, ZTS, XMM, YMM,
     &                            ZMM, SOURCE, NSITES)
#include <implicit.h>
      DIMENSION VCAV(*), XTS(*), YTS(*), ZTS(*),XMM(*),YMM(*),ZMM(*),
     &          SOURCE(*)
      DO IATOM = 1, NSITES
         XDIP = SOURCE(3*IATOM-2)
         YDIP = SOURCE(3*IATOM-1)
         ZDIP = SOURCE(3*IATOM)
         DO ITS = 1, NTS
            XL = XTS(ITS) - XMM(IATOM)
            YL = YTS(ITS) - YMM(IATOM)
            ZL = ZTS(ITS) - ZMM(IATOM)
            SCALAR = XL * XD + YL * YD + ZL * ZD
            R2 = XL ** 2 + YL ** 2 + ZL ** 2
            R3 = (DSQRT(R2)) ** 3
            VCAV(ITS) = VCAV(ITS) + SCALAR/R3
         ENDDO
      ENDDO
      RETURN
      END
