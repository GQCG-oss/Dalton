
      SUBROUTINE WAVPCM_INIT()

#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <maxorb.h>
#include <wavpcm.h>
#include <infpar.h>

      CALL TIMER('START ',TIMSTR,TIMEND)

C Write coordinates for the cavity in Angstroms
      AU2ANG = 0.529177249D0

      OPEN(101,file='cavity.inp')
      WRITE(101,*) NCTOT
      DO I = 1, NCTOT
         WRITE(101,*) CORD(1,I)*AU2ANG, CORD(2,I)*AU2ANG, 
     &        CORD(3,I)*AU2ANG, WAVSPH(I)
      END DO
      CLOSE(101)

      OPEN(101,file='nuclei.dat')
      WRITE(101,*) NCTOT, 'au'
      DO I = 1, NCTOT
         WRITE(101,*) CORD(1,I), CORD(2,I), CORD(3,I), CHARGE(I)
      END DO
      CLOSE(101)

      IF(NODTOT .GT. 1) THEN
         CALL MPIXBCAST(47,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)
      END IF

      WRITE(LUPRI,*) 'INITINITINIT'
      CLOSE(LUPRI)
      CALL WAVPCM_INITFORREAL()
      
      OPEN(UNIT=LUPRI,FILE='DALTON.OUT',STATUS='OLD',ACCESS='APPEND')
      WRITE(LUPRI,*) ''
      CALL TIMER('OPERAT',TIMSTR,TIMEND)
      
      END

      
      SUBROUTINE WAVPCM_INITFORREAL()
      
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <maxorb.h>
#include <wavpcm.h>
#include <infpar.h>

C     
C Initialize the wavelet BEM operator and PCM stuff.
C     
      
      IF(NODTOT .GT. 1) THEN
         CALL MPIXBCAST(NBAST, 1, 'INTEGER', MASTER)
         CALL MPIXBCAST(NOCCT, 1, 'INTEGER', MASTER)
         CALL MPIXBCAST(NVIRT, 1, 'INTEGER', MASTER)
      END IF

C Make sure that every slave agrees
      WAVPCM = .TRUE.

#ifdef USE_WAVPCM      
      CALL pcm_init(NBAST,NOCCT,NVIRT)      
#endif      
      END
      

C
      SUBROUTINE WAVPCM_SCF(DCAO, DVAO, FSOL, ESOLT, CMO, WRK, LWRK)
C
C DCAO(in) Density matrix for closed shells
C DVAO(in) Density matrix for open shells in ROHF (or something)
C FSOL(out) Solvation contribution to the Fock matrix
C ESOLT(out) Solvation energy
C CMO(in) Molecular orbitals
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <maxorb.h>
#include <wavpcm.h>
#include <dftcom.h>
#include <infpar.h>
C
      DIMENSION DCAO(*), DVAO(*), FSOL(*), CMO(*), WRK(LWRK)

C
C First save the density matrix to be used in the potential evaluations. 
C If ROHF stuff is needed, then do something with the DVAO also.
C
      CALL PKSYM1(WAVDEN,DCAO,NBAS,NSYM,-1)
C
C Then call PCM stuff
C
      IF(NODTOT .GT. 1) THEN
         CALL MPIXBCAST(48,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)
      END IF

      WRITE(LUPRI,*) ''
      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL WAVPCM_SCFFORREAL()
      CALL TIMER('WAVPOL',TIMSTR,TIMEND)
      WRITE(LUPRI,*) ''

C Get the results from the parallel part
      ESOLT = WAVSCFE
      DO I = 1, NNORBT
         FSOL(I) = WAVFCK(I)
      END DO

C Save the Fock matrix (in MO basis) for future use in solvent calculations
      DO I = 1, NNORBT
         WAVFCK(I) = 0D0
      END DO
      CALL UTHU(FSOL, WAVFCK, CMO, WRK, NBAS(1), NORB(1))

      OPEN(101,FILE='dmat.dat')
      WRITE(101,*) WAVDEN(1:NNORBT)
      CLOSE(101)

      END      


      SUBROUTINE WAVPCM_SCFFORREAL()

#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <maxorb.h>
#include <wavpcm.h>
#include <dftcom.h>
#include <infpar.h>

      DIMENSION FSOL(MXORBT)

      CALL MPIXBCAST(WAVDEN,NNORBT,'DOUBLE',MASTER)

      ESOLT = 0D0
#ifdef USE_WAVPCM      
      CALL pcm_polarization(ESOLT, FSOL)
#endif
      WAVSCFE = ESOLT
      DO I = 1, NNORBT
         WAVFCK(I) = FSOL(I)
      END DO

      END




C
C Linear response
C
      SUBROUTINE WAVPCM_LINRSP(NOSIM, ZYOVEC, XINDX, UDV, 
     *        EVECS, WRK, LWRK1)
C
C NOSIM(in) - Number of trial vectors
C ZYOVEC(KZYVAR,NOSIM)(in) - Trial vector(s)
C EVECS(KZYVAR,NOSIM)(inout) - Results
C XINDX - ?
C UDV - ?
C
#include <implicit.h>
#include <maxorb.h>
#include <infpri.h>
#include <infrsp.h>
#include <inftap.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <wavpcm.h>
#include <wrkrsp.h>
#include <infpar.h>

      DIMENSION ZYOVEC(*)
      DIMENSION EVECS(*)
      DIMENSION XINDX(*), UDV(*), WRK(*)

C Write down the orbitals and orbital energies
      KWORK = 1
      KMO = KWORK
      KWORK = KWORK + 2*N2BAST
      CALL READMO(WRK(KMO), 9, WRK(KWORK), LWRK1-KWORK)
C Orbital energies should be in WAVOE
      OPEN(101,FILE='MO.dat')
      KOEIND = 0
      DO I = 1, NOCC(1)+NVIR(1)
         WRITE(101,*) WAVOE(I)
         DO J = 1, NBAST
            WRITE(101,*) WRK(KMO+(I-1)*NBAST+J-1)
         END DO
      END DO
      CLOSE(101) 

      IF(NODTOT .GT. 1) THEN

         CALL MPIXBCAST(49,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)         

         CALL MPIXBCAST(NOSIM, 1, 'INTEGER', MASTER)
         CALL MPIXBCAST(KZYVAR, 1, 'INTEGER', MASTER)

         CALL MPIXBCAST(ZYOVEC, KZYVAR*NOSIM, 'DOUBLE', MASTER)
         CALL MPIXBCAST(EVECS, KZYVAR*NOSIM, 'DOUBLE', MASTER)
         
         CALL MPIXBCAST(TRPLET, 1, 'LOGICAL', MASTER)

      END IF

C Assume that density matrices & MO coefficients from the last SCF iteration
C are still valid
      CALL TIMER('START ',TIMSTR,TIMEND)
#ifdef USE_WAVPCM      
      IF (TRPLET) THEN
         call pcm_linrsp(KZYVAR, NOSIM, ZYOVEC, EVECS, 1)
      ELSE 
         call pcm_linrsp(KZYVAR, NOSIM, ZYOVEC, EVECS, 0)
      END IF
#endif
      CALL TIMER('WAVRSP',TIMSTR,TIMEND)

      END

C
C Linear response
C
      SUBROUTINE WAVPCM_LINRSP_SLAVE(WORK, LWORK)
C
C NOSIM(in) - Number of trial vectors
C ZYOVEC(KZYVAR,NOSIM)(in) - Trial vector(s)
C EVECS(KZYVAR,NOSIM)(inout) - Results
C XINDX - ?
C UDV - ?
C
#include <implicit.h>
#include <maxorb.h>
#include <infpri.h>
#include <infrsp.h>
#include <inftap.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <wavpcm.h>
#include <wrkrsp.h>
#include <infpar.h>

      DIMENSION WORK(*)

      INTEGER ZYOVEC, EVECS

      CALL MPIXBCAST(NOSIM, 1, 'INTEGER', MASTER)
      CALL MPIXBCAST(KZYVAR, 1, 'INTEGER', MASTER)

      ZYOVEC = 1
      EVECS = KZYVAR*NOSIM + 1

      CALL MPIXBCAST(WORK(ZYOVEC), KZYVAR*NOSIM, 'DOUBLE', MASTER)
      CALL MPIXBCAST(WORK(EVECS), KZYVAR*NOSIM, 'DOUBLE', MASTER)
         
      CALL MPIXBCAST(TRPLET, 1, 'LOGICAL', MASTER)
      
C Assume that density matrices & MO coefficients from the last SCF iteration
C are still valid
#ifdef USE_WAVPCM      
      IF (TRPLET) THEN
         call pcm_linrsp(KZYVAR, NOSIM, WORK(ZYOVEC), WORK(EVECS), 1)
      ELSE 
         call pcm_linrsp(KZYVAR, NOSIM, WORK(ZYOVEC), WORK(EVECS), 0)
      END IF
#endif
      END




C
C Driver for J1INT, calculates the electronic potential at some position
C
      SUBROUTINE DALTON_PCM_ELECTRON_POTENTIAL(XPOS, YPOS, ZPOS, RES, 
     &     WORK, LWORK)
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <codata.h>
#include <inforb.h>
#include <pcmnuclei.h>
#include <pcm.h>
#include <pcmlog.h>
#include <maxorb.h>
#include <wavpcm.h>
C     
      DIMENSION WORK(*)
C
      NTSIRROLD = NTSIRR
      NTSIRR = 1
      XTMP = XTSCOR(1)
      YTMP = YTSCOR(1)
      ZTMP = ZTSCOR(1)
      XTSCOR(1) = XPOS
      YTSCOR(1) = YPOS
      ZTSCOR(1) = ZPOS
C
      CALL J1INT(RES,.TRUE.,WAVDEN,1,.FALSE.,'NPETES ',1,
     &     WORK,LWORK)
C
      XTSCOR(1) = XTMP
      YTSCOR(1) = YTMP
      ZTSCOR(1) = ZTMP
      NTSIRR = NTSIRROLD
C
      RETURN
      END
C
C Driver for J1INT, calculates the electronic potential integral
C
      SUBROUTINE DALTON_PCM_INTEGRALS(XPOS, YPOS, ZPOS, RES, WORK,LWORK)
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <codata.h>
#include <inforb.h>
#include <pcmnuclei.h>
#include <pcm.h>
#include <pcmlog.h>
#include <maxorb.h>
#include <wavpcm.h>
C
      DIMENSION WORK(*)
      DIMENSION PCMTMP(1)
C     
      NTSIRROLD = NTSIRR
      NTSIRR = 1
      NTS = 1
      XTMP = XTSCOR(1)
      YTMP = YTSCOR(1)
      ZTMP = ZTSCOR(1)
      XTSCOR(1) = XPOS
      YTSCOR(1) = YPOS
      ZTSCOR(1) = ZPOS
      PCMTMP = 1D0
C
      CALL J1INT(PCMTMP,.FALSE.,RES,1,.FALSE.,'NPETES ',1,
     &     WORK,LWORK)
C
      XTSCOR(1) = XTMP
      YTSCOR(1) = YTMP
      ZTSCOR(1) = ZTMP
      NTSIRR = NTSIRROLD
C
      RETURN
      END


      SUBROUTINE WAVPCM_2EL(FOCK, DMAT, WRK, LWRK)

#include <implicit.h>
#include <maxorb.h>
#include <infpri.h>
#include <infrsp.h>
#include <inftap.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <wavpcm.h>
#include <wrkrsp.h>
#include <infpar.h>

      DIMENSION FOCK(*)
      DIMENSION DMAT(*)
      DIMENSION WRK(*)

C Slaves up
      IF(NODTOT .GT. 1) THEN
         CALL MPIXBCAST(50,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)         
      END IF

      CALL MPIXBCAST(DMAT, N2BASX, 'DOUBLE', MASTER)

#ifdef USE_WAVPCM      
      CALL pcm_ao_rsp_2elfock(FOCK, DMAT)
#endif      
      RETURN
      END


      SUBROUTINE WAVPCM_2EL_SLAVE(WORK, LWORK)

#include <implicit.h>
#include <maxorb.h>
#include <infpri.h>
#include <infrsp.h>
#include <inftap.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <wavpcm.h>
#include <wrkrsp.h>
#include <infpar.h>

      DIMENSION WORK(*)
      
      IDMAT = 1
      IFOCK = IDMAT + N2BASX
      CALL MPIXBCAST(WORK(IDMAT), N2BASX, 'DOUBLE', MASTER)

#ifdef USE_WAVPCM      
      CALL pcm_ao_rsp_2elfock(WORK(IFOCK), WORK(IDMAT))
#endif
      END
