C/* deck MAGPCMFIRST */
      SUBROUTINE MAGPCMFIRST(VEC2,NCOMP,IPRINT,INTPRI,WORK,LWORK)
C
c This routine add the PCM solvent contribution 
c to the 1st order (NCOMP=3) 
c magnetic field perturbation (see Cammi,JCP,109:3185).
c By Kenneth Ruud&Domenico Marchesan - Feb 2004
c
c VEC2 (OUTPUT)  : contribution to be added to perturbed Fock matrix
c NCOMP (INPUT)  : number of independent perturbation (3)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <pcmdef.h>
#include <orgcom.h>
#include <pcm.h>
#include <pcmlog.h>
#include <inforb.h>
#include <infpri.h>
#include <infvar.h>
C
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT
      DIMENSION WORK(*),VEC2(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      INTEGER BAS,NCOMP
c     mxcent : maximum number of atoms

c     perturbed potential
      KEXPVL = 1
      KLAST  = KEXPVL + NCOMP*NTS
      LWRK   = LWORK - KLAST + 1
      CALL DZERO(VEC2,3*N2BASX)
C
C     Send in three copies of the sum of electronic and nuclear charges
C
      CALL DCOPY(NTS,QSE,1,WORK(KEXPVL),1)
      CALL DSCAL(NTS,-1.0D0,WORK(KEXPVL),1)
      CALL DAXPY(NTS,-1.0D0,QSN,1,WORK(KEXPVL),1)
      CALL DCOPY(NTS,WORK(KEXPVL),1,WORK(KEXPVL +   NTS),1)
      CALL DCOPY(NTS,WORK(KEXPVL),1,WORK(KEXPVL + 2*NTS),1)
      CALL J1INT(WORK(KEXPVL),.FALSE.,VEC2,NCOMP,.FALSE.,'PCMBSOL',
     &           1,WORK(KLAST),LWRK)

c     **print out section**

      IF (IPRINT .GT. 5) THEN
         CALL AROUND(
     &        'First order solvent contributions to gradient in MAGPCM')
         WRITE (LUPRI,'(2X,A)') 'X coordinate'
         CALL OUTPUT(VEC2(1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE (LUPRI,'(2X,A)') 'Y coordinate'
         CALL OUTPUT(VEC2(1+N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,1,
     &        LUPRI)
         WRITE (LUPRI,'(2X,A)') 'Z coordinate'
         CALL OUTPUT(VEC2(1+2*N2BASX),1,NBAST,1,NBAST,
     &        NBAST,NBAST,1,LUPRI)
      END IF
C
      RETURN
      END


C/* Deck MAGPCMSEC */
      SUBROUTINE MAGPCMSEC(VEC2,NCOMP,IPRINT,INTPRI,WORK,LWORK,DENMAT
     &     ,POTPER,QSEPER)
c
c This routine add the PCM solvent contribution 
c to the the 2nd order (NCOMP=6) 
c magnetic field perturbation (see Cammi,JCP,109:3185).
c Domenico Marchesan&Kenneth Ruud - Feb 2004
c 
c VEC2 (OUTPUT)  : contribution to be added to perturbed Fock matrix
c NCOMP (INPUT)  : number of independent perturbation (6)
c DENMAT (INPUT) : density matrix
c POTPER (OUTPUT): perturbed potential. NOT used in ABAMAG
c QSEPER (OUTPUT): perturbed charges. NOT used in ABAMAG

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <pcmdef.h>
#include <orgcom.h>
#include <pcmnuclei.h>
#include <pcm.h>
#include <pcmlog.h>
#include <inforb.h>
#include <infpri.h>
#include <infvar.h>
#include <inftap.h>

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT
      DOUBLEPRECISION VEC2,DENMAT,QSEPER,POTPER
      DIMENSION WORK(*),VEC2(*),DENMAT(*)
      DIMENSION POTPER(*),QSEPER(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION NPTCON(MXTS)
      INTEGER NCOMP,INUC

c     MXTS   : maximum tesserae number
c     MXCENT : maximum atoms number


      IF (NCOMP.EQ.6) THEN
 
c     unpertubed  electronic potentials on tessera 
         KJ1AO = 1
c     perturbed electronic potential
         KINT = KJ1AO + NNBASX
c     D-1 inverse matrix Q=D-1*V
         KINVMT = KINT + NNBASX*NCOMP
         KLAST = KINVMT + NTSIRR*NTS
         LWRK = LWORK - KLAST + 1  
 
         CALL DZERO(VEC2,NNBASX*NCOMP)
         CALL DZERO(POTPER,NTSIRR*NCOMP)
         CALL DZERO(QSEPER,NTSIRR*NCOMP)
         CALL DZERO(WORK(KINVMT),NTSIRR*NTS)

      ELSE
         CALL QUIT("Number of components not recognized in MAGPCMSEC")
      END IF


c     read D-1 inverse matrix from file

      CALL GPOPEN(LUPCMD,'PCMDATA','UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &     IDUMMY,.FALSE.)
      REWIND(LUPCMD)
      READ(LUPCMD) (WORK(KINVMT+I-1), I = 1, NTSIRR*NTS)
      CALL GPCLOSE(LUPCMD,'KEEP')
 

c     save actual dipole origin
      XI = DIPORG(1)
      YI = DIPORG(2)
      ZI = DIPORG(3)

c     loop over tesserae
      DO  ITS = 1, NTSIRR
         
c     use dipole origin to pass tessera position 
c     to integration program 

         DIPORG(1) = XTSCOR(ITS)
         DIPORG(2) = YTSCOR(ITS)
         DIPORG(3) = ZTSCOR(ITS)

         NREAD = 2
         NTESP = 1
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         KPATOM = 0
         L=1

c     Get the WORK(KINT) pertubed electronic potentials on tessera

         CALL GET1IN(WORK(KINT),'PCMB2SL',NREAD,WORK(KLAST),LWRK,
     &        LABINT,INTREP,INTADR,L,TOFILE,KPATOM,TRIMAT,DUMMY,
     &        .FALSE.,DUMMY,INTPRI)      
         
c     1)Contribution due multiplication of pert potential times 
c     two times nuclear unperturbed charges QSN
c     plus electronic unp charges QSE 

         QTOTUN=-2*QSN(ITS)-QSE(ITS)

         DO ICOMP = 1, NCOMP

            CALL DAXPY(NNBASX,.5*qtotun,
     &           WORK(KINT+(ICOMP - 1)*NNBASX),1,VEC2(1+(ICOMP-1
     &           )*NNBASX),1)


c     Contract the six perturbed potentials with density matrix to
c     get total pertubed electronic potential on the tesserae  

            POTPER(ITS+(ICOMP-1)*NTSIRR)=DDOT(NNBASX,DENMAT,1,
     &           WORK(KINT+(ICOMP - 1)*NNBASX),1)

         END DO

      END DO


c     Get the total electronic perturbed charges Qe[per] on
c     the tesserae

      DO ICOMP=1,NCOMP

         CALL V2QPER(WORK(KINVMT),POTPER(1+(ICOMP-1)*NTSIRR),
     &        QSEPER(1+(ICOMP-1)*NTSIRR))

      END DO

c     2)Contribution due to Ve[unp](mu,nu) *  Qe[per] 

c     open file containing unpertubed potential
      CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',
     &     'SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUPROP)

      DO ITS=1,NTSIRR

c     Get the uncontracted unpertubed potential

         CALL REAPCM('J1-PCMIN','PCMFCK  ',LUPROP,WORK(KJ1AO),NNBASX)

         DO ICOMP=1,6
               
c     perform Ve[unp]*Qe[per] 

            CALL DAXPY(NNBASX,.5*QSEPER(ITS+(ICOMP-1)*NTSIRR),
     &           WORK(KJ1AO),1,VEC2(1+(ICOMP-1
     &           )*NNBASX),1)

         END DO

      END DO


      IF (IPRINT.GT.5) THEN
         CALL AROUND(
     &        'Second order solvent contributions in MAGPCM')
         WRITE (LUPRI,'(2X,A)') 'XX coordinate'
         CALL OUTPAK(VEC2(1),NBAST,1,LUPRI)
         WRITE (LUPRI,'(2X,A)') 'XY coordinate'
         CALL OUTPAK(VEC2(1+NNBASX),NBAST,1,LUPRI)
         WRITE (LUPRI,'(2X,A)') 'XZ coordinate'
         CALL OUTPAK(VEC2(1+2*NNBASX),NBAST,1,LUPRI)
         WRITE (LUPRI,'(2X,A)') 'YY coordinate'
         CALL OUTPAK(VEC2(1+3*NNBASX),NBAST,1,LUPRI)
         WRITE (LUPRI,'(2X,A)') 'YZ coordinate'
         CALL OUTPAK(VEC2(1+4*NNBASX),NBAST,1,LUPRI)
         WRITE (LUPRI,'(2X,A)') 'ZZ coordinate'
         CALL OUTPAK(VEC2(1+5*NNBASX),NBAST,1,LUPRI)
      END IF


c     restoring dipole origin
      DIPORG(1) = XI
      DIPORG(2) = YI
      DIPORG(3) = ZI


      CALL GPCLOSE(LUPROP,'KEEP')
         RETURN
         END


C/* Deck V2QPER */
      SUBROUTINE V2QPER(CINVMT,TJ2PER,QSPER)
C
C     Transform a given electronic potential on the tesserae 
c     into iduced electronic charges
c     Based on V2Q 
c     See ...
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <symmet.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION CINVMT(NTS*NTSIRR), TJ2PER(NTS),QSPER(NTS)
      LOGICAL QCHECK
#include <pcm.h>
#include <pcmlog.h>

      SQTNOP = DBLE(MAXREP + 1)
      DO ISYM = 0, MAXREP

         ISTART = ISYM * NTSIRR ** 2 + 1
         JSTART = ISYM * NTSIRR + 1

c     multiply D-1*TJ2PER to get the charges

         CALL DGEMV('N',NTSIRR,NTSIRR,-1.0D0,CINVMT(ISTART),NTSIRR,
     $        TJ2PER(JSTART),1,D0,QSPER(JSTART),1)

      ENDDO

c     multiply the area a(ITS) factor to the charges

      DO ITS = 1, NTSIRR

         QSPER(ITS) = QSPER(ITS)*AS(ITS)/SQTNOP

      END DO
      
      RETURN
      END
