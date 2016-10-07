C************************************************************************************
C  /* Deck calc_mep */
      SUBROUTINE CALC_MEP(CMO,WRK,KFRSAV,LFRSAV)
!
! Jacob Kongsted, calculate Molecular Electrostatic Potential
!                 in (x,y,z) points listed on MEP.INP file;
!                 values are written to MEP.OUT file.

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
#include "gnrinf.h"
#ifdef VAR_MPI
#include "infpar.h"
#include "mpif.h"
#endif
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WRK(LFRSAV)
      DIMENSION CMO(*)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LSAVE
      LOGICAL LOCDEB
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*2 UNITS
      PARAMETER ( D1 = 1.0D0 )
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)
      INTEGER NCENTERS, NLOOP, MLOOP, IPRTYP, IERR

      CALL QENTER('CALC_MEP')

C     Only scf for the moment 

      KFREE  = KFRSAV
      LFREE  = LFRSAV

      LOCDEB = .FALSE.

      LUMEP = -1
      CALL GPOPEN(LUMEP,'MEP.INP','OLD',' ',
     &           'FORMATTED',IDUMMY,.FALSE.)
      REWIND(LUMEP)
      READ(LUMEP,*) NCENTERS
      READ(LUMEP,'(A2)') UNITS
      CALL UPCASE(UNITS)

      IF (UNITS .EQ. 'AA') THEN
        SCAL = XTANG
      ELSE IF (UNITS .EQ. 'AU') THEN
        SCAL = 1.0D0
      ELSE
         CALL QUIT('Unknown units in POTENTIAL.INP')
      ENDIF
   
      CALL MEMGET('REAL',KDV   ,NNASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDENC ,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDENV ,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDENS ,NNBASX,WRK,KFREE,LFREE)

      CALL DZERO(WRK(KDV),NNASHX)
      CALL DZERO(WRK(KDENC),N2BASX)
      CALL DZERO(WRK(KDENV),N2BASX)
      CALL DZERO(WRK(KDENS),NNBASX)

      IF (NASHT .EQ. 1) THEN
         WRK(KDV) = D1
      ELSE IF (HSROHF) THEN
         DO I = 1, NASHT
            II = I*(I+1)/2
            WRK(KDV+II-1) = D1
         END DO
      ENDIF

      CALL FCKDEN((NISHT.GT.0),(NASHT.GT.0),WRK(KDENC),WRK(KDENV),
     &            CMO,WRK(KDV),WRK(KFREE),LFREE)
      CALL DAXPY(N2BASX,1.0D0,WRK(KDENV),1,WRK(KDENC),1)
      CALL DZERO(WRK(KDENV),N2BASX)
      CALL DGEFSP(NBAST,WRK(KDENC),WRK(KDENV))
      CALL PKSYM1(WRK(KDENV),WRK(KDENS),NBAS,NSYM,1)

      CALL MEMGET('REAL',KXMEP ,NCENTERS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KYMEP ,NCENTERS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KZMEP ,NCENTERS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KMEP  ,NCENTERS,WRK,KFREE,LFREE)

      CALL MEMGET('REAL',KMEFX ,NCENTERS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KMEFY ,NCENTERS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KMEFZ ,NCENTERS,WRK,KFREE,LFREE)

      DO 101 J = 1,NCENTERS
        READ(LUMEP,*) WRK(KXMEP+J-1),WRK(KYMEP+J-1),WRK(KZMEP+J-1)
        WRK(KXMEP+J-1) = WRK(KXMEP+J-1)/SCAL
        WRK(KYMEP+J-1) = WRK(KYMEP+J-1)/SCAL
        WRK(KZMEP+J-1) = WRK(KZMEP+J-1)/SCAL
  101 CONTINUE

      CLOSE(LUMEP)

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

#ifdef VAR_MPI
      IF (NODTOT .GE. 1) THEN
       NLOOP = NCENTERS / (NODTOT + 1)
       IF ((NODTOT + 1) * NLOOP .LT. NCENTERS) THEN
         MLOOP = NCENTERS - ((NODTOT + 1) * NLOOP)
         KXMEPN = KXMEP + ((NODTOT + 1) * NLOOP)
         KYMEPN = KYMEP + ((NODTOT + 1) * NLOOP)
         KZMEPN = KZMEP + ((NODTOT + 1) * NLOOP)
         KMEPN = KMEP + ((NODTOT + 1) * NLOOP)
         KMEFXN = KMEFX + ((NODTOT + 1) * NLOOP)
         KMEFYN = KMEFY + ((NODTOT + 1) * NLOOP)
         KMEFZN = KMEFZ + ((NODTOT + 1) * NLOOP)
       ELSE
         MLOOP = 0
       END IF

       IPRTYP = MEP_WORK
       ISOMETHING = -1
       CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
       CALL MPIXBCAST(ISOMETHING,1,'INTEGER',MASTER)

       CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
       CALL MPIXBCAST(NLOOP,1,'INTEGER',MASTER)
       CALL MPIXBCAST(WRK(KDENS),NNBASX,'DOUBLE',MASTER)

       CALL MPI_SCATTER(WRK(KXMEP),NLOOP,MPI_DOUBLE_PRECISION,
     &                  WRK(KXMEP),NLOOP,MPI_DOUBLE_PRECISION,MASTER,
     &                  MPI_COMM_WORLD,IERR)
       CALL MPI_SCATTER(WRK(KYMEP),NLOOP,MPI_DOUBLE_PRECISION,
     &                  WRK(KYMEP),NLOOP,MPI_DOUBLE_PRECISION,MASTER,
     &                  MPI_COMM_WORLD,IERR)
       CALL MPI_SCATTER(WRK(KZMEP),NLOOP,MPI_DOUBLE_PRECISION,
     &                  WRK(KZMEP),NLOOP,MPI_DOUBLE_PRECISION,MASTER,
     &                  MPI_COMM_WORLD,IERR)

       CALL MEP_LOOP(WRK(KXMEP),WRK(KYMEP),WRK(KZMEP),NLOOP,
     &                 WRK(KDENS),WRK(KFREE),LFREE,
     &                 WRK(KMEP),WRK(KMEFX),WRK(KMEFY),WRK(KMEFZ))

       IF (MLOOP .NE. 0) THEN
        CALL MEP_LOOP(WRK(KXMEPN),WRK(KYMEPN),WRK(KZMEPN),MLOOP,
     &                  WRK(KDENS),WRK(KFREE),LFREE,
     &                  WRK(KMEPN),WRK(KMEFXN),WRK(KMEFYN),WRK(KMEFZN))
       END IF

       CALL MPI_GATHER(WRK(KMEP),NLOOP,MPI_DOUBLE_PRECISION,WRK(KMEP),
     &                 NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                 IERR)
       CALL MPI_GATHER(WRK(KMEFX),NLOOP,MPI_DOUBLE_PRECISION,WRK(KMEFX),
     &                 NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                 IERR)
       CALL MPI_GATHER(WRK(KMEFY),NLOOP,MPI_DOUBLE_PRECISION,WRK(KMEFY),
     &                 NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                 IERR)
       CALL MPI_GATHER(WRK(KMEFZ),NLOOP,MPI_DOUBLE_PRECISION,WRK(KMEFZ),
     &                 NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                 IERR)
      ELSE
#endif
       CALL MEP_LOOP(WRK(KXMEP),WRK(KYMEP),WRK(KZMEP),NCENTERS,
     &                 WRK(KDENS),WRK(KFREE),LFREE,WRK(KMEP),WRK(KMEFX),
     &                 WRK(KMEFY),WRK(KMEFZ))
#ifdef VAR_MPI
      ENDIF
#endif
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      LUMEP = -1
      CALL GPOPEN(LUMEP,'MEP.OUT','NEW',' ',
     &           'FORMATTED',IDUMMY,.FALSE.)
      REWIND(LUMEP)
      WRITE(LUMEP,*) NCENTERS
      WRITE(LUMEP,'(A2)') UNITS

      DO 106 J = 1,NCENTERS
        WRITE(LUMEP,1000) WRK(KXMEP+J-1)*SCAL,WRK(KYMEP+J-1)*SCAL,
     &                    WRK(KZMEP+J-1)*SCAL,WRK(KMEP+J-1),
     &                    WRK(KMEFX+J-1),WRK(KMEFY+J-1),WRK(KMEFZ+J-1)
  106 CONTINUE

      CLOSE(LUMEP)

      CALL MEMREL('CALC_MEP',WRK,KFRSAV,KFRSAV,KFREE,LFREE)

 1000 FORMAT(7(F15.10,2X))

      CALL QEXIT('CALC_MEP')

      END
C************************************************************************************
#ifdef VAR_MPI
C  /* Deck MEP_SLAVE */
      SUBROUTINE MEP_SLAVE(WRK,LWRK)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
#include "gnrinf.h"
#ifdef VAR_MPI
#include "mtags.h"
#include "infpar.h"
#include "mpif.h"
#endif

      DIMENSION WRK(LWRK)
      INTEGER NLOOP

      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NLOOP,1,'INTEGER',MASTER)

      KDENS = 1
      KMEP = KDENS + NNBASX
      KMEFX = KMEP + NLOOP
      KMEFY = KMEFX + NLOOP
      KMEFZ = KMEFY + NLOOP
      KXMEP = KMEFZ + NLOOP
      KYMEP = KXMEP + NLOOP
      KZMEP = KYMEP + NLOOP
      KWRK = KZMEP + NLOOP

      IF (KWRK .GE. LWRK) CALL QUIT('NOT ENOUGH MEMORY')

      CALL DZERO(WRK(KDENS),NNBASX)
      CALL DZERO(WRK(KMEP),NLOOP)
      CALL DZERO(WRK(KMEFX),NLOOP)
      CALL DZERO(WRK(KMEFY),NLOOP)
      CALL DZERO(WRK(KMEFZ),NLOOP)
      CALL DZERO(WRK(KXMEP),NLOOP)
      CALL DZERO(WRK(KYMEP),NLOOP)
      CALL DZERO(WRK(KZMEP),NLOOP)

      CALL MPIXBCAST(WRK(KDENS),NNBASX,'DOUBLE',MASTER)

      CALL MPI_SCATTER(WRK(KXMEP),NLOOP,MPI_DOUBLE_PRECISION,WRK(KXMEP),
     &                 NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                 IERR)
      CALL MPI_SCATTER(WRK(KYMEP),NLOOP,MPI_DOUBLE_PRECISION,WRK(KYMEP),
     &                 NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                 IERR)
      CALL MPI_SCATTER(WRK(KZMEP),NLOOP,MPI_DOUBLE_PRECISION,WRK(KZMEP),
     &                 NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                 IERR)

      CALL MEP_LOOP(WRK(KXMEP),WRK(KYMEP),WRK(KZMEP),NLOOP,
     &                WRK(KDENS),WRK(KWRK),LWRK,
     &                WRK(KMEP),WRK(KMEFX),WRK(KMEFY),WRK(KMEFZ))

      CALL MPI_GATHER(WRK(KMEP),NLOOP,MPI_DOUBLE_PRECISION,WRK(KMEP),
     &                NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                IERR)
      CALL MPI_GATHER(WRK(KMEFX),NLOOP,MPI_DOUBLE_PRECISION,WRK(KMEFX),
     &                NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                IERR)
      CALL MPI_GATHER(WRK(KMEFY),NLOOP,MPI_DOUBLE_PRECISION,WRK(KMEFY),
     &                NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                IERR)
      CALL MPI_GATHER(WRK(KMEFZ),NLOOP,MPI_DOUBLE_PRECISION,WRK(KMEFZ),
     &                NLOOP,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,
     &                IERR)

      END
#endif
C************************************************************************************
C  /* Deck MEP_LOOP */
      SUBROUTINE MEP_LOOP(XMEP,YMEP,ZMEP,NCENTERS,DENS,
     &                      WRK,LWRK,MEP,MEFX,MEFY,MEFZ)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
#include "gnrinf.h"
#ifdef VAR_MPI
#include "mtags.h"
#include "infpar.h"
#include "mpif.h"
#endif

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LSAVE, LOCDEB
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DOUBLE PRECISION XMEP, YMEP, ZMEP, MEP, MEFX, MEFY, MEFZ
      DIMENSION XMEP(NCENTERS),YMEP(NCENTERS),ZMEP(NCENTERS)
      DIMENSION MEP(NCENTERS)
      DIMENSION MEFX(NCENTERS),MEFY(NCENTERS),MEFZ(NCENTERS)
      DIMENSION DENS(NNBASX)
      DIMENSION WRK(LWRK)

      LOCDEB=.FALSE.
      KMAT = 1
      KMATEF = KMAT + NNBASX
      KWRK = KMATEF + 3*NNBASX

      IF (KWRK .GE. LWRK) CALL QUIT('NOT ENOUGH MEMORY')

      CALL DZERO(WRK(KMAT),NNBASX)
      CALL DZERO(WRK(KMATEF),3*NNBASX)

      DO 102 J = 1,NCENTERS
        KPATOM = 0
        NOSIM = 1
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.
        IPRINT = 0
        DIPORG(1) = XMEP(J)
        DIPORG(2) = YMEP(J)
        DIPORG(3) = ZMEP(J)

        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMAT),'NPETES ',NOSIM,WRK(KWRK),
     &              LWRK,LABINT,INTREP,INTADR,J,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.

        IF (LOCDEB)  THEN
          WRITE (LUPRI,'(/A)') 'Pot. energy matrix in CALC_MEP'
          CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)
        ENDIF

        MEP(J) = DDOT(NNBASX,DENS,1,WRK(KMAT),1)

        ECHCHL = 0.0D0
        DO 103 I = 1,NUCIND
          XDIS   = XMEP(J) - CORD(1,I)
          YDIS   = YMEP(J) - CORD(2,I)
          ZDIS   = ZMEP(J) - CORD(3,I)
          DIST2  = XDIS**2+YDIS**2+ZDIS**2
          DIST   = SQRT(DIST2)
          ECHCHL = ECHCHL + CHARGE(I)/DIST
  103   CONTINUE

        MEP(J) = MEP(J) + ECHCHL

        NOSIM = 3

        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMATEF),'NEFIELD',NOSIM,WRK(KWRK),
     &              LWRK,LABINT,INTREP,INTADR,J,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.

        IF (LOCDEB)  THEN
           WRITE (LUPRI,'(/A)') ' Ex matrix:'
           CALL OUTPAK(WRK(KMATEF),NBAST,1,LUPRI)

           WRITE (LUPRI,'(/A)') ' Ey matrix:'
           CALL OUTPAK(WRK(KMATEF+NNBASX),NBAST,1,LUPRI)

           WRITE (LUPRI,'(/A)') ' Ez matrix:'
           CALL OUTPAK(WRK(KMATEF+2*NNBASX),NBAST,1,LUPRI)
        END IF

        MEFX(J) = DDOT(NNBASX,DENS,1,WRK(KMATEF),1)
        MEFY(J) = DDOT(NNBASX,DENS,1,WRK(KMATEF+NNBASX),1)
        MEFZ(J) = DDOT(NNBASX,DENS,1,WRK(KMATEF+2*NNBASX),1)

        ENUCX = 0.0D0
        ENUCY = 0.0D0
        ENUCZ = 0.0D0
        DO 104 I = 1,NUCIND
          XDIS   = XMEP(J) - CORD(1,I)
          YDIS   = YMEP(J) - CORD(2,I)
          ZDIS   = ZMEP(J) - CORD(3,I)
          DIST2  = XDIS**2+YDIS**2+ZDIS**2
          DIST   = SQRT(DIST2)
          DIST3  = DIST**3
          ENUCX  = ENUCX + CHARGE(I)*XDIS/DIST3
          ENUCY  = ENUCY + CHARGE(I)*YDIS/DIST3
          ENUCZ  = ENUCZ + CHARGE(I)*ZDIS/DIST3
  104   CONTINUE

        MEFX(J) = MEFX(J) + ENUCX
        MEFY(J) = MEFY(J) + ENUCY
        MEFZ(J) = MEFZ(J) + ENUCZ

  102 CONTINUE

      END

