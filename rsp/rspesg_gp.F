!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
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
      SUBROUTINE D0FCDIAG(D0DIAG,FCDIAG,FC)
#include <implicit.h>
#include <priunit.h> 
#include <infrsp.h>
#include <inforb.h>

      DIMENSION D0DIAG(*),FCDIAG(*),FC(*)

      DO 300 ISYM=1, NSYM  
        IOFF=IIORB(ISYM)+1
        
        DO 100 I=1, NISH(ISYM)
           D0DIAG(IORB(ISYM)+I) =2.0D0
 100    CONTINUE


        DO 200 I=1,NORB(ISYM)
          FCDIAG(IORB(ISYM)+I)=FC(IOFF)
          IOFF = IOFF + I + 1
 200   CONTINUE
 300  CONTINUE

      IF (IPRRSP.GE.10) THEN    
        CALL PRINTVEC('FCDIAG:        ',1,NORBT,FCDIAG)
        CALL PRINTVEC('D0DIAG:        ',1,NORBT,D0DIAG)
 
      END IF 
      RETURN
      END 
C
C END OF D0FCDIAG
C

      SUBROUTINE GETH1MO(H1MO,CMO,WRK,LWRK)
#include <implicit.h>
#include <priunit.h> 
#include <infrsp.h>
#include <inforb.h>

      DIMENSION H1MO(NORBT,NORBT),CMO(*),WRK(*)

C
C     Calls sirh1 to get the one electron integrals 
C     in AO basis (h1AO) and transforms the triangular 
C     output to the full symmetrical matrix NORBT*NORBT
C

      KH1AO  = 1
      KH1AOT = KH1AO + N2BASX

      KWRK1 = KH1AOT + NNORBX+1
      LWRK1 = LWRK - KWRK1 

      CALL SIRH1(WRK(KH1AOT),WRK(KWRK1),LWRK1)

      K    = KH1AOT  
      KOFF = KH1AO

      DO 300 ISYM=1,NSYM
        DO 200 I=1,NBAS(ISYM)
          DO 100 J=1,I
             WRK(KOFF+(I-1)*NBAST+J-1)=WRK(K)
             WRK(KOFF+(J-1)*NBAST+I-1)=WRK(K)
             K=K+1
 100      CONTINUE
 200    CONTINUE
        KOFF = KOFF + NBAS(ISYM)*NBAST + NBAS(ISYM)
 300  CONTINUE

      IF (IPRRSP .GE. 10) THEN
        CALL PRINTMAT('H1AO            ',1,NBAST,WRK(KH1AO))
      END IF

C     Transform the integrals to MO basis

      CALL AO2MO(WRK(KH1AO),H1MO,CMO,WRK(KWRK1),LWRK1)

      RETURN
      END

      SUBROUTINE GETESG_DENMAT_FOCMAT(D,F,NBAST,NNBAST)
#include <implicit.h>
#include "rspprp.h"
#include <esg.h>
#include <dummy.h>
#include <priunit.h>

      DIMENSION D(NNBAST), F(NNBAST)

      WRITE (LUPRI,'(A,/)') 
     &        'Reading density and fock matrices in 1-el. part'

      CALL FLSHFO(LUPRI)

      READ (LUESG) D
      READ (LUESG) F

      CALL GPCLOSE(LUESG,'KEEP')

      RETURN
      END

      SUBROUTINE GETESG_DENMAT(DMAT,NDMAT,WRK,LWRK)
#include <implicit.h>
#include <mxcent.h>
#include "rspprp.h"
#include <esg.h>
#include <dummy.h>
#include <inforb.h>
#include <abainf.h>
#include <priunit.h>

      DIMENSION DMAT(NBAST,NBAST,NDMAT),WRK(*)
      PARAMETER (MXDMAT=50)
      DIMENSION ISYMDM(MXDMAT)

      WRITE(LUPRI,'(A,/)') 'Reading density matrices in 2-el. part'

      KDSO = 1
      KWRK = KDSO + N2BASX
  
      DO 100 I=1,NDMAT
        ISYMDM(I) = 0
 100  CONTINUE
      ISYMDM(3) = ISYME - 1
      ISYMDM(4) = ISYME - 1 

      IF (KWRK.GT.LWRK) CALL STOPIT('GETESG_DENMAT',' ',KWRK,LWRK)

      CALL GPOPEN(LUESG2,'ESG_DMAT',' ',' ',' ',IDUMMY,.FALSE.)
      REWIND (LUESG2)

      DO IDMAT=2,NDMAT 
         CALL READT(LUESG2,N2BASX,WRK)
         CALL DSOTAO(WRK,DMAT(1,1,IDMAT),NBAST,ISYMDM(IDMAT),IPRESG)
      END DO

      CALL GPCLOSE(LUESG2,'KEEP')

      IF ( IPRESG .GE. 10 ) THEN
        WRITE(LUPRI,*) 'Reading density matrices in 2-el. part'
        CALL PRINTMAT('D0AO:          ',1,NBAST,DMAT(1,1,1))
        CALL PRINTMAT('D2AO:          ',1,NBAST,DMAT(1,1,2))
        CALL PRINTMAT('DXAO:          ',1,NBAST,DMAT(1,1,3))
        CALL PRINTMAT('DXSAO:         ',1,NBAST,DMAT(1,1,4))
      END IF

      RETURN
      END

      SUBROUTINE INVE2VEC(CMO,UDV,PV,FC,FV,FCAC,
     & H2AC,NSIM,INPVECS,OUTVECS,XINDX,WRK,LWRK)
C
C CALCULATES INV( E2 ) * VECTOR 
C
#include <implicit.h>
#include <dummy.h>
#include <codata.h>
#include <priunit.h>
#include <infopt.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include "rspprp.h"
#include <infpp.h>
#include <inflr.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <inftap.h>
#include <infsop.h>

      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
      DOUBLE PRECISION INPVECS(*),OUTVECS(*)
      CHARACTER*8 BLANK
C
      PARAMETER ( MAXSIM = 15, BLANK = '        ' )

      PARAMETER ( D2 = 2.0D0, D100 = 100.0D0, D0 = 0.0D0)
      PARAMETER ( D2R3 = (2.0D0/3.0D0))
C
C     space allocation for reduced E(2) and reduced S(2)
C
      KREDE  = 1
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MAXRM
      KEIVEC = KRESID + MAXRM
      KREDGD = KEIVEC + MAXRM*MAXRM
      KWRK1  = KREDGD + MAXRM*MAXRM
      LWRK1  = LWRK + 1 - KWRK1

      CALL DZERO(WRK,KWRK1)

      IF (IPRPP .GT. 2 .OR. LWRK1 .LT. 2*KZYVAR) THEN
         WRITE(LUPRI,'(/A)') ' >>> IN INVE2VEC:'
         WRITE(LUPRI,*)' THCPP,MAXRM ',THCPP,MAXRM
         WRITE(LUPRI,*)' KSYMOP,NGPPP(KSYMOP) ',KSYMOP,NGPPP(KSYMOP)
         WRITE(LUPRI,*)' LWRK ,LWRK1 ',LWRK,LWRK1
      END IF
C
C ALLOCATE WORK SPACE FOR EIGENVECTORS AND TRANSITION MOMENTS
C
      LNEED = 100 + 2*KZYVAR
C
C MAXIMUM NUMBER OF SIMULTANEOUS SOLUTION VECTORS
C
C      NSIM = MIN(KEXCNV, MAXSIM, (LWRK1-LNEED)/KZYVAR )
      IF (IPRPP .GT. 2 .OR. NSIM .LE. 0) THEN
         LWRK2 = KWRK1 + LNEED + KZYVAR
C        ... need at least space for one KZYVAR (NSIM = 1)
         WRITE (LUPRI,*) ' KEXCNV,NSIM,LWRK2 ',KEXCNV,NSIM,LWRK2
         IF (NSIM.LE.0) CALL ERRWRK('RSPPP work space',-LWRK2,LWRK)
      END IF
    
      KWRK2  = KWRK1 + NSIM*KZYVAR

      LWRK2  = LWRK   - KWRK2

      THCRSP = THCPP
      IPRRSP = IPRPP
      MAXIT  = MAXITP
      
      CALL DZERO(WRK(KEIVAL),MAXRM)
      CALL DZERO(WRK(KEIVEC),MAXRM)


      DO 1100 ISIM=1,NSIM
        CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,.TRUE.,BLANK,
     *            BLANK,INPVECS((ISIM-1)*KZYVAR+1),WRK(KREDGD),
     *            WRK(KREDE),WRK(KREDS),
     *            WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),
     *            XINDX,WRK(KWRK1),LWRK1)
        WRITE(LUPRI,*) 'Finished RSPCTL for ISIM=', ISIM
        CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     &            OUTVECS((ISIM-1)*KZYVAR+1), 
     &            WRK(KWRK1),NSIM,ISIM-1)
 1100  CONTINUE

      RETURN
      END
C
C *** END OF INVE2VEC
C
   

      SUBROUTINE FOLD(AIN,AOUT,N)
C
C  This subroutine folds the square matrix AIN into the triangular 
C   matrix AOUT 
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION AIN(N,N), AOUT(*)

      IJ = 0
      DO I=1,N
       DO J=1,I-1
        IJ = IJ + 1
        AOUT(IJ) = AIN(I,J) + AIN(J,I)
       END DO
       IJ = IJ + 1
       AOUT(IJ) = AIN(I,I)
      END DO 

      RETURN
      END 

      SUBROUTINE FMAT2VEC(NSIM,ZYVEC,ZYMAT)
#include <implicit.h>
#include <priunit.h>
      DIMENSION ZYVEC(KZYWOP,*), ZYMAT(NORBT,NORBT,*)
C
#include <maxorb.h>
#include <maxash.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      PARAMETER ( D1 = 1.0D0 )
      CALL DZERO(ZYVEC,NSIM*KZYVAR)
      DO 100 ISIM=1,NSIM
         DO 200 IG=1,KZWOPT
            I=JWOP(1,IG)
            J=JWOP(2,IG)
            ZYVEC(IG,ISIM)=ZYMAT(J,I,ISIM)-ZYMAT(I,J,ISIM)
            ZYVEC(IG+KZWOPT,ISIM)=-ZYVEC(IG,ISIM)
  200    CONTINUE
  100 CONTINUE

      RETURN
      END

      SUBROUTINE SYMMETRIZE(MINP,MOUT,NSIM,NORBT)
#include <implicit.h>
#include <priunit.h> 

      DIMENSION MINP(NORBT,NORBT,*),MOUT(NORBT,NORBT,*)
      DOUBLE PRECISION MINP, MOUT

      DO 300 ISIM=1,NSIM
        DO 200 I=1,NORBT
          DO 100 J=1,I
            MOUT(I,J,ISIM) = 0.5D0*(MINP(I,J,ISIM)
     &                             +MINP(J,I,ISIM))
            MOUT(J,I,ISIM) = 0.5D0*(MINP(I,J,ISIM)
     &                             +MINP(J,I,ISIM))
 100      CONTINUE      
 200    CONTINUE      
 300  CONTINUE      

      RETURN
      END

C
C  Old AO2MO for inputs without the symmetry
C
C      SUBROUTINE AO2MO(MAO,MMO,CMO,NORBT,WRK,LWRK)
C#include <implicit.h>
C#include <priunit.h> 
C
C
C     Purpose: transforms a two-index matrix from 
C     AO to MO basis (using WRK as scratch )
C
C     MMO = CMO^T * MAO * CMO 
C
C
C      DIMENSION MAO(*),MMO(*),CMO(*),WRK(*)
C
C      CALL DGEMM('T','N',NORBT,NORBT,NORBT,1.D0,
C     &                   CMO,NORBT,MAO,NORBT,0.D0,
C     &                   WRK,NORBT)
C
C      CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
C     &                   WRK,NORBT,CMO,NORBT,0.D0,
C     &                   MMO,NORBT)
C
C      RETURN
C      END
C
C
C     End of AO2MO_OLD
C

      SUBROUTINE AO2MO(AAO,BMO,CMO,WRK,LWRK)
#include <implicit.h>
#include <priunit.h> 
#include <inforb.h> 

C
C     Purpose: transforms a two-index matrix from 
C     AO to MO basis (using WRK as scratch )
C
C     BMO = CMO^T * AAO * CMO 
C
C  CMO is in rectangular symmetry blocks, 
C  AAO, BMO are full without any symmetry reduction 
C

      DIMENSION AAO(*),BMO(*),CMO(*),WRK(*)

      CALL DZERO(WRK,N2ORBX)
      CALL DZERO(BMO,N2ORBX)

C      CALL PRINTMAT('AO2MO: beg.    ',1,NORBT,AAO)

      ICOFF = 1 
      DO 100 ISYM=1,NSYM
        IF ( NBAS(ISYM) .GT. 0 ) THEN
C          IOFF =  IORB(ISYM)*NORBT+IORB(ISYM) +1
          IOFF =  IORB(ISYM) + 1
          CALL DGEMM('T','N',
C     &      NORB(ISYM),NORB(ISYM),NBAS(ISYM),1.D0,
     &      NORB(ISYM),NORBT,NBAS(ISYM),1.D0,
     &      CMO(ICOFF),NBAS(ISYM),
     &      AAO(IOFF),NORBT,1.0D0,
     &      WRK(IOFF),NORBT)
          ICOFF = ICOFF + N2ORB(ISYM)
        END IF
 100  CONTINUE

      ICOFF = 1
      DO 200 ISYM=1,NSYM 
        IF ( NBAS(ISYM) .GT. 0 ) THEN
C          IOFF =  IORB(ISYM)*NORBT+IORB(ISYM)+1
          IOFF =  IORB(ISYM)*NORBT+1
          CALL DGEMM('N','N',
C     &           NORB(ISYM),NORB(ISYM),NBAS(ISYM),1.0D0,
     &           NORBT,NORB(ISYM),NBAS(ISYM),1.0D0,
     &           WRK(IOFF),NORBT,
     &           CMO(ICOFF),NBAS(ISYM),1.0D0,
     &           BMO(IOFF),NORBT)
          ICOFF = ICOFF + N2ORB(ISYM)
        END IF
 200  CONTINUE

C      CALL PRINTMAT('AO2MO: final.  ',1,NORBT,BMO)
      
      RETURN
      END

C
C     End of AO2MO
C
      SUBROUTINE MO2AO(AMO,AAO,CMO,WRK,LWRK)
#include <implicit.h>
#include <priunit.h> 
#include <inforb.h>

C
C     Purpose: transforms a two-index matrix from 
C     AO to MO basis (using WRK as scratch )
C
C     AAO = CMO * AMO * CMO^T
C
C  AMO, AAO are full matrix without a symmetry reduction
C
C We make a loop in which we multiply CMO block 
C of certain symmetry with the MO matrix
C
      DIMENSION AAO(*),AMO(*),CMO(*),WRK(*)
      INTEGER COFF

      CALL DZERO(WRK,NORBT**2)
      CALL DZERO(AAO,NORBT**2)

      DO ISYM=1,NSYM
        IF ( NBAS(ISYM) .GT. 0 ) THEN
          CALL DGEMM('N','N',
     &          NBAS(ISYM),NORBT,NORB(ISYM),1.D0,
     &          CMO(ICMO(ISYM) + 1),NBAS(ISYM),
     &          AMO(IORB(ISYM) + 1),NORBT,0.D0, 
     &          WRK(IBAS(ISYM)+1),NBAST)
        END IF
      END DO

      DO ISYM=1,NSYM 
        IF ( NBAS(ISYM) .GT. 0 ) THEN
          CALL DGEMM('N','T',
     &             NORBT,NBAS(ISYM),NORB(ISYM),1.D0,
     &             WRK(IORB(ISYM)*NBAST+1),NBAST,
     &             CMO(ICMO(ISYM) + 1),NBAS(ISYM),0.D0,
     &             AAO(IBAS(ISYM)*NBAST+1),NBAST)
        END IF 
      END DO

      RETURN
      END

C
C This version was working for non-symmetry calc
C
C
C      SUBROUTINE MO2AO_OLD(MMO,MAO,CMO,NORBT,WRK,LWRK)
C#include <implicit.h>
C#include <priunit.h> 
C
C
C     Purpose: transforms a two-index matrix from 
C     AO to MO basis (using WRK as scratch )
C
C     MAO = CMO * MMO * CMO^T
C
C
C      DIMENSION MAO(*),MMO(*),CMO(*),WRK(*)
C      DOUBLE PRECISION MAO,MMO,CMO,WRK
C
C      CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
C     &                   CMO,NORBT,MMO,NORBT,0.D0,
C     &                   WRK,NORBT)
C
C      CALL DGEMM('N','T',NORBT,NORBT,NORBT,1.D0,
C     &                   WRK,NORBT,CMO,NORBT,0.D0,
C     &                   MAO,NORBT)
C
C      RETURN
C      END

C
C     End of MO2AO
C

      SUBROUTINE ESG_PVAL(DODFT,HFXFAC,PVAL,DMAT,A,B,C,D,F)
#include <implicit.h>
#include <inforb.h>
      INTEGER A, B, C, D
      PARAMETER (DP25 = 0.25 D00, DP5 = 0.5 D00,
     &           D1 = 1.0 D00, D2 = 2.0 D00, ZERADD = 1.D-15,
     &           D0 = 0.0 D00, DP125 = 0.125D00, D3=3.0D0, D4=4.0D0)
      DIMENSION PVAL(*),DMAT(NBAST,NBAST,*)
      DOUBLE PRECISION PVAL,DMAT
      LOGICAL DODFT

      IF (.NOT. DODFT) HFXFAC = 1.0D0
      PVAL(1) = D2*F*(DMAT(A,B,2)*DMAT(D,C,1)+DMAT(A,B,1)*DMAT(D,C,2)
     & - DP25*HFXFAC*(DMAT(C,A,2)*DMAT(D,B,1)+DMAT(D,A,2)*DMAT(C,B,1)+
     &         DMAT(C,A,1)*DMAT(D,B,2)+DMAT(D,A,1)*DMAT(C,B,2))
     & + ( DMAT(A,B,4)*DMAT(D,C,4) 
     & - DP125*HFXFAC*(DMAT(A,D,3)*DMAT(B,C,3)+DMAT(D,A,3)*DMAT(C,B,3)+
     &            DMAT(A,C,3)*DMAT(B,D,3)+DMAT(C,A,3)*DMAT(D,B,3))))

      RETURN
      END

C----------------------------------------
C     TEMPORARY SUBROUTINES FOR DEBUGING: 
C----------------------------------------
      SUBROUTINE PRINTVEC( TEXT, NSIM, KZYVAR, V )   
#include <implicit.h>
#include <priunit.h> 

      DIMENSION V(*)
      CHARACTER*15 TEXT

      WRITE (LUPRI,'(A)') TEXT

      DO 200 I=1,NSIM 
        WRITE( LUPRI, * ) 'VECTOR NUMBER : ', I
        DO 100 J=1,KZYVAR
          WRITE( LUPRI, 300 ) I, J,
     &           V( (I-1)*KZYVAR + J )
 100    CONTINUE
 200  CONTINUE

 300  FORMAT(3X,' V',I1,'(',I2,')=',F12.6)

      RETURN
      END
C
C *** END OF PRINTVEC
C
      SUBROUTINE PRINTVEC2( TEXT, NSIM, KZVAR, V )   
#include <implicit.h>
#include <priunit.h> 

      DIMENSION V(*)
      CHARACTER*15 TEXT

      WRITE (LUPRI,'(A)') TEXT

      DO ISIM=1,NSIM
         IF (NSIM.GT.1) WRITE( LUPRI, * ) 'MATRIX NUMBER : ', ISIM
         CALL OUTPUT(V,1,KZVAR,1,2,KZVAR,2,1,LUPRI)
      END DO

c      DO 200 I=1,NSIM 
c        WRITE( LUPRI, * ) 'VECTOR NUMBER : ', I
c        DO 100 J=1,KZYVAR
c          WRITE( LUPRI, 300 ) I, J,
c     &           V( (I-1)*KZYVAR + J )
c 100    CONTINUE
c 200  CONTINUE

c 300  FORMAT(3X,' V',I1,'(',I2,')=',F12.6)

      RETURN
      END
C
C *** END OF PRINTVEC
C


      SUBROUTINE PRINTMAT( TEXT, NSIM, NORBT, A )   
#include <implicit.h>
#include <priunit.h> 


      DIMENSION A(NORBT,NORBT,*)
      CHARACTER*15 TEXT 

      WRITE (LUPRI,'(A)') TEXT
      DO ISIM=1,NSIM
         IF (NSIM.GT.1) WRITE( LUPRI, * ) 'MATRIX NUMBER : ', ISIM
         CALL OUTPUT(A,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END DO
c      DO 200 ISIM=1,NSIM 
c        WRITE( LUPRI, * ) 'MATRIX NUMBER : ', ISIM
c          DO 100 I=1,NORBT
c            WRITE( LUPRI, 400 )
c     &             (A(I,J,ISIM),J=1,NORBT)
c 100      CONTINUE
c 200  CONTINUE
c
c 400  FORMAT(6(3X,F12.6))

      RETURN
      END
C
C *** END OF PRINTVEC
C
