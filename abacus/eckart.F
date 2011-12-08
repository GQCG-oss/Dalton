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
#include <single.h>
C
C /* Deck ecktrn */
      SUBROUTINE ECKTRN(GRDCAR,HESCAR,DIP0,SUSTOT,GTRAN,QUADT,TMAT,
     &            GTRANT,POLARS,POLDD,DIANQC,SPNTOT,
     &            MAGSUS,MOLGFA,QUADRU,SHIELD,SPINRO,POLAR,ALFA,NQCC,
     &            SPNSPN,NCART,NCRTOT,MXFR,NFRVAL,FRVAL,GEOM,
     &            AMASS,NATTYP,ECK,GEOM2,IPRINT,HESSIA,WORK,LWORK)
C     
C     Define an Eckart frame, and transform properties to this frame.
C     Based on the program ECKART by Juha Vaara:
C Juha Vaara  130398
C last change 310398
C
C     Modified for use in Dalton by Kenneth Ruud, San Diego April 1999
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>      
#include <nuclei.h>
C
      PARAMETER (NUCMAX = 15,NSTMAX = 10,NTRIAL = 1500,
     &     TOLX = 1.0D-12,TOLF = 1.0D-12)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      LOGICAL MAGSUS, MOLGFA, QUADRU, SHIELD, SPINRO, POLAR, ALFA,
     &        NQCC, SPNSPN, PLANAR, LINEAR, HESSIA
C
      DIMENSION GEOM(3,MXCENT), AMASS(MXCENT), NATTYP(MXCENT),
     &          ECK(3,MXCENT), GEOM2(NUCDEP,3), EQCM(3), XN(3)
      DIMENSION CMAT(3,3), EIGVAL(3), ANGMOM(3), TINERT(3,3), 
     &          OMEGA(3,3), BVEC(3), XJMAT(3,3), ECKP(3)
      DIMENSION GRDCAR(NCRTOT), HESCAR(NCART,NCART), DIP0(3), 
     &          SUSTOT(3,3), GTRAN(3,3), QUADT(3,3), TMAT(3,3,MXCENT),
     &          GTRANT(3,3,MXCENT), POLARS(3,3), POLDD(3,3,MXFR),
     &          DIANQC(3,3,MXCENT), SPNTOT(MXCOOR,MXCOOR), 
     &          FRVAL(NFRVAL), WORK(LWORK)
C
#include <orgcom.h>
#include <cbiwlk.h>
C
C     Define instantaneous geometry
C     
      CALL CMMASS(GEOM,AMASS,NATTYP,WORK,IPRINT)
      IF (IPRINT .GE. 4) THEN
         CALL HEADER('Instantaneous geometry:',0)
         DO I = 1, NUCDEP
            WRITE(LUPRI,'(I2,3F15.10)') I,GEOM(1,I),GEOM(2,I),
     &                                    GEOM(3,I)
         END DO
      END IF
C
C     Center of mass for equilibrium geometry
C
      IF (IPRINT .GE. 4) THEN
         CALL HEADER('Equilibrium geometry',0)
         DO I = 1, NUCDEP
            WRITE(LUPRI,'(I2,3F15.8)') I,ECKGEO(1,I),ECKGEO(2,I),
     &                                   ECKGEO(3,I)
         END DO
      END IF
      TOTMAS = D0
      CALL DZERO(EQCM,3)
      DO I = 1, NUCDEP
         RMASS = AMASS(I)
         TOTMAS = TOTMAS + RMASS
         DO J = 1, 3
            EQCM(J) = EQCM(J) + RMASS*ECKGEO(J,I)
         END DO
      END DO
      DO J = 1, 3
         EQCM(J) = EQCM(J)/TOTMAS
      END DO
C
C     We now start to build the Eckart frame
C
      IF (IPRINT .GE. 2) THEN
         WRITE(LUPRI,'(A,3F15.6)') ' equilibrium CM  : ',EQCM(1),
     &                                               EQCM(2),EQCM(3)
         WRITE(LUPRI,'(A,3F15.6)') ' instantaneous CM: ',CMXYZ(1),
     &                                               CMXYZ(2),CMXYZ(3)
      END IF
C
C construct B vector and J matrix for this isotopomer
C
      CALL DZERO(BVEC,3)
      CALL DZERO(XJMAT,9)
      DO 16 J = 1, NUCDEP
         GEOM(1,J)   = GEOM(1,J) - CMXYZ(1)
         GEOM(2,J)   = GEOM(2,J) - CMXYZ(2)
         GEOM(3,J)   = GEOM(3,J) - CMXYZ(3)
         ECK(1,J) = ECKGEO(1,J) - EQCM(1)
         ECK(2,J) = ECKGEO(2,J) - EQCM(2)
         ECK(3,J) = ECKGEO(3,J) - EQCM(3)
         CALL CROSSP(ECK(1,J),GEOM(1,J),ANGMOM)
         BVEC(1) = BVEC(1) + AMASS(J)*ANGMOM(1)
         BVEC(2) = BVEC(2) + AMASS(J)*ANGMOM(2)
         BVEC(3) = BVEC(3) + AMASS(J)*ANGMOM(3)
         TMP = DDOT(3,ECK(1,J),1,GEOM(1,J),1)
         XJMAT(1,1) = XJMAT(1,1) + AMASS(J)*(TMP -ECK(1,J)*GEOM(1,J))
         XJMAT(2,2) = XJMAT(2,2) + AMASS(J)*(TMP -ECK(2,J)*GEOM(2,J))
         XJMAT(3,3) = XJMAT(3,3) + AMASS(J)*(TMP -ECK(3,J)*GEOM(3,J))
         XJMAT(1,2) = XJMAT(1,2) - AMASS(J)*ECK(1,J)*GEOM(2,J)
         XJMAT(1,3) = XJMAT(1,3) - AMASS(J)*ECK(1,J)*GEOM(3,J)
         XJMAT(2,1) = XJMAT(2,1) - AMASS(J)*ECK(2,J)*GEOM(1,J)
         XJMAT(2,3) = XJMAT(2,3) - AMASS(J)*ECK(2,J)*GEOM(3,J)
         XJMAT(3,1) = XJMAT(3,1) - AMASS(J)*ECK(3,J)*GEOM(1,J)
         XJMAT(3,2) = XJMAT(3,2) - AMASS(J)*ECK(3,J)*GEOM(2,J)
 16   CONTINUE
C
C     Mismatch in the ordering of atom and coordinate for the geometry in
C     CMMASS and WLKDIN
C
      DO I = 1, NUCDEP
         DO J = 1, 3
            GEOM2(I,J) = GEOM(J,I)
         END DO
      END DO
C
C     Moment of inertia for instantaneous geometry
C
      CALL WLKDIN(GEOM2,AMASS,NUCDEP,ANGMOM,TINERT,OMEGA,
     &            EIGVAL,CMAT,.TRUE.,PLANAR,LINEAR)
C
C     We will need rotational g tensors and spin-rotation constants in 
C     space-fixed xyz axis system
C
      IF (MOLGFA .AND. .NOT. LINEAR) THEN
         CALL DGEMM('N','N',3,3,3,1.D0,
     &              CMAT,3,
     &              GTRAN,3,0.D0,
     &              OMEGA,3)
         CALL DGEMM('N','T',3,3,3,1.D0,
     &              OMEGA,3,
     &              CMAT,3,0.D0,
     &              GTRAN,3)
         IF (IPRINT .GE. 4) THEN
            CALL HEADER('Rotational g tensor in Cartesian '//
     &                  'input frame',0)
            CALL OUTPUT(GTRAN,1,3,1,3,3,3,1,LUPRI)
         END IF
      END IF
C
      IF (SPINRO) THEN
         DO I = 1, NUCDEP
            CALL DGEMM('N','N',3,3,3,1.D0,
     &                 CMAT,3,
     &                 GTRANT(1,1,I),3,0.D0,
     &                 OMEGA,3)
            CALL DGEMM('N','T',3,3,3,1.D0,
     &                 OMEGA,3,
     &                 CMAT,3,0.D0,
     &                 GTRANT(1,1,I),3)
         END DO
         IF (IPRINT .GE. 4) THEN
            DO I = 1, NUCDEP
               CALL HEADER('Spin-rotation tensor in Cartesian '//
     &                     'input frame',0)
               CALL OUTPUT(GTRANT(1,1,I),1,3,1,3,3,3,1,LUPRI)
            END DO
         END IF
      END IF
C
      IF (IPRINT .GE. 3) THEN
         CALL HEADER('Equilibrium geometry in the CM frame:',0)
         DO 17 J = 1, NUCDEP
            WRITE(LUPRI,'(I2,3F15.10)') J,ECK(1,J),ECK(2,J),
     &                                    ECK(3,J)
 17      CONTINUE
         CALL HEADER('Instantaneous geometry in the CM frame:',0)
         DO 18 J = 1, NUCDEP
            WRITE(LUPRI,'(I2,3F15.10)') J,GEOM(1,J),GEOM(2,J),GEOM(3,J)
 18      CONTINUE
      END IF
      IF (IPRINT .GE. 2) THEN
         CALL HEADER('B-vector: ',0)
         WRITE(LUPRI,'(A,F14.6,2F15.6)') 'B',(BVEC(I), I=1,3)
         CALL HEADER('J-matrix:',0)
         CALL OUTPUT(XJMAT,1,3,1,3,3,3,1,LUPRI)
      END IF
C
      ECKP(1) = 1.0D0
      ECKP(2) = 1.0D0
      ECKP(3) = 1.0D0
C
      CALL ECK_MNEWT(NTRIAL,ECKP,BVEC,XJMAT,TOLX,TOLF)
      XN(1) = SIN(ECKP(1))*COS(ECKP(2))
      XN(2) = SIN(ECKP(1))*SIN(ECKP(2))
      XN(3) = COS(ECKP(1))
      AN = ECKP(3)
C     
      IF (IPRINT .GE. 2) THEN
         CALL HEADER('Parameters of the Eckart frame:',0)
         WRITE(LUPRI,'(A,F13.6,2F15.6)') 'NV',(XN(I), I = 1, 3)
         WRITE(LUPRI,'(A,F13.6)') 'AN',AN
      END IF
C
C
      IF (ABS(XN(1)).LT.1.0D-7 .AND. ABS(XN(2)).LT.1.0D-7 .AND.
     &        ABS(D1-XN(3)).LT.1.0D-7 ) THEN
         WRITE(LUPRI,'(A)') ' (rotation in the xy-plane only;'
     &           //' zeroing a redundant parameter)'
         ECKP(2) = D0
         ECKP(1) = D0
      ENDIF
C
      CALL CLCGLD(ECKP,CMAT)
      IF (IPRINT .GE. 2) THEN
         CALL HEADER('C-matrix:',0)
         CALL OUTPUT(CMAT,1,3,1,3,3,3,1,LUPRI)
      END IF
C
      IF (IPRINT .GE. 4) 
     &     CALL HEADER('Instantaneous geometry in the Eckart frame:',0)
      DO 20 J = 1, NUCDEP
         CALL DGEMM('N','N',3,1,3,1.D0,
     &              CMAT,3,
     &              GEOM(1,J),3,0.D0,
     &              ANGMOM,3)
         CALL DCOPY(3,ANGMOM,1,GEOM(1,J),1)
      IF (IPRINT .GE. 4) WRITE(LUPRI,'(A,I2,3F15.10)') 
     &        'IN',J, GEOM(1,J), GEOM(2,J), GEOM(3,J)
 20   CONTINUE
C
C     Transform all properties to Eckart frame
C
      CALL DGEMM('N','N',3,1,3,1.D0,
     &           CMAT,3,
     &           DIP0,3,0.D0,
     &           ANGMOM,3)
      CALL DCOPY(3,ANGMOM,1,DIP0,1)
      IF (IPRINT .GE. 2) THEN
         CALL HEADER('Dipole moment in Eckart frame',0)
         CALL OUTPUT(DIP0,1,3,1,1,3,1,1,LUPRI)
      END IF
C
      IF (MAGSUS) THEN
         CALL DGEMM('N','N',3,3,3,1.D0,
     &              CMAT,3,
     &              SUSTOT,3,0.D0,
     &              OMEGA,3)
         CALL DGEMM('N','T',3,3,3,1.D0,
     &              OMEGA,3,
     &              CMAT,3,0.D0,
     &              SUSTOT,3)
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Magnetizability tensor '
     &           //'in the Eckart frame:',0)
            CALL OUTPUT(SUSTOT,1,3,1,3,3,3,1,LUPRI)
         ENDIF
      END IF
C
      IF (MOLGFA) THEN
         CALL DGEMM('N','N',3,3,3,1.D0,
     &              CMAT,3,
     &              GTRAN,3,0.D0,
     &              OMEGA,3)
         CALL DGEMM('N','T',3,3,3,1.D0,
     &              OMEGA,3,
     &              CMAT,3,0.D0,
     &              GTRAN,3)
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Molecular rotational g-tensor '
     &           //'in the Eckart frame:',0)
            CALL OUTPUT(GTRAN,1,3,1,3,3,3,1,LUPRI)
         ENDIF
      END IF
C
      IF (QUADRU) THEN
         CALL DGEMM('N','N',3,3,3,1.D0,
     &              CMAT,3,
     &              QUADT,3,0.D0,
     &              OMEGA,3)
         CALL DGEMM('N','T',3,3,3,1.D0,
     &              OMEGA,3,
     &              CMAT,3,0.D0,
     &              QUADT,3)
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Quadrupole moment tensor '
     &           //'in the Eckart frame:',0)
            CALL OUTPUT(QUADT,1,3,1,3,3,3,1,LUPRI)
         ENDIF
      END IF
C
C
      IF (HESSIA) THEN
C
C     Gradient
C
         DO I = 1, 3*NUCDEP, 3
            CALL DGEMM('N','N',3,3,3,1.0D0,
     &                 CMAT,3,
     &           GRDCAR(I),3,0.0D0,
     &           OMEGA,3)
            CALL DGEMM('N','T',3,3,3,1.0D0,
     &                 OMEGA,3,
     &                 CMAT,3,0.0D0,
     &                 GRDCAR(I),3)
         END DO
C
C     Hessian
C
         DO I = 1, NCART, 3
         DO J = 1, NCART, 3
            CALL DGEMM('N','N',3,3,3,1.0D0,
     &                 CMAT,3,
     &           HESCAR(I,J),3,0.0D0,
     &           OMEGA,3)
            CALL DGEMM('N','T',3,3,3,1.0D0,
     &                 OMEGA,3,
     &                 CMAT,3,0.0D0,
     &                 HESCAR(I,J),3)
         END DO
         END DO
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Molecular Hessian in the Eckart frame:',0)
            CALL OUTPAK(HESCAR,NCART,1,LUPRI)
         END IF
      END IF
C
      IF (SHIELD) THEN
         DO I =1 , NUCDEP
            CALL DGEMM('N','N',3,3,3,1.D0,
     &                 CMAT,3,
     &                 TMAT(1,1,I),3,0.D0,
     &                 OMEGA,3)
            CALL DGEMM('N','T',3,3,3,1.D0,
     &                 OMEGA,3,
     &                 CMAT,3,0.D0,
     &                 TMAT(1,1,I),3)
         END DO
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Shielding tensors '
     &              //'in the Eckart frame:',0)
            DO I = 1, NUCDEP
               CALL OUTPUT(TMAT(1,1,I),1,3,1,3,3,3,1,LUPRI)
            END DO
         ENDIF
      END IF
C
      IF (SPINRO) THEN
         DO I =1 , NUCDEP
            CALL DGEMM('N','N',3,3,3,1.D0,
     &                 CMAT,3,
     &                 GTRANT(1,1,I),3,0.D0,
     &                 OMEGA,3)
            CALL DGEMM('N','T',3,3,3,1.D0,
     &                 OMEGA,3,
     &                 CMAT,3,0.D0,
     &                 GTRANT(1,1,I),3)
         END DO
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Nuclear spin-rotation tensors '
     &              //'in the Eckart frame:',0)
            DO I = 1, NUCDEP
               CALL OUTPUT(GTRANT(1,1,I),1,3,1,3,3,3,1,LUPRI)
            END DO
         ENDIF
      END IF
C
      IF (POLAR) THEN
         CALL DGEMM('N','N',3,3,3,1.D0,
     &              CMAT,3,
     &              POLARS,3,0.D0,
     &              OMEGA,3)
         CALL DGEMM('N','T',3,3,3,1.D0,
     &              OMEGA,3,
     &              CMAT,3,0.D0,
     &              POLARS,3)
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Polarzability tensor '
     &           //'in the Eckart frame:',0)
            CALL OUTPUT(POLARS,1,3,1,3,3,3,1,LUPRI)
         ENDIF
      END IF
C
      IF (ALFA) THEN
         DO IFR = 1, NFRVAL
            CALL DGEMM('N','N',3,3,3,1.D0,
     &                 CMAT,3,
     &                 POLDD(1,1,IFR),3,0.D0,
     &                 OMEGA,3)
            CALL DGEMM('N','T',3,3,3,1.D0,
     &                 OMEGA,3,
     &                 CMAT,3,0.D0,
     &                 POLDD(1,1,IFR),3)
         END DO
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Polarizability tensors '
     &           //'in the Eckart frame:',0)
            DO IFR = 1, NFRVAL
               CALL OUTPUT(POLDD(1,1,IFR),1,3,1,3,3,3,1,LUPRI)
            END DO
         ENDIF
      END IF
C
      IF (NQCC) THEN
         DO I = 1, NUCDEP
            CALL DGEMM('N','N',3,3,3,1.D0,
     &                 CMAT,3,
     &                 DIANQC(1,1,I),3,0.D0,
     &                 OMEGA,3)
            CALL DGEMM('N','T',3,3,3,1.D0,
     &                 OMEGA,3,
     &                 CMAT,3,0.D0,
     &                 DIANQC(1,1,I),3)
         END DO
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Nuclear quadrupole moment tensors '
     &           //'in the Eckart frame:',0)
            DO I = 1, NUCDEP
               CALL OUTPUT(DIANQC(1,1,I),1,3,1,3,3,3,1,LUPRI)
            END DO
         ENDIF
      END IF
C
      RETURN
      END
C
      SUBROUTINE CROSSP(A,B,C)
C
C calculates (AX,AY,AZ) x (BX,BY,BZ) = CZ,CY,CZ)
C
#include <implicit.h>
      DIMENSION A(3), B(3), C(3)
C
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
      RETURN
      END
      SUBROUTINE ECK_USRFUN(ECKP,FVEC,FJAC,BX,BY,BZ,
     &                      JXX,JXY,JXZ,JYX,JYY,JYZ,JZX,JZY,JZZ)
C
C the system of equations to be solved for obtaining the Eckart frame
C is specified with this subroutine
C
      IMPLICIT NONE
C
      DOUBLE PRECISION ECKP(3),FVEC(3),FJAC(3,3),SINK,COSK,
     &     SINL,COSL,SINM,COSM,JXX,JXY,JXZ,JYX,JYY,JYZ,JZX,JZY,JZZ,
     &     BX,BY,BZ
C
      SINK    = SIN(ECKP(1))
      COSK    = COS(ECKP(1))
      SINL    = SIN(ECKP(2))
      COSL    = COS(ECKP(2))
      SINM    = SIN(ECKP(3))
      COSM    = COS(ECKP(3))
C
      FVEC(1) = BX*COSM + (1 - COSM)*(JZY*COSK**2 +
     &     JXY*COSK*COSL*SINK + JYY*COSK*SINK*SINL 
     &     -JZZ*COSK*SINK*SINL - JXZ*COSL*SINK**2*SINL 
     &     -JYZ*SINK**2*SINL**2) 
     &     +(JXZ*COSK + JXX*COSL*SINK + JXY*SINK*SINL)*SINM
      FVEC(2) = BY*COSM + (1 - COSM)*(-(JZX*COSK**2) - 
     &     JXX*COSK*COSL*SINK + JZZ*COSK*COSL*SINK 
     &     +JXZ*COSL**2*SINK**2 - JYX*COSK*SINK*SINL + 
     &     JYZ*COSL*SINK**2*SINL) + (JYZ*COSK + 
     &     JYX*COSL*SINK + JYY*SINK*SINL)*SINM
      FVEC(3) = BZ*COSM + (1 - COSM)*(-(JZY*COSK*COSL*SINK) 
     &     - JXY*COSL**2*SINK**2 + JZX*COSK*SINK*SINL + 
     &     JXX*COSL*SINK**2*SINL - JYY*COSL*SINK**2*SINL + 
     &     JYX*SINK**2*SINL**2) + (JZZ*COSK + JZX*COSL*SINK 
     &     +JZY*SINK*SINL)*SINM
C
      FJAC(1,1) = (1 - COSM)*(JXY*COSK**2*COSL - 
     &     2*JZY*COSK*SINK - JXY*COSL*SINK**2 + 
     &     JYY*COSK**2*SINL - JZZ*COSK**2*SINL - 
     &     2*JXZ*COSK*COSL*SINK*SINL - 
     &     JYY*SINK**2*SINL + JZZ*SINK**2*SINL - 
     &     2*JYZ*COSK*SINK*SINL**2) + 
     &     (JXX*COSK*COSL - JXZ*SINK + JXY*COSK*SINL)*SINM
      FJAC(2,1) = (1 - COSM)*(-(JXX*COSK**2*COSL) + 
     &     JZZ*COSK**2*COSL + 2*JZX*COSK*SINK + 
     &     2*JXZ*COSK*COSL**2*SINK + JXX*COSL*SINK**2 
     &     - JZZ*COSL*SINK**2 - 
     &     JYX*COSK**2*SINL + 2*JYZ*COSK*COSL*SINK*SINL 
     &     + JYX*SINK**2*SINL) + 
     &     (JYX*COSK*COSL - JYZ*SINK + JYY*COSK*SINL)*SINM
      FJAC(3,1) = (1 - COSM)*(-(JZY*COSK**2*COSL) - 
     &     2*JXY*COSK*COSL**2*SINK + JZY*COSL*SINK**2 + 
     &     JZX*COSK**2*SINL + 
     &     2*JXX*COSK*COSL*SINK*SINL - 
     &     2*JYY*COSK*COSL*SINK*SINL - 
     &     JZX*SINK**2*SINL + 2*JYX*COSK*SINK*SINL**2) + 
     &     (JZX*COSK*COSL - JZZ*SINK + JZY*COSK*SINL)*SINM
      FJAC(1,2) = (1 - COSM)*(JYY*COSK*COSL*SINK - 
     &     JZZ*COSK*COSL*SINK - JXZ*COSL**2*SINK**2 - 
     &     JXY*COSK*SINK*SINL - 
     &     2*JYZ*COSL*SINK**2*SINL + JXZ*SINK**2*SINL**2) + 
     &     (JXY*COSL*SINK - JXX*SINK*SINL)*SINM
      FJAC(2,2) = (1 - COSM)*(-(JYX*COSK*COSL*SINK) 
     &     + JYZ*COSL**2*SINK**2 + JXX*COSK*SINK*SINL - 
     &     JZZ*COSK*SINK*SINL - 2*JXZ*COSL*SINK**2*SINL 
     &     - JYZ*SINK**2*SINL**2) + 
     &     (JYY*COSL*SINK - JYX*SINK*SINL)*SINM
      FJAC(3,2) = (1 - COSM)*(JZX*COSK*COSL*SINK + 
     &     JXX*COSL**2*SINK**2 - JYY*COSL**2*SINK**2 + 
     &     JZY*COSK*SINK*SINL + 
     &     2*JXY*COSL*SINK**2*SINL + 2*JYX*COSL*SINK**2*SINL - 
     &     JXX*SINK**2*SINL**2 + JYY*SINK**2*SINL**2) + 
     &     (JZY*COSL*SINK - JZX*SINK*SINL)*SINM
      FJAC(1,3) =  COSM*(JXZ*COSK + JXX*COSL*SINK + 
     &     JXY*SINK*SINL) - BX*SINM + 
     &     (JZY*COSK**2 + JXY*COSK*COSL*SINK + 
     &     JYY*COSK*SINK*SINL - JZZ*COSK*SINK*SINL - 
     &     JXZ*COSL*SINK**2*SINL - JYZ*SINK**2*SINL**2)*SINM
      FJAC(2,3) = COSM*(JYZ*COSK + JYX*COSL*SINK + 
     &     JYY*SINK*SINL) - BY*SINM + 
     &     (-(JZX*COSK**2) - JXX*COSK*COSL*SINK + 
     &     JZZ*COSK*COSL*SINK + 
     &     JXZ*COSL**2*SINK**2 - JYX*COSK*SINK*SINL 
     &     + JYZ*COSL*SINK**2*SINL)*SINM
      FJAC(3,3) = COSM*(JZZ*COSK + JZX*COSL*SINK + 
     &     JZY*SINK*SINL) - BZ*SINM + 
     &     (-(JZY*COSK*COSL*SINK) - JXY*COSL**2*SINK**2 + 
     &     JZX*COSK*SINK*SINL + 
     &     JXX*COSL*SINK**2*SINL - 
     &     JYY*COSL*SINK**2*SINL + JYX*SINK**2*SINL**2)*SINM
C
      RETURN
      END
C
C
      SUBROUTINE CLCGLD(ECKP,C)
C
C calculates the matrix required for transformation
C from the input to the Eckart frame 
C The 'Goldstein' way (Classical Mechanics, p.165) and
C thus avoiding unnecessary interchange of axes!
C
      IMPLICIT NONE
C
      DOUBLE PRECISION ECKP(3),C(3,3),NX,NY,NZ,SINPH,M1COS,COSPH
C
      NX    = SIN( ECKP(1) ) * COS( ECKP(2) )
      NY    = SIN( ECKP(1) ) * SIN( ECKP(2) )
      NZ    = COS( ECKP(1) )
      COSPH = COS( ECKP(3) )
      M1COS = 1.D0 - COSPH
      SINPH = SIN( ECKP(3) )
C
      C(1,1) = COSPH + NX * NX * M1COS 
      C(2,2) = COSPH + NY * NY * M1COS 
      C(3,3) = COSPH + NZ * NZ * M1COS 
C
      C(1,2) =         NX * NY * M1COS - NZ * SINPH
      C(1,3) =         NX * NZ * M1COS + NY * SINPH
      C(2,1) =         NY * NX * M1COS + NZ * SINPH
      C(2,3) =         NY * NZ * M1COS - NX * SINPH
      C(3,1) =         NZ * NX * M1COS - NY * SINPH
      C(3,2) =         NZ * NY * M1COS + NX * SINPH
C
      RETURN
      END
C
C
      SUBROUTINE ECK_mnewt(ntrial,x,bvec,xjmat,tolx,tolf)
C
      INTEGER ntrial
      DOUBLE PRECISION tolf,tolx,x(3)
      INTEGER i,k,indx(3)
      DOUBLE PRECISION errf,errx,fjac(3,3),fvec(3),p(3),
     &     bvec(3), xjmat(3,3)
      do k=1,ntrial               
         call ECK_usrfun(x,fvec,fjac,bvec(1),bvec(2),bvec(3),
     &               xjmat(1,1),xjmat(1,2),xjmat(1,3),xjmat(2,1),
     &               xjmat(2,2),xjmat(2,3),xjmat(3,1),xjmat(3,2),
     &               xjmat(3,3))
         errf =0.0d0
         do i=1,3
            errf=errf+abs(fvec(i))
         enddo
         if(errf.le.tolf)return
         do i=1,3                       
            p(i)=-fvec(i)
         enddo 
         call ECK_lucdmpHJ(fjac,indx)    
         call ECK_lubksb(fjac,indx,p)  
         errx=0.0d0
         do i=1,3
            errx=errx+abs(p(i))
            x(i)=x(i)+p(i)
         enddo
         if(errx.le.tolx)return
      enddo
      return
      END
C
C
      SUBROUTINE ECK_lucdmpHJ(a,indx)
C
C Rev. Mar. 2006 hjaaj so it can treat singular
C  matrices (which occur for linear molecules)
C
      INTEGER indx(3)
      DOUBLE PRECISION a(3,3),TINY
      PARAMETER (TINY=1.0d-20) 
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(3)
#include <priunit.h>
      call output(a,1,3,1,3,3,3,1,LUPRI)
      do i=1,3
         aamax=0.0d0
         do j=1,3
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         enddo
CHJ      if (aamax.eq.0.0d0) call quit('singular matrix in ludcmp')
CHJ      vv(i)=1.d0/aamax                 
         if (aamax.le.TINY) then
            a(i,i)=TINY
            vv(i)=1.d0/TINY
         else
            vv(i)=1.d0/aamax                 
         end if
CHJ end
      enddo
      do j=1,3
         do i=1,j-1                
            sum=a(i,j)
            do k=1,i-1
               sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
         enddo 
         aamax=0.0d0
         imax=j
         do i=j,3
            sum=a(i,j)                 
            do k=1,j-1
               sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
            dum=vv(i)*abs(sum)         
            if (dum.ge.aamax) then     
               imax=i
               aamax=dum
            endif
         enddo
         if (j.ne.imax)then             
            do k=1,3
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            enddo
            vv(imax)=vv(j)             
         endif
         indx(j)=imax
         if(a(j,j).eq.0.0d0)a(j,j)=TINY
         if(j.ne.3)then
            dum=1.d0/a(j,j)
            do i=j+1,3
               a(i,j)=a(i,j)*dum
            enddo
         endif
      enddo                        
      return
      END
C
C
      SUBROUTINE ECK_lubksb(a,indx,b)
C
C
      INTEGER indx(3)
      DOUBLE PRECISION a(3,3),b(3)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0               
      do i=1,3
         ll=indx(i)      
         sum=b(ll)        
         b(ll)=b(i)
         if (ii.ne.0)then
            do j=ii,i-1
               sum=sum-a(i,j)*b(j)
            enddo
         else if (sum .ne. 0.0d0) then
            ii=i           
         endif               
         b(i)=sum
      enddo
      do i=3,1,-1         
         sum=b(i)
         do j=i+1,3
            sum=sum-a(i,j)*b(j)
         enddo
         b(i)=sum/a(i,i)     
      enddo 
      return                   
      END
