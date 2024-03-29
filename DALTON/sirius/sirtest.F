!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
!#define noSPECIAL_LINTST
C
C  /* Deck lintst */
      SUBROUTINE LINTST(CMO,CREF,EACTVN,GORB,DV,PV,FC,FV,FCAC,H2AC,
     *                  INDXCI,WRK,LFREE)
C
C     24-Oct-1989 pj
C
C     Calculate full L matrix explicitly for test purposes.
C
#include "implicit.h"
      DIMENSION CMO(*),CREF(*),GORB(*),DV(*),PV(*),FC(*),FV(*)
      DIMENSION FCAC(*),H2AC(*),INDXCI(*),WRK(LFREE)
C
C Used from common blocks:
C   INFVAR : NCONF,NWOPT,JWOPSY
C
#include "maxorb.h"
#include "priunit.h"
#include "infvar.h"
#include "infpri.h"
C
C
      PARAMETER ( DTOL = 1.0D-8 )
      PARAMETER ( D1 = 1.0D0 , D0 = 0.0D0 )
C
      CALL QENTER('LINTST')
      CALL HEADER('flag(61): construct MC Hessian explicitly',-1)
      KL2   = 1
      KBVEC = KL2    + NVAR*NVAR
      KWRK1 = KBVEC  + NVAR
      LWRK1 = LFREE  - KWRK1
C
      CALL DZERO(WRK(KBVEC),NVAR)

#ifdef SPECIAL_LINTST
      ! modify to include the two columns you want details about
      icol1 = 3
      icol2 = 14
      call dzero(WRK(KL2),NVAR*NVAR)
      ipri4 = 999
      ipri6 = 999
      iprsir = 999
      do i = 1,npflag
         p6flag(i) = .true.
      end do
      do i = 1,npflag
         p4flag(i) = .true.
      end do
      call dzero(cref,nconf)
      cref(1) = 1.0d0 ! perhaps modify
#endif

      DO 100 I = 1,NVAR
#ifdef SPECIAL_LINTST
      if (i .ne. icol1 .and. i .ne. icol2) cycle
#endif
         IF (I.LE.NCONF) THEN
            NCSIM = 1
            NOSIM = 0
            IOFF  = I
         ELSE
            NCSIM = 0
            NOSIM = 1
            IOFF  = I - NCONF
         ENDIF
         WRK(KBVEC-1+IOFF) = D1
         IF (( NCSIM.GT.0 ).AND.( JWOPSY.EQ.1)) THEN
            IF (NCONF .EQ. 1) THEN
               CALL DZERO(WRK(KL2+(I-1)*NVAR),NVAR)
               GO TO 100
            ELSE
               CALL DAXPY(NCONF,-CREF(I),CREF,1,WRK(KBVEC),1)
            END IF
         END IF
         IF (IPRI6.GT.90) THEN
            IF(NOSIM.GT.0)  THEN
               KDIM = NWOPT
               WRITE(LUPRI,'(/A)')' LINTST ORBITAL TRIAL VECTOR'
            END IF
            IF(NCSIM.GT.0) THEN
               KDIM = NCONF
               WRITE(LUPRI,'(/A)')' LINTST CONFIGURATION TRIAL VECTOR'
            END IF
            CALL OUTPUT(WRK(KBVEC),1,KDIM,1,1,KDIM,1,1,LUPRI)
         END IF
         CALL LINTRN(NCSIM,NOSIM,WRK(KBVEC),WRK(KBVEC),CMO,CREF,EACTVN,
     *               GORB,DV,PV,FC,FV,FCAC,H2AC,
     *               INDXCI,WRK(KWRK1),1,LWRK1)
C        CALL LINTRN(NCSIM,NOSIM,BCVECS,BOVECS,
C    *               CMO,CREF,EACTVN,GORB,DV,PV,
C    *               FC,FV,FCAC,H2AC, INDXCI,WRK,KFRSAV,LFRSAV)
C
C PROJECT OUT RERERENCE STATE COMPONENTS FROM LINEAR TRANSFORMED
C L2 VECTOR
C
         IF ((NCONF.GE.1).AND.( JWOPSY.EQ.1)) THEN
            XL2OVL = DDOT(NCONF,CREF,1,WRK(KWRK1),1)
            CALL DAXPY(NCONF,-XL2OVL,CREF,1,WRK(KWRK1),1)
         END IF
         CALL DCOPY(NVAR,WRK(KWRK1),1,WRK(KL2+(I-1)*NVAR),1)
         IF (( NCSIM.GT.0 ).AND.( JWOPSY.EQ.1)) THEN
            CALL DZERO(WRK(KBVEC),NCONF)
         ELSE
            WRK(KBVEC-1+IOFF) = D0
         END IF
 100  CONTINUE
      IF (IPRI6 .GT. 10 .OR. NVAR.LE.30) THEN
         WRITE(LUPRI,'(A,I8)')' E[2] MATRIX : DIMENSION',NVAR
         CALL OUTPUT(WRK(KL2),1,NVAR,1,NVAR,NVAR,NVAR,-1,LUPRI)
      END IF
C
C CHECK LARGEST DEVIATION BETWEEN L2(I,J) AND L2(J,I) FOR
C THE VARIOUS COMPONENTS
C
      COMX   = D0
      IC     = 0
      JC     = 0
      ORMX   = D0
      IO     = 0
      JO     = 0
      ORCOMX = D0
      IOC    = 0
      JOC    = 0
      DO 150 I = 1,NVAR
#ifdef SPECIAL_LINTST
      if (i .ne. icol1 .and. i .ne. icol2) cycle
#endif
         DO 160 J=1,I-1
#ifdef SPECIAL_LINTST
      if (j .ne. icol1 .and. j .ne. icol2) cycle
#endif
            IF (I.LE.NCONF) THEN
               XL2IJ = WRK(KL2-1+(J-1)*NVAR+I)
               XL2JI = WRK(KL2-1+(I-1)*NVAR+J)
               CODV = ABS(XL2IJ-XL2JI)
               IF (CODV.GT.COMX) THEN
                  IC = I
                  JC = J
                  COMX = CODV
               END IF
            ELSE IF (J.LE.NCONF) THEN
               XL2IJ = WRK(KL2-1+(J-1)*NVAR+I)
               XL2JI = WRK(KL2-1+(I-1)*NVAR+J)
               ORCODV = ABS(XL2IJ-XL2JI)
               IF (ORCODV.GT.ORCOMX) THEN
                  IOC = I
                  JOC = J
                  ORCOMX = ORCODV
               END IF
            ELSE
               XL2IJ = WRK(KL2-1+(J-1)*NVAR+I)
               XL2JI = WRK(KL2-1+(I-1)*NVAR+J)
               ORCODV = ABS(XL2IJ-XL2JI)
               IF (ORCODV.GT.ORCOMX) THEN
                  IOC = I
                  JOC = J
                  ORCOMX = ORCODV
               END IF
            END IF
 160     CONTINUE
 150  CONTINUE
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' L2(I,J)-conf block : MAXIMUM DEVIATION, IC='
     *,IC,' JC=',JC,' DEV=',COMX
      IF ( COMX.GT.DTOL)
     *WRITE(LUPRI,'(/2(/A,I5,A,I5,A,1P,G16.8))')
     *' IC=',JC,'  JC=',IC,' L2(IC,JC)',WRK(KL2-1+(IC-1)*NVAR+JC),
     *' IC=',IC,'  JC=',JC,' L2(IC,JC)',WRK(KL2-1+(JC-1)*NVAR+IC)
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' L2(I,J)-orbital block : MAXIMUM DEVIATION, IO='
     *,IO,' JO=',JO,' DEV=',ORMX
      IF ( ORMX.GT.DTOL)
     *WRITE(LUPRI,'(/2(/A,I5,A,I5,A,1P,G16.8))')
     *' IO=',JO,'  JO=',IO,' L2(IO,JO)',WRK(KL2-1+(IO-1)*NVAR+JO),
     *' IO=',IO,'  JO=',JO,' L2(IO,JO)',WRK(KL2-1+(JO-1)*NVAR+IO)
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' L2(I,J) conf-orbital block : MAXIMUM DEVIATION, IOC='
     *,IOC,' JOC=',JOC,' DEV=',ORCOMX
      IF ( ORCOMX.GT.DTOL)
     *WRITE(LUPRI,'(/2(/A,I5,A,I5,A,1P,G16.8))')
     *'IOC=',JOC,' JOC=',IOC,' L2(IOC,JOC)',WRK(KL2-1+(IOC-1)*NVAR+JOC),
     *'IOC=',IOC,' JOC=',JOC,' L2(IOC,JOC)',WRK(KL2-1+(JOC-1)*NVAR+IOC)
C
C     End of LINTST.
C
      CALL QUIT('End of LINTST (flag(61))')
      CALL QEXIT('LINTST')
      RETURN
      END
