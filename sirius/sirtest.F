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
C$Id: sirtest.F,v 1.2 2001-05-31 06:38:37 vebjornb Exp $
C
#include <single.h>
C  /* Deck lintst */
      SUBROUTINE LINTST(CMO,CREF,EACTVN,GORB,DV,PV,FC,FV,FCAC,H2AC,
     *                  INDXCI,WRK,LFREE)
C
C     24-Oct-1989 pj
C
C     Calculate full L matrix explicitly for test purposes.
C
#include <implicit.h>
      DIMENSION CMO(*),CREF(*),GORB(*),DV(*),PV(*),FC(*),FV(*)
      DIMENSION FCAC(*),H2AC(*),INDXCI(*),WRK(LFREE)
C
C Used from common blocks:
C   INFVAR : NCONF,NWOPT,JWOPSY
C
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include <infpri.h>
C
C
      PARAMETER ( DTOL = 1.0D-8 )
      PARAMETER ( D1 = 1.0D0 , D0 = 0.0D0 )
C
      KL2   = 1
      KBVEC = KL2    + NVAR*NVAR
      KWRK1 = KBVEC  + NVAR
      LWRK1 = LFREE  - KWRK1
C
      CALL DZERO(WRK(KBVEC),NVAR)
      DO 100 I = 1,NVAR
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
               WRITE(LUPRI,'(/A)')' ORBITAL TRIAL VECTOR'
            END IF
            IF(NCSIM.GT.0) THEN
               KDIM = NCONF
               WRITE(LUPRI,'(/A)')' CONFIGURATION TRIAL VECTOR'
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
      IF (IPRI6 .GT. 10) THEN
         WRITE(LUPRI,'(A,I8)')' L(2) MATRIX : DIMENSION ',NVAR
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
         DO 160 J=1,I-1
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
     *' IC=',IC,' JC=',JC,' L2(IC,JC)',WRK(KL2-1+(IC-1)*NVAR+JC),
     *' IC=',JC,' JC=',IC,' L2(IC,JC)',WRK(KL2-1+(JC-1)*NVAR+IC)
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' L2(I,J)-orbital block : MAXIMUM DEVIATION, IO='
     *,IO,' JO=',JO,' DEV=',ORMX
      IF ( ORMX.GT.DTOL)
     *WRITE(LUPRI,'(/2(/A,I5,A,I5,A,1P,G16.8))')
     *' IO=',IO,' JO=',JO,' L2(IO,JO)',WRK(KL2-1+(IO-1)*NVAR+JO),
     *' IO=',JO,' JO=',IO,' L2(IO,JO)',WRK(KL2-1+(JO-1)*NVAR+IO)
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' L2(I,J) conf-orbital block : MAXIMUM DEVIATION, IOC='
     *,IOC,' JOC=',JOC,' DEV=',ORCOMX
      IF ( ORCOMX.GT.DTOL)
     *WRITE(LUPRI,'(/2(/A,I5,A,I5,A,1P,G16.8))')
     *'IOC=',IOC,'JOC=',JOC,' L2(IOC,JOC)',WRK(KL2-1+(IOC-1)*NVAR+JOC),
     *'IOC=',JOC,'JOC=',IOC,' L2(IOC,JOC)',WRK(KL2-1+(JOC-1)*NVAR+IOC)
C
C     End of LINTST.
C
      RETURN
      END
