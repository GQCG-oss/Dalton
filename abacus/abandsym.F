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
C     /*Deck fcsini*/
      SUBROUTINE FCSINI
C     ************************************************************************
C     *** This routine initializes the variables needed for a more general ***
C     ***                          symmetry treatment.                     ***
C     ************************************************************************
#include <implicit.h>
C
#include <fcsym.h>
C
      IF (FCLASS .EQ. 'C1 ') THEN
C
C     *** Not supposed to do anything symmetry related. ***
C
      ELSE IF ((FCLASS(2:2).EQ.'i') .OR. (FCLASS(2:2).EQ.'s')) THEN
         IF (FCLASS(2:2).EQ.'s') VPLANE = .TRUE.
         IF (FCLASS(2:2).EQ.'i') ICNTR  = .TRUE.
         NCORDR = 1
         NGORDR = 2
         NGVERT = 2
         NCVERT = 2
         NUMELM = 1
         NMTRX  = 1
         NONEDI = 1
         NDEGIR = 0
         N1DIME = 2
         N2DIME = 0
      ELSE
         IF ((FCLASS(1:1).EQ.'S').OR.((FCLASS(1:1).EQ.'D').AND.
     &                                (FCLASS(3:3).EQ.'d'))) THEN
            ROTARE = .TRUE.
            READ (FCLASS(2:2),'(I1)') NCORDR
            NONEDI = 2
            IF (FCLASS(3:3).EQ.'d') THEN
               NCORDR = 2*NCORDR
            END IF
            NDEGIR = (NCORDR-NONEDI)/2
            SEPDEG = .TRUE.
         ELSE
            MROTAX = .TRUE.
            READ (FCLASS(2:2),'(I1)') NCORDR
            NONEDI = 1
            IF (MOD(NCORDR,2).EQ.0) NONEDI = 2
            NDEGIR = (NCORDR-NONEDI)/2
         END IF
C
         IF ((FCLASS(1:1).NE.'D').AND.(FCLASS(3:3) .EQ. ' ')) THEN
            NGORDR = NCORDR
            NGVERT = NONEDI + 4*NDEGIR
            NCVERT = NONEDI +   NDEGIR
            NUMELM = 0
            N1DIME = NONEDI
            N2DIME = NDEGIR
         ELSE IF ((FCLASS(1:1).EQ.'C').AND.(FCLASS(3:3).EQ.'h')) THEN
            NGORDR = 2*NCORDR
            NGVERT = 2*(NONEDI + 4*NDEGIR)
            NCVERT = 2*(NONEDI +   NDEGIR)
            NUMELM = 1
            N1DIME = 2*NONEDI
            N2DIME = 2*NDEGIR
            SEPDEG = .TRUE.
         ELSE
            NGORDR = 2*NCORDR
            NUMELM = 1
            N1DIME = 2*NONEDI
            N2DIME = NDEGIR
            IF ((FCLASS(1:1).EQ.'D').AND.(FCLASS(3:3).NE.'d')) THEN
               NGORDR = NGORDR*2
               NUMELM = 2
               N1DIME = 2*N1DIME
               N2DIME = 2*N2DIME
            END IF
            NGVERT = NGORDR
            NCVERT = N1DIME + N2DIME
         END IF
C
         NMTRX = NUMELM + 1
      END IF
C
      NIREP = N1DIME + N2DIME
C
      IF ((FCLASS(1:1).EQ.'D').AND.(FCLASS(3:3).NE.'d')) ROTAX2 = .TRUE.
      IF  (FCLASS(3:3).EQ.'v')                           VPLANE = .TRUE.
      IF  (FCLASS(3:3).EQ.'h')                           HPLANE = .TRUE.
C
      RETURN
      END
C
C
C     /*Deck grpchr*/
      SUBROUTINE GRPCHR(COOR,SYMCOR,GRIREP,CHRCTR,WORK,ICRIRP,LWORK,
     &                  IPRINT)
C
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0)
#include <nuclei.h>
#include <trkoor.h>
#include <fcsym.h>
C
      DIMENSION COOR(NCOOR), SYMCOR(NCOOR,NCOOR), GRIREP(NGORDR,NGVERT),
     &          CHRCTR(NGORDR,NCVERT), WORK(LWORK), ICRIRP(NCOOR,2)
C
      IF (FCLASS .NE. 'C1 ') THEN
         CALL HEADER('Making symmetry adapted coor in '//FCLASS(1:3),-1)
C
         KSIRRP = 1
         KLAST  = KSIRRP + NCORDR**2
         LWRK   = LWORK  - KLAST
         CALL CYCLIC(WORK(KSIRRP),GRIREP,CHRCTR,WORK(KLAST),LWRK,IPRINT)
C
         KATMTR = 1
         KCORTR = KATMTR +   NUCDEP**2*NMTRX
         KTRMTX = KCORTR + 9         *NMTRX
         KTMPSM = KTRMTX + 9*NUCDEP**2*NGORDR
         KTMPTR = KTMPSM +   NCOOR
         KTMPVC = KTMPTR + 2*NCOOR**2
         KNSTBC = KTMPVC +   NCOOR
         KLAST  = KNSTBC +   NUCDEP
         LWRK   = LWORK - KLAST
         CALL MKSYMC(COOR,SYMCOR,WORK(KATMTR),WORK(KCORTR),WORK(KTRMTX),
     &               WORK(KTMPSM),GRIREP,WORK(KTMPTR),WORK(KTMPVC),
     &               WORK(KLAST),ICRIRP,WORK(KNSTBC),LWRK,IPRINT)
      ELSE
         CALL HEADER('There is only C1-symmetry, using cart. coor',-1)
C
         KDIM = NCOOR**2
         CALL DZERO(SYMCOR,KDIM)
         DO 100 ICOOR = 1, NCOOR
            SYMCOR(ICOOR,ICOOR) = D1
            ICRIRP(ICOOR,1) = 1
            ICRIRP(ICOOR,2) = 0
 100     CONTINUE
         N1DIME = 1
      END IF
C
      RETURN
      END
C
C     /*Deck cyclic*/
      SUBROUTINE CYCLIC(SIRREP,GRIREP,CHRCTR,WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <pi.h>
#include <priunit.h>
      PARAMETER (NDEG = 2)
      PARAMETER (D1 = 1.0D0)
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
C
      DIMENSION SIRREP(NCORDR,NCORDR), GRIREP(NGORDR,NGVERT),
     &          CHRCTR(NGORDR,NCVERT), WORK(LWORK)
#include <fcsym.h>
C
      DO 100 J = 0, NONEDI-1
      DO 100 I = 1, NCORDR
         SIRREP(I,J+1) = DBLE((-1)**(J*(I-1)))
 100  CONTINUE
C
      IORDR  = NONEDI
      DO 200 N = 1, NDEGIR
      DO 200 M = 1, NDEG
         IORDR = IORDR + 1
         DO 300 K = 1, NCORDR
            PHASE = D1
            IF (M .EQ. 2) PHASE = -D1
C
            SIRREP(K,IORDR) = PHASE*DBLE(N*(K-1))
 300     CONTINUE
 200  CONTINUE
C
      CONST = D2*PI/DBLE(NCORDR)
C
      IF ((FCLASS(1:1).EQ.'S').OR.((FCLASS(1:1).EQ.'D').AND.
     &                             (FCLASS(3:3).EQ.'d'))) THEN
         KRDREP = 1
         CALL REDRST(WORK(KRDREP),SIRREP,NDEG,IPRINT)
      ELSE
         KRDREP = 1
         CALL REDTRA(WORK(KRDREP),SIRREP,NDEG,IPRINT)
      END IF
C
      IF ((FCLASS(1:1).EQ.'C').AND.(FCLASS(3:3).EQ.'h')) THEN
         KRDREP = 1
         CALL CHREP(GRIREP,WORK(KRDREP),SIRREP,NDEG,IPRINT)
      ELSE
         KRDREP = 1
         KTMPRD = KRDREP + NCORDR*NDEGIR*NDEG**2
         CALL GENREP(GRIREP,WORK(KRDREP),SIRREP,WORK(KTMPRD),NDEG,
     &               IPRINT)
      END IF
C
C     *** Calculate character table. ***
C
      CALL CALCHR(GRIREP,CHRCTR)
C
C     *** Print section. ***
C
      IF (IPRINT .GT. 7) THEN
         WRITE (LUPRI,'(A)') '                                      '
         WRITE (LUPRI,'(A)') 'The char. for the cyclic "basis" '//FCLASS
         DO J = 1, NONEDI
            WRITE (LUPRI,'(I4,A,10F6.1)')
     &              J, '. irep', (SIRREP(I,J),I=1,NCORDR)
            WRITE (LUPRI,'(A)')  '                   '
         END DO
C
         IORDR = NONEDI
         DO K = 1, NDEGIR
            DO J = 1, NDEG
               IORDR = IORDR + 1
               WRITE (LUPRI,'(I4,A,10F6.1)')  NONEDI+K,
     &                            '. irep', (SIRREP(I,IORDR),I=1,NCORDR)
            END DO
            WRITE (LUPRI,'(A)')  '                   '
         END DO
C
      END IF
C
      CALL HEADER('Group matrix elements for the group '//FCLASS(1:3),
     &                                                              -1)
C
      DO ID1 = 1, N1DIME
         WRITE (LUPRI,'(I2,A)') ID1, '. irep:'
         WRITE (LUPRI,'(24F6.2)') (GRIREP(I,ID1), I = 1, NGORDR)
      END DO
C
      IORDR = N1DIME
      DO ID2 = 1, N2DIME
         JEXT = N1DIME + ID2
         WRITE (LUPRI,'(I2,A)') JEXT, '. irep:'
         DO I = 1, NDEG**2
            IORDR = IORDR + 1
            WRITE (LUPRI,'(24F6.2)')
     &                       (GRIREP(J,IORDR), J = 1, NGORDR)
         END DO
      END DO
      WRITE (LUPRI,'(A)') '                                      '
C
      CALL HEADER('Characters for the group '// FCLASS(1:3),0)
      DO IREP = 1, N1DIME+N2DIME
         WRITE (LUPRI,'(I2,A)') IREP, '. irep:'
         WRITE (LUPRI,'(24F5.1)') (CHRCTR(I,IREP), I = 1, NGORDR)
      END DO
      WRITE (LUPRI,'(A)') '                                      '
C
      RETURN
      END
C
C
C     /*Deck chrep*/
      SUBROUTINE CHREP(GRIREP,REDREP,SIRREP,NDEG,IPRINT)
#include <implicit.h>
      PARAMETER (D1 = 1.0D0)
C
#include <fcsym.h>
      DIMENSION GRIREP(NGORDR,NGVERT), SIRREP(NCORDR,NCORDR),
     &          REDREP(NDEG,NDEG,NCORDR,NDEGIR)
C
      DO 100 J = 1, NONEDI
      DO 100 I = 1, NCORDR
         JNM = 2*(J-1) + 1
         GRIREP(I       ,JNM  ) =  SIRREP(I,J)
         GRIREP(I+NCORDR,JNM  ) =  SIRREP(I,J)
         GRIREP(I       ,JNM+1) =  SIRREP(I,J)
         GRIREP(I+NCORDR,JNM+1) = -SIRREP(I,J)
 100  CONTINUE
C
      DO 200 M     = 1, 2
      DO 200 L     = 1, NDEGIR
      DO 200 IORDR = 1, NCORDR
         PHASE = D1
         IF (M.EQ.2) PHASE = -D1
         IPLACE = 2*NONEDI + 4*((M-1)*NDEGIR + (L-1))
C
         DO 300 J = 1, NDEG
         DO 300 I = 1, NDEG
            IPLACE = IPLACE + 1
            GRIREP(IORDR       ,IPLACE) =       REDREP(I,J,IORDR,L)
            GRIREP(IORDR+NCORDR,IPLACE) = PHASE*REDREP(I,J,IORDR,L)
 300     CONTINUE
 200  CONTINUE
C
      RETURN
      END
C
C
C     /*Deck redtra*/
      SUBROUTINE REDTRA(REDREP,SIRREP,NDEG,IPRINT)
C
#include <implicit.h>
#include <pi.h>
#include <priunit.h>
C
      PARAMETER (D2 = 2.0D0)
#include <fcsym.h>
      DIMENSION REDREP(NDEG,NDEG,NCORDR,NDEGIR), SIRREP(NCORDR,NCORDR)
C
      CONST = D2*PI/DBLE(NCORDR)
C
      DO 100 J = 1, NDEGIR
      DO 100 I = 1, NCORDR
         JEXT = 2*J-1 + NONEDI
         ANGLE = CONST*SIRREP(I,JEXT)
C
         IF (IPRINT .GT. 20) THEN
            WRITE (LUPRI,'(A)') 'Main rotation angle:'
            DEGANG = 180.0D0*ANGLE/PI
            WRITE (LUPRI,'(F10.6,F12.6)') ANGLE, DEGANG
         END IF
C
         REDREP(1,1,I,J) =  COS(ANGLE)
         REDREP(1,2,I,J) = -SIN(ANGLE)
         REDREP(2,1,I,J) =  SIN(ANGLE)
         REDREP(2,2,I,J) =  COS(ANGLE)
 100  CONTINUE
C
      IF (IPRINT .GT. 7) THEN
         WRITE (LUPRI,'(A)') '                            '
         WRITE (LUPRI,'(A)') 'The transformed cyclic group' //
     &                       'degenerate representations'
         DO 200 L = 1, NDEGIR
         DO 200 J = 1, NDEG
         DO 200 I = 1, NDEG
            WRITE (LUPRI,'(10F8.4)') (REDREP(I,J,K,L),K=1,NCORDR)
 200     CONTINUE
      END IF
C
      RETURN
      END
C
C
C     /*Deck redrst*/
      SUBROUTINE REDRST(REDREP,SIRREP,NDEG,IPRINT)
C
#include <implicit.h>
#include <pi.h>
#include <priunit.h>
C
      PARAMETER (D2 = 2.0D0, D4 = 4.0D0)
#include <fcsym.h>
      DIMENSION REDREP(NDEG,NDEG,NCORDR,NDEGIR), SIRREP(NCORDR,NCORDR)
C
      CONST = D2*PI/DBLE(NCORDR)
C
      DO 100 J = 1, NDEGIR
      DO 100 I = 1, NCORDR
         JEXT = 2*J-1 + NONEDI
         ANGLE = CONST*SIRREP(I,JEXT)
C
         IF (IPRINT .GT. 20) THEN
            WRITE (LUPRI,'(A)') 'Main rotation angle:'
            DEGANG = 180.0D0*ANGLE/PI
            WRITE (LUPRI,'(F10.6,F12.6)') ANGLE, DEGANG
         END IF
C
         REDREP(1,1,I,J) =  COS(ANGLE)
         REDREP(1,2,I,J) = -SIN(ANGLE)
         REDREP(2,1,I,J) =  SIN(ANGLE)
         REDREP(2,2,I,J) =  COS(ANGLE)
 100  CONTINUE
C
      IF (IPRINT .GT. 7) THEN
         WRITE (LUPRI,'(A)') '                            '
         WRITE (LUPRI,'(A)') 'The transformed improper rotational' //
     &                       'group degenerate representations'
         DO 200 L = 1, NDEGIR
            WRITE (LUPRI,'(A,i4)') 'Degenerate irep', L
            DO 300 J = 1, NDEG
            DO 300 I = 1, NDEG
               WRITE (LUPRI,'(10F8.4)') (REDREP(I,J,K,L),K=1,NCORDR)
 300        CONTINUE
 200     CONTINUE
      END IF
C
      RETURN
      END
C
C
C
      SUBROUTINE GENREP(GRIREP,REDREP,SIRREP,TMPRED,NDEG,IPRINT)
C
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
#include <fcsym.h>
C
      DIMENSION GRIREP(NGORDR,NGVERT), SIRREP(NCORDR,NCORDR),
     &          REDREP(NDEG,NDEG,NCORDR,NDEGIR),
     &          TMPRED(2*NCORDR,2*NCORDR)
C
      IF (NUMELM .EQ. 0) THEN
         DO 10 J = 1, NONEDI
         DO 10 I = 1, NCORDR
            GRIREP(I,J) = SIRREP(I,J)
 10      CONTINUE
C
         IPLACE = NONEDI
         DO 20 L = 1, NDEGIR
         DO 20 K = 1, NDEG
         DO 20 J = 1, NDEG
            IPLACE = IPLACE + 1
            DO 30 I = 1, NCORDR
               GRIREP(I,IPLACE) = REDREP(J,K,I,L)
 30         CONTINUE
 20      CONTINUE
      ELSE
         DO 100 L = 1, 2
         DO 100 K = 1, 2
            ISTART = (K-1)*NCORDR
            DO 200 J = 1, NONEDI
            DO 200 I = 1, NCORDR
               PHASE = D1
               IF ((K.EQ.2).AND.(L.EQ.2)) PHASE = -D1
C
               INM = ISTART + I
               JNM = NONEDI*(L-1) + J
               GRIREP(INM,JNM) = PHASE*SIRREP(I,J)
 200        CONTINUE
 100     CONTINUE
C
         DO 300 M = 1, 2
            ISTART = (M-1)*NCORDR
            IPLACE =  2*NONEDI
            DO 400 L = 1, NDEGIR
            DO 400 K = 1, NDEG
            DO 400 J = 1, NDEG
               IPLACE = IPLACE + 1
               PHASE  = D1
               IF ((M.EQ.2).AND.(K.EQ.2)) PHASE = -D1
C
               DO 500 I = 1, NCORDR
                  INM = ISTART + I
                  GRIREP(INM,IPLACE) = PHASE*REDREP(J,K,I,L)
 500           CONTINUE
 400        CONTINUE
 300     CONTINUE
      END IF
C
      IF (NUMELM .EQ. 2) THEN

         DO 600 J = 1, 2*NCORDR
         DO 600 I = 1, 2*NCORDR
            TMPRED(I,J) = GRIREP(I,J)
 600     CONTINUE
C
         DO 700 L = 1, 2
         DO 700 K = 1, 2
C
            PHASE = D1
            IF ((K.EQ.2).AND.(L.EQ.2)) PHASE = -D1
C
            DO 800 J = 1, 2*NONEDI
            DO 800 I = 1, 2*NCORDR
               INM = (L-1)*2*NCORDR+I
               JNM = (K-1)*2*NONEDI+J
               GRIREP(INM,JNM) = PHASE*TMPRED(I,J)
 800        CONTINUE
 700     CONTINUE
C
         ISTART = 2*NONEDI
         DO 900 M = 1, 2
         DO 900 L = 1, 2
            JVAL = ISTART
            DO 1100 K = 1, NDEGIR
            DO 1100 J = 1, NDEG**2
               JVAL = JVAL + 1
               PHASE = D1
               IF ((L.EQ.2).AND.(M.EQ.2)) PHASE = -D1
C
               DO 1100 I = 1, 2*NCORDR
                  INM = (M-1)*2*NCORDR + I
                  JNM = (L-1)*4*NDEGIR + 2*NONEDI + JVAL
                  GRIREP(INM,JNM) = PHASE*TMPRED(I,JVAL)
 1100          CONTINUE
C
 900     CONTINUE
C
      END IF
      RETURN
      END
C
C
C     /*Deck calchr*/
      SUBROUTINE CALCHR(GRIREP,CHRCTR)
#include <implicit.h>
#include <priunit.h>
C
#include <fcsym.h>
C
      DIMENSION GRIREP(NGORDR,NGVERT), CHRCTR(NGORDR,NCVERT)
C
      DO 100 IREP  = 1, N1DIME
      DO 100 IORDR = 1, NGORDR
         CHRCTR(IORDR,IREP) = GRIREP(IORDR,IREP)
 100  CONTINUE
C
      IGPLAC = N1DIME + 1
      DO 200 I = 1, N2DIME
         IREP = N1DIME + I
         DO 300 IORDR = 1, NGORDR
            CHRCTR(IORDR,IREP) = GRIREP(IORDR,IGPLAC  )
     &                         + GRIREP(IORDR,IGPLAC+3)
 300     CONTINUE
         IGPLAC = IGPLAC + 4
 200  CONTINUE
C
      RETURN
      END
C
C
C     /*Deck atmtra*/
      SUBROUTINE TRMTRX(COOR,ATMTRA,CORTRA,IPRINT)
C
#include <implicit.h>
#include <pi.h>
#include <mxcent.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0)
      PARAMETER (DTHR = 1.0D-6)
#include <nuclei.h>
#include <trkoor.h>
#include <fcsym.h>
C
      DIMENSION COOR(NCOOR), ATMTRA(NUCDEP,NUCDEP,NMTRX),
     &          CORTRA(3,3,NMTRX), TMPCOR(3)
C
      INMTRX = 0
C
C     *** Check how atoms are connected by the main rotation ***
C
      IF (MROTAX) THEN
         INMTRX = INMTRX + 1
         CALL DZERO(CORTRA(1,1,INMTRX),9)
C
         ROTM = D2*PI/DBLE(NCORDR)
C
         CORTRA(1,1,INMTRX) =  COS(ROTM)
         CORTRA(2,1,INMTRX) = -SIN(ROTM)
         CORTRA(1,2,INMTRX) =  SIN(ROTM)
         CORTRA(2,2,INMTRX) =  COS(ROTM)
         CORTRA(3,3,INMTRX) =  D1
C
         CALL FNDATR(ATMTRA(1,1,INMTRX),COOR,CORTRA(1,1,INMTRX),TMPCOR,
     &               'main rotational axis.')
      ELSE IF (ROTARE) THEN
         INMTRX = INMTRX + 1
         ROTM = D2*PI/DBLE(NCORDR)
C
         CORTRA(1,1,INMTRX) =  COS(ROTM)
         CORTRA(2,1,INMTRX) = -SIN(ROTM)
         CORTRA(1,2,INMTRX) =  SIN(ROTM)
         CORTRA(2,2,INMTRX) =  COS(ROTM)
         CORTRA(3,3,INMTRX) = -D1
C
         CALL FNDATR(ATMTRA(1,1,INMTRX),COOR,CORTRA(1,1,INMTRX),TMPCOR,
     &               'improper rotational axis.')
C
      END IF
C
      IF (ROTAX2) THEN
         INMTRX = INMTRX + 1
         CALL DZERO(CORTRA(1,1,INMTRX),9)
C
         CORTRA(1,1,INMTRX) =  D1
         CORTRA(2,2,INMTRX) = -D1
         CORTRA(3,3,INMTRX) = -D1
C
         CALL FNDATR(ATMTRA(1,1,INMTRX),COOR,CORTRA(1,1,INMTRX),TMPCOR,
     &               'rotational axis.')
      END IF
C
      IF (HPLANE) THEN
         INMTRX = INMTRX + 1
         CALL DZERO(CORTRA(1,1,INMTRX),9)
C
         CORTRA(1,1,INMTRX) =  D1
         CORTRA(2,2,INMTRX) =  D1
         CORTRA(3,3,INMTRX) = -D1
C
         CALL FNDATR(ATMTRA(1,1,INMTRX),COOR,CORTRA(1,1,INMTRX),TMPCOR,
     &               'horizontal plane.')
      END IF
C
      IF (VPLANE) THEN
         INMTRX = INMTRX + 1
         CALL DZERO(CORTRA(1,1,INMTRX),9)
C
         CORTRA(1,1,INMTRX) =  D1
         CORTRA(2,2,INMTRX) = -D1
         CORTRA(3,3,INMTRX) =  D1
C
         CALL FNDATR(ATMTRA(1,1,INMTRX),COOR,CORTRA(1,1,INMTRX),TMPCOR,
     &               'vertical mirror plane.')
      END IF
C
      IF (ICNTR) THEN
         INMTRX = INMTRX + 1
         CALL DZERO(CORTRA(1,1,INMTRX),9)
C
         CORTRA(1,1,INMTRX) = -D1
         CORTRA(2,2,INMTRX) = -D1
         CORTRA(3,3,INMTRX) = -D1
C
         CALL FNDATR(ATMTRA(1,1,INMTRX),COOR,CORTRA(1,1,INMTRX),TMPCOR,
     &               'iversion center.')
      END IF
C
      IF (IPRINT .GT. 20) THEN
C
         INMTRX = 0
C
         IF (MROTAX) THEN
            INMTRX = INMTRX + 1
            WRITE (LUPRI,'(A)') 'Atom transformation matrix' //
     &                          ' for main rotation.'
            DO J = 1, NUCDEP
               WRITE (LUPRI,'(20F5.1)') (ATMTRA(I,J,INMTRX),I=1, NUCDEP)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
            WRITE (LUPRI,'(A)') 'Coordinate transformation matrix' //
     &                          ' for main rotation.'
            DO J = 1, 3
               WRITE (LUPRI,'(20F8.4)') (CORTRA(I,J,INMTRX),I=1,3)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
         END IF
C
         IF (ROTARE) THEN
            INMTRX = INMTRX + 1
            WRITE (LUPRI,'(A)') 'Atom transformation matrix' //
     &                          ' for improper rotation.'
            DO J = 1, NUCDEP
               WRITE (LUPRI,'(20F5.1)') (ATMTRA(I,J,INMTRX),I=1, NUCDEP)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
            WRITE (LUPRI,'(A)') 'Coordinate transformation matrix' //
     &                          ' for improper rotation.'
            DO J = 1, 3
               WRITE (LUPRI,'(20F8.4)') (CORTRA(I,J,INMTRX),I=1,3)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
         END IF
C
         IF (ROTAX2) THEN
            INMTRX = INMTRX + 1
            WRITE (LUPRI,'(A)') 'Atom transformation matrix for' //
     &                          ' 2. rotation axis'
            DO J = 1, NUCDEP
               WRITE (LUPRI,'(20F5.1)') (ATMTRA(I,J,INMTRX),I=1, NUCDEP)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
C
            WRITE (LUPRI,'(A)') 'Coordinate transformation matrix' //
     &                          ' for 2. rotational axis.'
            DO J = 1, 3
               WRITE (LUPRI,'(20F8.4)') (CORTRA(I,J,INMTRX),I=1,3)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
         END IF
C
         IF (HPLANE) THEN
            INMTRX = INMTRX + 1
            WRITE (LUPRI,'(A)') 'Atom transformation matrix for' //
     &                          ' horizontal mirror plane'
            DO J = 1, NUCDEP
               WRITE (LUPRI,'(20F5.1)') (ATMTRA(I,J,INMTRX),I=1, NUCDEP)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
C
            WRITE (LUPRI,'(A)') 'Coordinate transformation matrix' //
     &                          ' for horizontal mirror plane.'
            DO J = 1, 3
               WRITE (LUPRI,'(20F8.4)') (CORTRA(I,J,INMTRX),I=1,3)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
         END IF
C
         IF (VPLANE) THEN
            INMTRX = INMTRX + 1
            WRITE (LUPRI,'(A)') 'Atom transformation matrix for' //
     &                          ' vertical mirror plane.'
            DO J = 1, NUCDEP
               WRITE (LUPRI,'(20F5.1)') (ATMTRA(I,J,INMTRX),I=1, NUCDEP)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
C
            WRITE (LUPRI,'(A)') 'Coordinate transformation matrix' //
     &                          ' for vertical mirror plane.'
            DO J = 1, 3
               WRITE (LUPRI,'(20F8.4)') (CORTRA(I,J,INMTRX),I=1,3)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
         END IF
C
         IF (ICNTR) THEN
            INMTRX = INMTRX + 1
            WRITE (LUPRI,'(A)') 'Atom transformation matrix for' //
     &                          ' inversion center'
            DO J = 1, NUCDEP
               WRITE (LUPRI,'(20F5.1)') (ATMTRA(I,J,INMTRX),I=1, NUCDEP)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
C
            WRITE (LUPRI,'(A)') 'Coordinate transformation matrix' //
     &                          ' for inversion center.'
            DO J = 1, 3
               WRITE (LUPRI,'(20F8.4)') (CORTRA(I,J,INMTRX),I=1,3)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
         END IF

      END IF
C
      RETURN
      END
C
C
C     /*Deck fndatr*/
      SUBROUTINE FNDATR(AMAT,ATMCOR,TRAMAT,TMPCOR,TEXT)
C
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
C
      PARAMETER (D1 = 1.0D0, DTHR = 1.0D-4)
C
#include <nuclei.h>
#include <trkoor.h>
#include <fcsym.h>
      CHARACTER*20 TEXT
      LOGICAL FOUND
      DIMENSION AMAT(NUCDEP,NUCDEP), ATMCOR(NCOOR), TRAMAT(3,3),
     &          TMPTRA(3,3), TMPCOR(3)
C
      CALL DZERO(AMAT,NUCDEP**2)
C
      DO 200 IATOM2 = 1, NUCDEP
         CALL DZERO(TMPCOR,3)
C
         ICOOR2 = 3*(IATOM2-1)
         DO 300 IX2 = 1, 3
         DO 300 IX1 = 1, 3
            TMPCOR(IX2) = TMPCOR(IX2)
     &                  + TRAMAT(IX1,IX2)*ATMCOR(ICOOR2+IX1)
 300     CONTINUE
C
         FOUND = .FALSE.
         DO 400 IATOM1 = 1, NUCDEP
            ICOOR1 = 3*(IATOM1-1)
            IF ((ABS(TMPCOR(1)-ATMCOR(ICOOR1+1)).LT.DTHR) .AND.
     &          (ABS(TMPCOR(2)-ATMCOR(ICOOR1+2)).LT.DTHR) .AND.
     &          (ABS(TMPCOR(3)-ATMCOR(ICOOR1+3)).LT.DTHR)) THEN
               FOUND = .TRUE.
               AMAT(IATOM2,IATOM1) = D1
            END IF
 400     CONTINUE
         IF (.NOT. FOUND) THEN

            WRITE (LUPRI,'(/A/A,I5)')
     &        'Error in transformation matrix for ' // TEXT,
     &        'Atom ',IATOM2
            CALL OUTPUT(TRAMAT,1,3,1,3,3,3,1,LUPRI)
            WRITE (LUPRI,'(/A,3F20.10//A)')
     &        'Transformed coordinates:',(TMPCOR(IX2),IX2=1,3),
     &        'All atom coordinates:'
            CALL OUTPUT(ATMCOR,1,3,1,NUCDEP,3,NUCDEP,1,LUPRI)
              
            CALL QUIT('You may have entered wrong group information.')
         END IF
 200  CONTINUE
C
      RETURN
      END
C
C
C     /*Deck mksymc*/
      SUBROUTINE MKSYMC(COOR,SYMCOR,ATMTRA,CORTRA,TRAMTX,TMPSYM,GRIREP,
     &                  TMPTRA,TMPVEC,WORK,ICRIRP,NSTBCT,LWORK,IPRINT)
C
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0, DTHR = 1.0D-6)
C
#include <nuclei.h>
#include <trkoor.h>
#include <fcsym.h>
      LOGICAL EXIST
      DIMENSION COOR(NCOOR), SYMCOR(NCOOR,NCOOR),
     &          ATMTRA(NUCDEP,NUCDEP,NMTRX), CORTRA(3,3,NMTRX),
     &          TRAMTX(NCOOR,NCOOR,NGORDR), TMPSYM(NCOOR),
     &          GRIREP(NGORDR,NGVERT), TMPTRA(NCOOR,NCOOR,2),
     &          TMPVEC(NCOOR), WORK(LWORK)
      DIMENSION NSTBCT(NUCDEP), ICRIRP(NCOOR,2)
C
      CALL DZERO(SYMCOR,NCOOR**2       )
      CALL DZERO(TRAMTX,NCOOR**2*NGORDR)
      CALL DZERO(ATMTRA,NUCDEP*NUCDEP*NMTRX)
      CALL DZERO(CORTRA,9*NMTRX)
C
      CALL TRMTRX(COOR,ATMTRA,CORTRA,IPRINT)
C
      CALL FNDSTB(COOR,ATMTRA,TRAMTX,TMPTRA,TMPVEC,NSTBCT,KINDCT,IPRINT)
C
      CALL DZERO(TRAMTX,NCOOR**2*NGORDR)
C
      DO 100 IORDR = 1, NCOOR
         TRAMTX(IORDR,IORDR,1) = D1
 100  CONTINUE
C
      IF (.NOT. MROTAX .AND. .NOT. ROTARE) THEN
C
         DO 200 INMTRX = 1, NMTRX
            DO 300 IATOM2 = 1, NUCDEP
            DO 300 IATOM1 = 1, NUCDEP
C
               DO 400 IX2 = 1, 3
               DO 400 IX1 = 1, 3
                  IORDR1 = 3*(IATOM1-1) + IX1
                  IORDR2 = 3*(IATOM2-1) + IX2
C
                  TRAMTX(IORDR1,IORDR2,INMTRX+1)
     &             = CORTRA(IX1,IX2,INMTRX)*ATMTRA(IATOM1,IATOM2,INMTRX)
 400           CONTINUE
 300        CONTINUE
 200     CONTINUE
C
      ELSE
         DO 500 IATOM2 = 1, NUCDEP
         DO 500 IATOM1 = 1, NUCDEP
C
            DO 600 IX2 = 1, 3
            DO 600 IX1 = 1, 3
               IORDR1 = 3*(IATOM1-1) + IX1
               IORDR2 = 3*(IATOM2-1) + IX2
C
               TRAMTX(IORDR1,IORDR2,2)
     &                       = CORTRA(IX1,IX2,1)*ATMTRA(IATOM1,IATOM2,1)
 600        CONTINUE
 500     CONTINUE
C
         DO 700 IROT = 2, NCORDR-1
C
            DO 701 ICOOR2 = 1, NCOOR
            DO 701 ICOOR1 = 1, NCOOR
               TMPTRA(ICOOR1,ICOOR2,2) = TRAMTX(ICOOR1,ICOOR2,2)
 701        CONTINUE
C
            DO 800 ITMP = 2, IROT
               IF (MOD(ITMP,2).EQ.0) THEN
                  IM1 = 1
                  IM2 = 2
               ELSE
                  IM1 = 2
                  IM2 = 1
               END IF
C
               CALL DZERO(TMPTRA(1,1,IM1),NCOOR**2)
C
               DO 900 ICOOR3 = 1, NCOOR
               DO 900 ICOOR2 = 1, NCOOR
               DO 900 ICOOR1 = 1, NCOOR
                  TMPTRA(ICOOR1,ICOOR3,IM1) = TMPTRA(ICOOR1,ICOOR3,IM1)
     &                                      + TMPTRA(ICOOR1,ICOOR2,IM2)
     &                                       *TRAMTX(ICOOR2,ICOOR3,2  )
 900           CONTINUE
 800        CONTINUE
C
            DO 1100 ICOOR2 = 1, NCOOR
            DO 1100 ICOOR1 = 1, NCOOR
               TRAMTX(ICOOR1,ICOOR2,IROT+1) = TMPTRA(ICOOR1,ICOOR2,IM1)
 1100       CONTINUE
C
 700     CONTINUE
C
         DO 1200 IELM = 1, NUMELM
            MLTMAX = IELM*NCORDR
C
            DO 1300 IATOM2 = 1, NUCDEP
            DO 1300 IATOM1 = 1, NUCDEP
C
               DO 1400 IX2 = 1, 3
               DO 1400 IX1 = 1, 3
                  ICOOR1 = 3*(IATOM1-1) + IX1
                  ICOOR2 = 3*(IATOM2-1) + IX2
C
                  TRAMTX(ICOOR1,ICOOR2,MLTMAX+1)
     &                                    = CORTRA(IX1   ,IX2   ,IELM+1)
     &                                     *ATMTRA(IATOM1,IATOM2,IELM+1)
 1400          CONTINUE
 1300       CONTINUE
C
            DO 1500 IOPR = 2, MLTMAX
               CALL DZERO(TRAMTX(1,1,MLTMAX+IOPR),NCOOR**2)
C
               DO 1600 ICOOR3 = 1, NCOOR
               DO 1600 ICOOR2 = 1, NCOOR
               DO 1600 ICOOR1 = 1, NCOOR
                  TRAMTX(ICOOR1,ICOOR3,MLTMAX+IOPR)
     &                               = TRAMTX(ICOOR1,ICOOR3,MLTMAX+IOPR)
     &                               + TRAMTX(ICOOR1,ICOOR2,MLTMAX+1   )
     &                                *TRAMTX(ICOOR2,ICOOR3,       IOPR)
 1600          CONTINUE
 1500       CONTINUE
 1200     CONTINUE
      END IF
C
      ISYMCO = 0
      IPLACE = 1
      DO 2100 IREP   = 1, NIREP
         NSDIM = 1
         IF (IREP.GT.N1DIME) NSDIM = 2
      DO 2101 IPL = 1, NSDIM
      DO 2101 ICENT  = 1, KINDCT
      DO 2101 ICOOR2 = 3*(NSTBCT(ICENT)-1)+1, 3*(NSTBCT(ICENT)-1)+3
         CALL DZERO(TMPSYM,NCOOR)
         IF (IREP .LE. N1DIME) THEN
            IPLACE1 = IREP
         ELSE
            IF (IPL.EQ.1) THEN
               IPLACE1 = (N1DIME + 1) + 4*(IREP - N1DIME - 1)
               IPLACE2 = IPLACE1 + 1
            ELSE
               IPLACE1 = N1DIME + 4*(IREP - N1DIME)
               IPLACE2 = IPLACE1 - 1
            END IF
         END IF
C
         DO 2200 ICOOR1 = 1, NCOOR
         DO 2200 IGORDR = 1, NGORDR
            TMPSYM(ICOOR1) = TMPSYM(ICOOR1)
     &                     + GRIREP(IGORDR,IPLACE1)
     &                      *TRAMTX(ICOOR1,ICOOR2,IGORDR)
 2200    CONTINUE
C
         RNORM2 = 0.0D0
         DO 2300 ICOOR1 = 1, NCOOR
            RNORM2 = RNORM2 + (TMPSYM(ICOOR1))**2
 2300    CONTINUE
         IF (RNORM2 .LT. DTHR) GOTO 2700
C
         DO 2350 ICOOR1 = 1, NCOOR
            TMPSYM(ICOOR1) = TMPSYM(ICOOR1)/SQRT(RNORM2)
 2350    CONTINUE
C
C        *** Orthogonalize it. ***
C
         EXIST = .FALSE.
         CALL ORTSCP(TMPSYM,SYMCOR,ICRIRP,ISYMCO,IREP,IPL,IPRINT,EXIST)
C
         DO 2400 IPHASE = 1, 2
         DO 2400 ICOO2  = 1, ISYMCO
            IF (.NOT.EXIST) THEN
               EXIST = .TRUE.
               PHASE = D1
               IF (IPHASE.EQ.2) PHASE = -D1
C
               DO 2500 ICOO1  = 1, NCOOR
                  DIFF2 = (PHASE*TMPSYM(ICOO1)
     &                  -        SYMCOR(ICOO1,ICOO2))**2
C
                  IF (DIFF2 .GT. DTHR) EXIST = .FALSE.
 2500          CONTINUE
C
            END IF
 2400    CONTINUE
         IF (EXIST) GOTO 2700
C
C        *** We have now found a new symmetry coordinate, and find ***
C        *** an adress for it.                                     ***
C
         ISYMCO = ISYMCO + 1
C
C        *** If the projection operator for the second row is used  ***
C        *** then we need to put the symmetry coordinate 1 index on ***
C
         ICRIRP(ISYMCO,2) = IPL-1
         ICRIRP(ISYMCO,1) = IREP
         DO ICOOR1 = 1, NCOOR
            SYMCOR(ICOOR1,ISYMCO) = TMPSYM(ICOOR1)
         END DO
C
C        *** Using the shift function to create a partner function. ***
C
         IF (IREP .GT. N1DIME) THEN
            CALL DZERO(TMPSYM,NCOOR)
C
            DO 2800 IGORDR = 1, NGORDR
            DO 2800 ICOO2  = 1, NCOOR
               CALL DZERO(TMPVEC,NCOOR)
               DO ICOO1 = 1, NCOOR
                  TMPVEC(ICOO2) = TMPVEC(ICOO2)
     &                          + TRAMTX(ICOO1,ICOO2,IGORDR)
     &                           *SYMCOR(ICOO1,ISYMCO)
               END DO
C
               TMPSYM(ICOO2) = TMPSYM(ICOO2)
     &                       + GRIREP(IGORDR,IPLACE2)
     &                        *TMPVEC(ICOO2)
 2800       CONTINUE
C
            RNORM2 = 0.0D0
            DO ICOOR1 = 1, NCOOR
               RNORM2 = RNORM2 + TMPSYM(ICOOR1)**2
            END DO
C
C           *** Normalizing partner geometry. ***
C
            DO ICOOR1 = 1, NCOOR
               TMPSYM(ICOOR1) = TMPSYM(ICOOR1)/SQRT(RNORM2)
            END DO
C
C           *** Orthogonalizing partner coordinate. ***
C
            CALL ORTSCP(TMPSYM,SYMCOR,ICRIRP,ISYMCO,IREP,IPL,IPRINT,
     &                  EXIST)
C
            ISYMCO = ISYMCO + 1
C
            ICRIRP(ISYMCO,1) = IREP
            ICRIRP(ISYMCO,2) = 1
            IF (IPL.EQ.2) ICRIRP(ISYMCO,2) = 0
C
            DO ICOOR1 = 1, NCOOR
               SYMCOR(ICOOR1,ISYMCO) = TMPSYM(ICOOR1)
            END DO
C
         END IF
C
 2700    CONTINUE
 2101 CONTINUE
 2100 CONTINUE
C
C
      IF (N1DIME.LT.NIREP) THEN
         KTMPCR = 1
         KICTMP = KTMPCR + NCOOR**2
         CALL SRTSCR(SYMCOR,WORK(KTMPCR),ICRIRP,WORK(KICTMP),IPRINT)
      END IF
C
      NUMTIM = (NCOOR-1)/10 + 1
      CALL HEADER ('Symmetry coordinates',-1)
      DO I = 1, NUMTIM
         ILEFT = NCOOR - 10*(I-1)
         ISTRT = 10*(I-1) + 1
         IEND  =(ISTRT-1) + MIN(ILEFT,10)
         WRITE (LUPRI,'(A,I4,10I8)') 'Symmetry',
     &        (ICRIRP(ICOOR2,1),ICOOR2=ISTRT,IEND)
         DO ICOOR1 = 1, NCOOR
            WRITE (LUPRI,'(A,10F8.4)') '      ',
     &             (SYMCOR(ICOOR1,ICOOR2),ICOOR2=ISTRT,IEND)
         END DO
         WRITE (LUPRI,'(A)') '                             '
         WRITE (LUPRI,'(A)') '                             '
      END DO
      WRITE (LUPRI,'(A)') '                                '
C
      IF (IPRINT .GT. 20) THEN
         WRITE (LUPRI,'(A)') '                                '
         WRITE (LUPRI,'(A)') 'Transformation matrices for the ' //
     &                       'different symmetry operations.'
         WRITE (LUPRI,'(A)') '                                '
         DO IGORDR = 1, NGORDR
            WRITE (LUPRI,'(A,I4)') ' Matrix number: ', IGORDR
            DO I = 1, NUMTIM
               ILEFT = NCOOR - 18*(I-1)
               ISTRT = 18*(I-1) + 1
               IEND  = (ISTRT-1) + MIN(ILEFT,18)
               DO IORDR2 = 1, NCOOR
                  WRITE (LUPRI,'(18F8.4)')
     &               (TRAMTX(IORDR1,IORDR2,IGORDR),IORDR1=ISTRT,IEND)
               END DO
               WRITE (LUPRI,'(A)') '                              '
               WRITE (LUPRI,'(A)') '                              '
            END DO
            WRITE (LUPRI,'(A)') '                              '
         END DO
      END IF
C
      IF (ISYMCO.LT.NCOOR) CALL QUIT('Unable to make enough sym. coor')
C
C     *** Orthogonality test for coordinates. ***
C
      ITEST = 0
      DO ICOOR1 = 1, NCOOR
      DO ICOOR2 = 1, ICOOR1-1
         GCCNST = 0.0D0
         DO ICOOR3 = 1, NCOOR
            GCCNST = GCCNST + SYMCOR(ICOOR3,ICOOR1)
     &                       *SYMCOR(ICOOR3,ICOOR2)
         END DO
C
         IF (ABS(GCCNST) .GT. DTHR) THEN
            ITEST = ITEST + 1
            WRITE (LUPRI,'(5X,A,2I4)')
     &     'ERROR, Non orthogonal symmetry coordinates:', ICOOR1, ICOOR2
         END IF
      END DO
      END DO
      IF (ITEST.GT.0) CALL QUIT('Problems with symmetry coordinates,' //
     &     'numerical differentiation. Please contact dalton admin.')
C
      RETURN
      END
C
C
C     /*Deck fndstb*/
      SUBROUTINE FNDSTB(COOR,ATMTRA,TRAMTX,TMPTRA,TMPVEC,NSTBCT,KINDCT,
     &                  IPRINT)
C
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0, DTHR = 1.0D-6)
#include <nuclei.h>
#include <trkoor.h>
#include <fcsym.h>
C
      LOGICAL EXIST
      DIMENSION ATMTRA(NUCDEP,NUCDEP,NMTRX), TMPTRA(NCOOR,NCOOR,2),
     &          TRAMTX(NCOOR,NCOOR,NGORDR), NSTBCT(NUCDEP),
     &          TMPVEC(NCOOR), COOR(NCOOR)
C
      DO 100 IATOM = 1, NUCDEP
         TRAMTX(IATOM,IATOM,1) = D1
 100  CONTINUE
C
      IF (.NOT. MROTAX .AND. .NOT. ROTARE) THEN
C
         DO 200 INMTRX = 1, NMTRX
         DO 200 IATOM2 = 1, NUCDEP
         DO 200 IATOM1 = 1, NUCDEP
            TRAMTX(IATOM1,IATOM2,INMTRX+1)= ATMTRA(IATOM1,IATOM2,INMTRX)
 200     CONTINUE
C
      ELSE
         DO 300 IATOM2 = 1, NUCDEP
         DO 300 IATOM1 = 1, NUCDEP
            TRAMTX(IATOM1,IATOM2,2) = ATMTRA(IATOM1,IATOM2,1)
 300     CONTINUE
C
         DO 400 IROT = 2, NCORDR-1
C
            DO 500 IATOM2 = 1, NCOOR
            DO 500 IATOM1 = 1, NCOOR
               TMPTRA(IATOM1,IATOM2,2) = TRAMTX(IATOM1,IATOM2,2)
 500        CONTINUE
C
            DO 600 ITMP = 2, IROT
               IF (MOD(ITMP,2).EQ.0) THEN
                  IM1 = 1
                  IM2 = 2
               ELSE
                  IM1 = 2
                  IM2 = 1
               END IF
C
               CALL DZERO(TMPTRA(1,1,IM1),NCOOR**2)
C
               DO 700 IATOM3 = 1, NUCDEP
               DO 700 IATOM2 = 1, NUCDEP
               DO 700 IATOM1 = 1, NUCDEP
                  TMPTRA(IATOM1,IATOM3,IM1) = TMPTRA(IATOM1,IATOM3,IM1)
     &                                      + TMPTRA(IATOM1,IATOM2,IM2)
     &                                       *TRAMTX(IATOM2,IATOM3,2  )
 700           CONTINUE
 600        CONTINUE
C
            DO 800 IATOM2 = 1, NCOOR
            DO 800 IATOM1 = 1, NCOOR
               TRAMTX(IATOM1,IATOM2,IROT+1) = TMPTRA(IATOM1,IATOM2,IM1)
 800        CONTINUE
C
 400     CONTINUE
C
         DO 900 IELM = 1, NUMELM
            MLTMAX = IELM*NCORDR
C
            DO 1100 IATOM2 = 1, NUCDEP
            DO 1100 IATOM1 = 1, NUCDEP
C
               TRAMTX(IATOM1,IATOM2,MLTMAX+1)
     &                                    = ATMTRA(IATOM1,IATOM2,IELM+1)
 1100       CONTINUE
C
            DO 1200 IOPR = 2, MLTMAX
               CALL DZERO(TRAMTX(1,1,MLTMAX+IOPR),NCOOR**2)
C
               DO 1300 IATOM3 = 1, NUCDEP
               DO 1300 IATOM2 = 1, NUCDEP
               DO 1300 IATOM1 = 1, NUCDEP
                  TRAMTX(IATOM1,IATOM3,MLTMAX+IOPR)
     &                               = TRAMTX(IATOM1,IATOM3,MLTMAX+IOPR)
     &                               + TRAMTX(IATOM1,IATOM2,MLTMAX+1   )
     &                                *TRAMTX(IATOM2,IATOM3,       IOPR)
 1300          CONTINUE
 1200       CONTINUE
 900     CONTINUE
      END IF
C
      KINDCT = 0
      CALL IZERO(NSTBCT,NUCDEP)
      DO 1400 IATOM = 1, NUCDEP
C
         EXIST = .FALSE.
         DO 1500 IORDR  = 1, NGORDR
         DO 1500 IINDCT = 1, KINDCT
            IF (ABS(TRAMTX(NSTBCT(IINDCT),IATOM,IORDR)).GT.DTHR) THEN
               GOTO 1600
            END IF
 1500    CONTINUE
C
         KINDCT = KINDCT + 1
         NSTBCT(KINDCT) = IATOM
C
 1600    CONTINUE
 1400 CONTINUE
C
c      DO 1700 IINDCT = 1, KINDCT
c         IF ((FCLASS(1:1).EQ.'S').OR.(FCLASS(3:3).EQ.'d')) THEN
c            IYCOOR = 3*(NSTBCT(INDCT)-1) + 1
c         ELSE
c            IYCOOR = 3*(NSTBCT(IINDCT)-1) + 2
c         END IF
cC
c         IF (ABS(COOR(IYCOOR)) .LT. DTHR) THEN
cC
c            DO 1800 IORDR = 2, NGORDR
cC
c               DO 1900 IATOM = 1, NUCDEP
c                  TMPVEC(IATOM) = TRAMTX(NSTBCT(IINDCT),IATOM,IORDR)
c                  IF (ABS(TMPVEC(IATOM)).GT.DTHR)
c     &                                          NSTBCT(IINDCT) = IATOM
c 1900          CONTINUE
cC
c               IYCOOR = 3*(NSTBCT(IINDCT)-1) + 2
c               IF (ABS(COOR(IYCOOR)).GT.DTHR) GOTO 2100
c 1800       CONTINUE
c         END IF
cC
c 2100    CONTINUE
cC
c 1700 CONTINUE
C
      IF (IPRINT .GT. 20) THEN
         WRITE (LUPRI,'(A)') 'Atom transformation matrices:'
C
         DO IMTX = 1, NGORDR
            WRITE (LUPRI,'(A,I2)') 'Operation ', IMTX
            DO I = 1, NUCDEP
               WRITE (LUPRI,'(12F5.1)') (TRAMTX(I,J,IMTX),J=1, NUCDEP)
            END DO
            WRITE (LUPRI,'(A)') '                                '
         END DO
C
         WRITE (LUPRI,'(A)') 'Symmetry independent centers:'
         DO IINDCT = 1, KINDCT
            WRITE (LUPRI,'(A,I4)') '    ', NSTBCT(IINDCT)
         END DO
      END IF
C
      RETURN
      END
C
C
C     /*Deck fcscrn*/
      SUBROUTINE FCSCRN(GRIREP,WORK,KDPMTX,INDSTP,ICRIRP,IRPIND,IDXTMP,
     &                  IDDBTP,IRPDEG,IREPST,NPRTNR,LWORK,NLDPMX,LDPMTX,
     &                  IFRSTD,IORDR,IRSRDR,MAXINR,IINNER,NMPRTN,IPRINT,
     &                  CLNRGY,PRTNR,ALRCAL,SCND)
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
C
#include <trkoor.h>
#include <numder.h>
#include <fcsym.h>
#include <pvibav.h>
      LOGICAL CLNRGY, PRTNR, SCND, C2NRGY, DEPFC, FOUND, ALRCAL
      DIMENSION GRIREP(NGORDR,NGVERT), WORK(LWORK),
     &          INDSTP(NTORDR), KDPMTX(LDPMTX,NSTRDR,IFRSTD),
     &          IRPDEG(NMORDR), IREPST(NMORDR), ICRIRP(NCOOR,2),
     &          IRPIND(NMORDR), IDXTMP(NMORDR), IDDBTP(NMORDR),
     &          NPRTNR(MAXINR)
C
      IF ((FCLASS .NE. 'C1 ').AND.(NAORDR.LT.1).AND..NOT.CNMPRP) THEN
C
C        ***Check symmetry for the original component ***
C
         KKIRPD = 1
         KIDDEG = KKIRPD + NMORDR
         KIDEGI = KIDDEG + NMORDR
C
         C2NRGY = .FALSE.
         CALL SYMMLT(GRIREP,INDSTP,ICRIRP,IREPST,IRPDEG,WORK(KKIRPD),
     &               WORK(KIDDEG),WORK(KIDEGI),NMORDR,NTORDR,IRSRDR,
     &               IORDR,IPRINT,C2NRGY)
C
C        *** Checkicng whether this is a dependent force constant ***
C        *** (as a result of using groups with degenerate ireps)  ***
C
         DEPFC = .FALSE.
         IF (C2NRGY) THEN
            CALL FFCDEP(KDPMTX,INDSTP,IDXTMP,NLDPMX,LDPMTX,IFRSTD,
     &                  IRSRDR,IORDR,0,IPRINT,DEPFC)
         END IF
C
C        *** Final test for this component ***
C
         CLNRGY = C2NRGY.AND..NOT.DEPFC
C
C        *** Check if it is possible to construct a force-constant ***
C        *** with two additional equal and arbitary indices. This  ***
C        *** force constant will contain this energy value.        ***
C
         IF (((NMORDR-IORDR).GE.2).AND.(.NOT.CLNRGY)
     &                            .AND.(.NOT.SCND  )) THEN
            ISTRT = IRSRDR + 1
            INUM = (NMORDR-IORDR)/2
C
            NSTP = NDCOOR
            IDDBTP(1) = 0
            DO 100 II = 2, INUM
               NSTP = NSTP*(NDCOOR+II-1)
               IDDBTP(II) = 1
 100        CONTINUE
C
            DO 200 ISTP = 1, NSTP
               IF (.NOT.CLNRGY) THEN
                  FOUND = .FALSE.
                  DO 300 II = INUM-1, 1, -1
                     IF ((IDDBTP(II).GT.IDDBTP(II+1))
     &                                       .AND.(.NOT.FOUND)) THEN
                        FOUND = .TRUE.
                        IDDBTP(II+1) = IDDBTP(II+1) + 1
                     END IF
 300              CONTINUE
                  IF (.NOT.FOUND) IDDBTP(1) = IDDBTP(1) + 1
C
                  IDX = 0
                  DO 400 IJ = 1, INUM
                  DO 400 II = 1, 2
                     IDX = IDX + 1
                     INDSTP(ISTRT+IDX) = IDDBTP(IJ)
 400              CONTINUE
C
                  C2NRGY = .FALSE.
                  CALL SYMMLT(GRIREP,INDSTP,ICRIRP,IREPST,IRPDEG,
     &                        WORK(KKIRPD),WORK(KIDDEG),WORK(KIDEGI),
     &                        NMORDR,NTORDR,IRSRDR+2*INUM,IORDR+2*INUM,
     &                        IPRINT,C2NRGY)
C
C                  *** Checking if this is a dependent component. ***
C
                  DEPFC = .FALSE.
                  IF (C2NRGY) THEN
                     CALL FFCDEP(KDPMTX,INDSTP,IDXTMP,NLDPMX,LDPMTX,
     &                           IFRSTD,IRSRDR,IORDR,2*INUM,IPRINT,
     &                           DEPFC)
                  END IF
C
C                 *** Final test for this component. ***
C
                  CLNRGY = C2NRGY.AND..NOT.DEPFC
               END IF
 200        CONTINUE
C
         END IF
C
C        *** Check symmetry for use of component in higher derivatives ***
C        *** MIN(NMORDR-IORDR,IORDR) -> if calculating derivatives of  ***
C        *** order > 2*iordr then clnrgy = .true.                      ***
C
         IF ((IORDR.LT.NMORDR) .AND. .NOT.CLNRGY .AND. .NOT.SCND) THEN
            DO 500 IN = 1, NMORDR-IORDR
C
               NRP = 1
               DO 600 IRP = 1, IN
                  NRP = NRP*(IRSRDR + IRP)/IRP
 600           CONTINUE
C
C              *** Which components do we have to add ***
C
               CALL IZERO(IRPIND,NMORDR)
               IF (IORDR.EQ.1) THEN
                  IRPIND(1) = 1
               ELSE
                  IRPIND(1) = 0
               END IF
               DO 700 I = 2, IN
                  IRPIND(I) = 1
 700           CONTINUE
C
               DO 800 IRP = 1, NRP
C
                  IF (IORDR.GT.1) THEN
                     IF (IN.GT.1) THEN
                        FOUND = .FALSE.
                        DO 900 IIN = 1, IN
                           IF ((IRPIND(IIN).LE.IRPIND(IIN+1)).AND.
     &                                           (.NOT.FOUND)) THEN
                              FOUND = .TRUE.
                              IRPIND(IIN) = IRPIND(IIN) + 1
                           END IF
 900                    CONTINUE
                     ELSE
                        IRPIND(1) = IRPIND(1) + 1
                     END IF
                  END IF
C
C                 *** Assign the new components already in the array, ***
C                 ***               and screen again.                 ***
C
                  IF (.NOT. CLNRGY) THEN
                     ISTRT = IRSRDR + 1
                     DO 1100 IIND = 1, IN
                        INDSTP(ISTRT+IIND) = INDSTP(IRPIND(IIND))
 1100                CONTINUE
C
                     C2NRGY = .FALSE.
                     CALL SYMMLT(GRIREP,INDSTP,ICRIRP,IREPST,IRPDEG,
     &                           WORK(KKIRPD),WORK(KIDDEG),WORK(KIDEGI),
     &                           NMORDR,NTORDR,IRSRDR+IN,IORDR+IN,
     &                           IPRINT,C2NRGY)
C
C                    *** Checking if this is a dependent component. ***
C
                     DEPFC = .FALSE.
                     IF (C2NRGY) THEN
                        CALL FFCDEP(KDPMTX,INDSTP,IDXTMP,NLDPMX,LDPMTX,
     &                              IFRSTD,IRSRDR,IORDR,IN,IPRINT,DEPFC)
                     END IF
C
C                    *** Final test for this component. ***
C
                     CLNRGY = C2NRGY.AND..NOT.DEPFC
                  END IF
C
C                 *** Check if it is possible to construct a force-constant ***
C                 *** with two additional equal and arbitary indices. This  ***
C                 *** force constant will contain this energy value.        ***
C
                  IF (((NMORDR-(IORDR+IN)).GE.2).AND.(.NOT.CLNRGY)) THEN
C
                    ISTRT = IRSRDR + 1 + IN
                    INUM = (NMORDR-(IORDR+IN))/2
C
                    NSTP = NDCOOR
                    IDDBTP(1) = 0
                    DO 1200 II = 2, INUM
                       NSTP = NSTP*(NDCOOR+II-1)
                       IDDBTP(II) = 1
 1200               CONTINUE
C
                    DO 1300 ISTP = 1, NSTP
                       IF (.NOT.CLNRGY) THEN
                          FOUND = .FALSE.
                          DO 1400 II = INUM-1, 1, -1
                             IF ((IDDBTP(II).GT.IDDBTP(II+1))
     &                                      .AND.(.NOT.FOUND)) THEN
                                FOUND = .TRUE.
                                IDDBTP(II+1) = IDDBTP(II+1) + 1
                             END IF
 1400                     CONTINUE
                          IF (.NOT.FOUND) IDDBTP(1) = IDDBTP(1) + 1
C
                          IDX = 0
                          DO 1500 IJ = 1, INUM
                          DO 1500 II = 1, 2
                             IDX = IDX + 1
                             INDSTP(ISTRT+IDX) = IDDBTP(IJ)
 1500                     CONTINUE
C
                          C2NRGY = .FALSE.
                          CALL SYMMLT(GRIREP,INDSTP,ICRIRP,IREPST,
     &                                IRPDEG,WORK(KKIRPD),WORK(KIDDEG),
     &                                WORK(KIDEGI),NMORDR,NTORDR,
     &                                IRSRDR+IN+2*INUM,IORDR+IN+2*INUM,
     &                                IPRINT,C2NRGY)
C
C                         *** Checking if this is a dependent component. ***
C
                          DEPFC = .FALSE.
                          IF (C2NRGY) THEN
                             CALL FFCDEP(KDPMTX,INDSTP,IDXTMP,NLDPMX,
     &                                   LDPMTX,IFRSTD,IRSRDR,IORDR,
     &                                   IN+2*INUM,IPRINT,DEPFC)
                          END IF
C
C                         *** Final test for this component. ***
C
                          CLNRGY = C2NRGY.AND..NOT.DEPFC
                       END IF
 1300               CONTINUE
C
                 END IF
C
 800           CONTINUE
 500        CONTINUE
         END IF
C
C        *** Check if further simplyfication is possible through ***
C        *** the use of partner geometry's                       ***
C
         IF (CLNRGY.AND.(NAORDR.EQ.0).AND.(.NOT.SCND)) THEN
            KIDCTP = 1
            KINDTP = KIDCTP + NMORDR
            KLAST  = KINDTP + 5
            LWRK   = LWORK  - KLAST + 1
            CALL PRTNRP(GRIREP,WORK(KLAST),ICRIRP,INDSTP,NPRTNR,MAXINR,
     &                  IINNER,IORDR,IRSRDR,NMPRTN,WORK(KIDCTP),
     &                  WORK(KINDTP),LWRK,CLNRGY,PRTNR,ALRCAL)
         END IF
C
      ELSE
         CLNRGY = .TRUE.
      END IF
C
c      if (iordr.eq.2) then
c         if ((icrirp(indstp(1),1).eq.5).and.
c     &       (icrirp(indstp(1),2).eq.1).and.
c     &       (icrirp(indstp(2),1).eq.5).and.
c     &       (icrirp(indstp(2),2).eq.0))
c     &        stop ' '
c      end if
      RETURN
      END
C
C
C     /*Deck symmlt*/
      SUBROUTINE SYMMLT(GRIREP,INDSTP,ICRIRP,IREPST,IRPDEG,KIRPDG,
     &                  IDDEGI,IDEGID,NMORDR,NTORDR,IRSRDR,IORDR,IPRINT,
     &                  CLNRGY)
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, DTHR = 1.0D-12)
#include <trkoor.h>
#include <fcsym.h>
      LOGICAL CLNRGY
      DIMENSION GRIREP(NGORDR,NGVERT)
      DIMENSION INDSTP(NTORDR), KIRPDG(NMORDR), IRPDEG(NMORDR),
     &          IREPST(NMORDR), IDEGID(NMORDR), IDDEGI(NMORDR),
     &          ICRIRP(NCOOR,2)
C
C     *** Starting points for the relevant rows in the relevant  ***
C     *** irreps For instance, the starting point for the 1. row ***
C     *** of E in C3v is 2, while the second row is 4.           ***
C     *** IREPST -> Starting point of irrep                      ***
C     *** KIRPDG -> dimension of the irrep                       ***
C     *** IDEGGI -> coordinates that belong to 2 dim. irreps     ***
C
      NMAX = 0
      IC   = 0
      DO 100 IC2 = 1, IRSRDR+1
         NC1 = 1
         IF (IC2.EQ.1) NC1 = IORDR-IRSRDR
         DO 200 IC1 = 1, NC1
            IC = IC + 1
            IREPST(IC) = IRPSTR(GRIREP,ICRIRP,NCOOR,INDSTP(IC2),IPRINT)
            IF (ICRIRP(INDSTP(IC2),1) .GT. N1DIME) THEN
               NMAX = NMAX + 1
               KIRPDG(IC) = 2
               IDDEGI(NMAX) = IC
            ELSE
               KIRPDG(IC) = 1
            END IF
 200     CONTINUE
 100  CONTINUE
C
C     *** IRPDEG -> The row number we want use in our product -  ***
C     *** IREPST (irep-start). For two 2-dimensional irreps this ***
C     *** is 11, 12, 21 and 22                                   ***
C
      DO 300 IC = 1, IORDR
         IRPDEG(IC) = 1
 300  CONTINUE
C
C     *** Number of 2 dimensional irreps. ***
C
      DO 400 IMAX = 0, NMAX
         INUM = 1
         IDEN = 1
         DO 500 I = 1, IMAX
            INUM = INUM*(NMAX-I+1)
            IDEN = IDEN*I
 500     CONTINUE
         NPOSBL = INUM/IDEN
C
C        *** Number of possible permutations for IMAX values of 2 ***
C
         DO 600 IPOSBL = 1, NPOSBL
C
            DO 700 IC = 1, IORDR
               IRPDEG(IC) = 1
 700        CONTINUE
C
C           *** IMAX .gt. 0 -> at least 1 two dimensional irrep. ***
C           *** IDEGID -> which values of IDDEGID is 2           ***
C
            IF (IMAX .GT. 0) THEN
               IF (IPOSBL.EQ.1) THEN
                  DO 800 IRDR = 1, IMAX
                     IDEGID(IRDR) = IMAX - IRDR + 1
 800              CONTINUE
               ELSE
                  DO 900 IRDR = 1, IMAX
                     IF (IDEGID(IRDR).LE.NMAX-IRDR) THEN
                        IDEGID(IRDR) = IDEGID(IRDR) + 1
                        GOTO 1100
                     END IF
 900              CONTINUE
               END IF
            END IF
C
 1100       CONTINUE
C
C           *** Assigning the nesscecary values of 2 ***
C
            DO 1200 IRDR = 1, IMAX
               IRPDEG(IDDEGI(IDEGID(IRDR))) = 2
 1200       CONTINUE
C
C           *** Calculating the screening condition -> P. Taylor, ***
C           *** Lecture notes in quantum chemistry 56, Springer   ***
C
            SMTXEL = D0
            DO 2100 IORDER = 1, NGORDR
C
C              *** Finding the starting point for the  ***
C              *** screening needed for the analytical ***
C              *** derivative. Output: RTMPCO          ***
C
               CALL ANLSYM(GRIREP,ICRIRP,INDSTP,IREPST,RTMPCO,NCOOR,
     &                     IORDER,IORDR,IPRINT)
C
               DO 2200 IC = 1, IORDR
C
C                 *** Special care for groups based on separable  ***
C                 *** degenerate groups, Sn, Cnh Dnd. GOT doesn't ***
C                 *** apply, only LOT.                            ***
C
                  IF (SEPDEG) THEN
                     IFAC = 2
                     IF (ICRIRP(INDSTP(IC),2).EQ.1) IFAC = -2
                     FMULT = GRIREP(IORDER,IREPST(IC)+IRPDEG(IC)     )
     &                     + GRIREP(IORDER,IREPST(IC)+IRPDEG(IC)+IFAC)
                  ELSE
                     FMULT = GRIREP(IORDER,IREPST(IC)+IRPDEG(IC))
                  END IF
                  RTMPCO = RTMPCO*FMULT
 2200          CONTINUE
C
               SMTXEL = SMTXEL + RTMPCO
 2100       CONTINUE
C
            IF (ABS(SMTXEL) .GT. DTHR) CLNRGY = .TRUE.
 600     CONTINUE
 400  CONTINUE
C
      RETURN
      END
C
C
C     /* Deck irpstr */
      FUNCTION IRPSTR(GRIREP,ICRIRP,NCOOR,ICOOR,IPRINT)
C     ******************************************************
C     *** Function that takes a coordinate as input, and ***
C     *** returns the starting point of the irrep that   ***
C     *** span that irrep in GRIREP.                     ***
C     ******************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <fcsym.h>
      DIMENSION GRIREP(NGORDR,NGVERT)
      DIMENSION ICRIRP(NCOOR,2)
C
C     *** Finding the starting point. ***
C
      IF (ICRIRP(ICOOR,1) .GT. N1DIME) THEN
         ISTRTT = N1DIME + 4*(ICRIRP(ICOOR,1)-N1DIME-1)
         IF (ICRIRP(ICOOR,2).EQ.1) THEN
            ISTRTT = ISTRTT + 2
         END IF
      ELSE
         ISTRTT = ICRIRP(ICOOR,1)-1
      END IF
C
C     *** Assigning the value. ***
C
      IRPSTR = ISTRTT
C
      RETURN
      END
C
C
C     /* Deck anasym */
      SUBROUTINE ANLSYM(GRIREP,ICRIRP,INDSTP,IREPST,RTMPCO,NCOOR,IORDER,
     &                  IORDR,IPRINT)
C     *******************************************************
C     *** Subroutine that decides how the symmetry is for ***
C     *** analytical derivatives for this numerical       ***
C     *** distortion. Output parameter is RTMPCO, which is***
C     *** assigned a value according to the symmetry prop.***
C     *** of the analytical derivative.                   ***
C     *******************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0, D0 = 0.0D0)
#include <numder.h>
#include <fcsym.h>
      LOGICAL SBDELM, SMEIRP
      DIMENSION GRIREP(NGORDR,NGVERT)
      DIMENSION INDSTP(NTORDR), IREPST(NMORDR), ICRIRP(NCOOR,2)
C
      IF (NAORDR.EQ.0) THEN
         RTMPCO=D1
      ELSE IF (NAORDR.EQ.1) THEN
         RTMPCO = D0
         DO ICOOR = 1, NCOOR
            IPSTRT = IRPSTR(GRIREP,ICRIRP,NCOOR,ICOOR,IPRINT)
            SBDELM = SMEIRP(GRIREP,ICRIRP,INDSTP,IREPST,ICOOR,IPSTRT,
     &                      IORDR,NCOOR,IPRINT)
            IF (SBDELM) THEN
               NSTEP = 1
               IF (ICRIRP(ICOOR,1).GT.N1DIME) NSTEP = 2
               DO ISTP = 1, NSTEP
                  RTMPCO = RTMPCO + GRIREP(IORDER,IPSTRT+ISTP)
               END DO
            END IF
         END DO
      END IF
C
      RETURN
      END
C
C
C     /* Deck smeirp */
      LOGICAL FUNCTION SMEIRP(GRIREP,ICRIRP,INDSTP,IREPST,ICOOR,IPSTRT,
     &                        IORDR,NCOOR,IPRINT)
C     *********************************************************
C     *** Subroutine that checks if the coordinate ICOOR is ***
C     *** totally symmetric in the distorted geometry. The  ***
C     *** distortions is done along the INDSTP coordinates. ***
C     *********************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1=1.0D0, DTHRS=1.0D-12)
#include <fcsym.h>
#include <numder.h>
      LOGICAL EXSTEL, TSMEIR
      DIMENSION GRIREP(NGORDR,NGVERT)
      DIMENSION INDSTP(NTORDR), ICRIRP(NCOOR,2), IREPST(NMORDR)
C
      TSMEIR = .TRUE.
      DO IEL = 1, NGORDR
C
C        *** Checkin if this element is in the subgroup. ***
C
         EXSTEL = .TRUE.
         DO IC = 1, IORDR
            KSTP = ICRIRP(INDSTP(IC),2) + 1
            IF ((GRIREP(IEL,IREPST(IC)+KSTP)-D1)**2 .GT. DTHRS)
     &                                                  EXSTEL = .FALSE.
         END DO
C
C        *** If this element exist, check if irep is totally ***
C        *** symmetric.                                      ***
C
         KSTP = ICRIRP(ICOOR,2) + 1
         IF (EXSTEL.AND.((GRIREP(IEL,IPSTRT+KSTP)-D1)**2.GT.DTHRS))
     &                                                  TSMEIR = .FALSE.
C
      END DO
C
C     *** Assigning the result to the function. ***
C
      SMEIRP = TSMEIR
C
      RETURN
      END
C
C
C     /* Deck ffcdep */
      SUBROUTINE FFCDEP(KDPMTX,INDSTP,IDXTMP,NLDPMX,LDPMTX,IFRSTD,
     &                  IRSRDR,IORDR,IEIIN,IPRINT,DEPFC)
#include <implicit.h>
Chjaaj DEBUG priunit.h
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DTHR = 1.0D-12)
#include <fcsym.h>
#include <numder.h>
      LOGICAL DONE, DEPFC
      DIMENSION KDPMTX(LDPMTX,NSTRDR,IFRSTD), INDSTP(NTORDR),
     &          IDXTMP(NMORDR)
C
      IDXTMP(1) = INDSTP(1)
C
      DO 100 II = 2, IORDR-IRSRDR
         IDXTMP(II) = INDSTP(1)
 100  CONTINUE
C
      ISTRT = IORDR-IRSRDR
      DO 200 II = 1, IRSRDR
         IDXTMP(ISTRT+II) = INDSTP(II+1)
 200  CONTINUE
C
      DO 300 II = ISTRT+IRSRDR+1, ISTRT+IRSRDR+IEIIN
         DONE = .FALSE.
         DO 400 IJ = ISTRT, II-1
            IF ((IDXTMP(IJ).LT.INDSTP(II)).AND.(.NOT.DONE)) THEN
               DO 500 IK = II, IJ+1, -1
                  IDXTMP(IK) = IDXTMP(IK-1)
 500           CONTINUE
               IDXTMP(IJ) = INDSTP(II)
               DONE = .TRUE.
            END IF
 400     CONTINUE
         IF (.NOT.DONE) IDXTMP(II) = INDSTP(II)
 300  CONTINUE
C
      ITHDIM = ISTRT+IRSRDR+IEIIN
      if (ITHDIM+1 .gt. NSTRDR) then
         write(lupri,*) 'FFCDEP error: ITHDIM+1 .gt. NSTRDR'
         write(lupri,*)
     &   'FFCDEP error: ITHDIM,ISTRT,IRSRDR,IEIIN,NSTRDR',
     &    ITHDIM,ISTRT,IRSRDR,IEIIN,NSTRDR
         call quit('error in FFCDEP')
      end if
      MINTST = MIN(ITHDIM,NMORDR)
      DO 600 II = 1, NLDPMX
         IF (((KDPMTX(II,ITHDIM+1,1).EQ.0).OR.(MINTST.EQ.NMORDR)).AND.
     &                                                (.NOT.DEPFC)) THEN
            DEPFC = .TRUE.
            DO 700 IRDR = 1, ITHDIM
               DEPFC = DEPFC.AND.(KDPMTX(II,IRDR,1).EQ.IDXTMP(IRDR))
 700        CONTINUE
         END IF
 600  CONTINUE
C
      RETURN
      END
C
C
C     /*Deck ortscp*/
      SUBROUTINE ORTSCP(TRIAL,SYMCOR,ICRIRP,ISYMCO,KIREP,KIPL,IPRINT,
     &                  EXIST)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DTHR = 1.0D-10)
#include <trkoor.h>
      LOGICAL EXIST
      DIMENSION TRIAL(NCOOR), SYMCOR(NCOOR,NCOOR), ICRIRP(NCOOR,2)
C
      DO  ICOOR = 1, ISYMCO
         IREP = ICRIRP(ICOOR,1)
C
         IF ((IREP .EQ. KIREP) .AND. (ICRIRP(ICOOR,2).EQ.(KIPL-1))) THEN
C
            RLENGS = D0
            DO ICOOR1 = 1, NCOOR
               RLENGS = RLENGS + SYMCOR(ICOOR1,ICOOR)**2
            END DO
            RINVS = D1/RLENGS
C
            GCCNST = D0
            DO ICOOR1 = 1, NCOOR
               GCCNST = GCCNST + SYMCOR(ICOOR1,ICOOR)*TRIAL(ICOOR1)
            END DO
C
            IF (ABS(GCCNST) .GT. DTHR) THEN
C
               DO ICOOR1 = 1, NCOOR
                  TRIAL(ICOOR1) = TRIAL(ICOOR1)
     &                          - RINVS*GCCNST*SYMCOR(ICOOR1,ICOOR)
               END DO
C
            END IF
C
C           *** Normalizing. ***
C
            RLENGT = D0
            DO ICOOR1 = 1, NCOOR
               RLENGT = RLENGT+ TRIAL(ICOOR1)**2
            END DO
            RLENGT = SQRT(RLENGT)
            IF (RLENGT.GT.DTHR) THEN
               RINV = D1/RLENGT
               DO ICOOR1 = 1, NCOOR
                  TRIAL(ICOOR1) = TRIAL(ICOOR1)*RINV
               END DO
            END IF
         END IF
      END DO
C
C     *** Has it become the zero vector?***
C
      RLENGT = D0
      DO ICOOR1 = 1, NCOOR
         RLENGT = RLENGT+ TRIAL(ICOOR1)**2
      END DO
      IF (RLENGT.LT.DTHR) EXIST = .TRUE.
C
C     *** Test printing. ***
C
      IF ((IPRINT .GT. 20) .AND. .NOT.EXIST) THEN
         WRITE (LUPRI,'(A)') 'Orthoganlized symmetry coordinate:'
C
         WRITE (LUPRI,'(A,2I4)') 'Symmetry', KIREP, KIPL-1
         DO I = 1, NCOOR
            WRITE (LUPRI,'(10X,F8.4)') TRIAL(I)
         END DO
      END IF
C
      RETURN
      END
C
C
C     /* Deck srtscr */
      SUBROUTINE SRTSCR(SYMCOR,TMPCOR,ICRIRP,ICTMP,IPRINT)
C     ******************************************************************
C     ***** This subroutine sorts the symmetry adapted coordinates *****
C     ***** such that all the coordinates that transforms as 1.    *****
C     ***** of a 2 dimensional irrep (for a given irrep), lies     *****
C     ***** before the coordinates that transforms as 2. row of    *****
C     ***** that irrep.                                            *****
C     ******************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <fcsym.h>
#include <trkoor.h>
      DIMENSION SYMCOR(NCOOR,NCOOR), TMPCOR(NCOOR,NCOOR)
      DIMENSION ICRIRP(NCOOR,2), ICTMP(NCOOR,2)
      LOGICAL FND2D
C
C     *** Copying all the coordinates belonging to 1. dim ireps ***
C     *** over to TMPCOR                                        ***
C
      IC = 0
      DO 200 IC2 = 1, NCOOR
         IF (ICRIRP(IC2,1) .LE. N1DIME) THEN
            IC = IC + 1
            ICTMP(IC,1) = ICRIRP(IC2,1)
            ICTMP(IC,2) = 0
            DO 300 IC1 = 1, NCOOR
               TMPCOR(IC1,IC2) = SYMCOR(IC1,IC2)
 300        CONTINUE
         END IF
 200  CONTINUE
      N1DCOR = IC
C
C     *** Sorting the rest of the coordinates according ***
C     *** to the principle above.                       ***
C
      IC = N1DCOR
      ITMPC = N1DCOR
      DO 400 IIREP = N1DIME+1, NIREP
C
         DO ITIM   = 0, 1
         DO ICOOR1 = N1DCOR, NCOOR
            IF ((ICRIRP(ICOOR1,1).EQ.IIREP).AND.
     &          (ICRIRP(ICOOR1,2).EQ.ITIM )) THEN
               IC = IC + 1
               ICTMP(IC,1) = IIREP
               ICTMP(IC,2) = ITIM
               DO ICOOR2 = 1, NCOOR
                  TMPCOR(ICOOR2,IC) = SYMCOR(ICOOR2,ICOOR1)
               END DO
            END IF
         END DO
         END DO
 400  CONTINUE
C
C     *** Test print ***
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Symmetry coordinates before sorting',-1)
         DO 700 IC1 = 1, NCOOR
            WRITE (LUPRI,'(12F12.8)') (SYMCOR(IC1,IC2),IC2=1,NCOOR)
 700     CONTINUE
      END IF
C
C     *** Assigning the new coordinate order to SYMCOR ***
C
      DO IC2 = 1, NCOOR
         ICRIRP(IC2,1) = ICTMP(IC2,1)
         ICRIRP(IC2,2) = ICTMP(IC2,2)
         DO IC1 = 1, NCOOR
            SYMCOR(IC1,IC2) = TMPCOR(IC1,IC2)
         END DO
      END DO
C
C     *** Test print ***
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('The rows are placed as follows',-1)
         WRITE (LUPRI,'(15I6)') (ICRIRP(IC,1), IC = 1, NCOOR)
         WRITE (LUPRI,'(15I6)') (ICRIRP(IC,2), IC = 1, NCOOR)
         CALL HEADER('Symmetry coordinates after sorting',-1)
         DO 900 IC1 = 1, NCOOR
            WRITE (LUPRI,'(12F15.8)') (SYMCOR(IC1,IC2),IC2=1,NCOOR)
 900     CONTINUE
      END IF
C
      RETURN
      END
C
C
C     /*Deck fsdcst*/
      SUBROUTINE FSDCST(SYMCOR,GRIREP,DCOEFF,WORK,KDPMTX,NMIDPC,ICRIRP,
     &                  LDPMTX,IFRSTD,NLDPMX,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
#include <numder.h>
#include <fcsym.h>
      DIMENSION SYMCOR(NCOOR,NCOOR), GRIREP(NGORDR,NGVERT),
     &          DCOEFF(LDPMTX,IFRSTD), WORK(LWORK)
      DIMENSION ICRIRP(NCOOR,2), KDPMTX(LDPMTX,NSTRDR,IFRSTD),
     &          NMIDPC(LDPMTX)
C
      call flshfo(5)
      DO 100 IORDR = 1, NMORDR
         IF (IORDR .EQ. 2) THEN
            CALL SCNFCS(DCOEFF,KDPMTX,ICRIRP,LDPMTX,IFRSTD,NLDPMX,
     &                  IPRINT)
C
            DO 200 II = 1, NLDPMX
               NMIDPC(II) = 1
 200        CONTINUE
         END IF
C
         IF (IORDR .EQ. 3) THEN
            KEQUMT = 1
            KIDXTS = KEQUMT + 2**(2*IORDR)
            KLAST  = KIDXTS + 8*NMORDR
            LWRK = LWORK - KLAST
            CALL THRDFC(DCOEFF,GRIREP,WORK(KEQUMT),WORK(KLAST),KDPMTX,
     &                  NMIDPC,ICRIRP,WORK(KIDXTS),LDPMTX,IFRSTD,NLDPMX,
     &                  LWRK,IPRINT)
         END IF
C
         IF (IORDR .EQ. 4) THEN
            KEQUMT = 1
            KIDXTS = KEQUMT + 2**(2*IORDR)
            KLAST  = KIDXTS + 16*NMORDR
            LWRK   = LWORK - KLAST
            CALL FRTHFC(DCOEFF,GRIREP,WORK(KEQUMT),WORK(KLAST),KDPMTX,
     &                  NMIDPC,ICRIRP,WORK(KIDXTS),LDPMTX,IFRSTD,NLDPMX,
     &                  LWRK,IPRINT)
         END IF
 100  CONTINUE
C
      RETURN
      END
C
C
C     /* Deck thrdfc*/
      SUBROUTINE THRDFC(DCOEFF,GRIREP,EQUMTX,WORK,KDPMTX,NMIDPC,ICRIRP,
     &                  IDXTST,LDPMTX,IFRSTD,NLDPMX,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0=0.0D0, D1=1.0D0, DMTHR = 1.0D-6)
#include <trkoor.h>
#include <numder.h>
#include <fcsym.h>
      LOGICAL DPNDCY, EXIST, IEXIST, TRIVIA
      DIMENSION GRIREP(NGORDR,NGVERT), DCOEFF(LDPMTX,IFRSTD),
     &          TCOEFF(8,8), EQUMTX(8,8),WORK(LWORK)
      DIMENSION KDPMTX(LDPMTX,NSTRDR,IFRSTD), NMIDPC(LDPMTX),
     &          ICRIRP(NCOOR,2),IDEP(8,8), IDXTST(NMORDR,8),
     &          ICRWMX(8,3)
C
      NLDPST = NLDPMX
C
      DO 100 IC3 = 1, NDCOOR
      DO 100 IC2 = 1, IC3
      DO 100 IC1 = 1, IC2
C
         IF (ICRIRP(IC1,1).GT.N1DIME) THEN
            ISTII = N1DIME + 4*(ICRIRP(IC1,1)-N1DIME-1)
            IF (ICRIRP(IC1,2).EQ.0) THEN
               N1STR =  1
               N1END =  2
               N1STP =  1
            ELSE
               N1STR =  2
               N1END =  1
               N1STP = -1
            END IF
         ELSE
            ISTII = ICRIRP(IC1,1)-1
            N1STR = 1
            N1END = 1
            N1STP = 1
         END IF
C
         IF (ICRIRP(IC2,1).GT.N1DIME) THEN
            ISTIJ = N1DIME + 4*(ICRIRP(IC2,1)-N1DIME-1)
            IF (ICRIRP(IC2,2).EQ.0) THEN
               N2STR =  1
               N2END =  2
               N2STP =  1
            ELSE
               N2STR =  2
               N2END =  1
               N2STP = -1
            END IF
         ELSE
            ISTIJ = ICRIRP(IC2,1)-1
            N2STR = 1
            N2END = 1
            N2STP = 1
         END IF
C
         IF (ICRIRP(IC3,1).GT.N1DIME) THEN
            ISTIK = N1DIME + 4*(ICRIRP(IC3,1)-N1DIME-1)
            IF (ICRIRP(IC3,2).EQ.0) THEN
               N3STR =  1
               N3END =  2
               N3STP =  1
            ELSE
               N3STR =  2
               N3END =  1
               N3STP = -1
            END IF
         ELSE
            ISTIK = ICRIRP(IC3,1)-1
            N3STR = 1
            N3END = 1
            N3STP = 1
         END IF
C
C        *** Number of coordinates that transform as the  ***
C        *** same row of the same irep (which is equal to ***
C        *** the spacing between function 1 and 2 of the  ***
C        *** same mode in a 2-dimensional irep)           ***
C
         N1ICOR = 0
         N2ICOR = 0
         N3ICOR = 0
         DO 200 IC = 1, NDCOOR
            IF ((ICRIRP(IC,1).EQ.ICRIRP(IC1,1)).AND.
     &          (ICRIRP(IC,2).EQ.            0)) THEN
               N1ICOR = N1ICOR + 1
            END IF
C
            IF ((ICRIRP(IC,1).EQ.ICRIRP(IC2,1)).AND.
     &          (ICRIRP(IC,2).EQ.            0)) THEN
               N2ICOR = N2ICOR + 1
            END IF
C
            IF ((ICRIRP(IC,1).EQ.ICRIRP(IC3,1)).AND.
     &          (ICRIRP(IC,2).EQ.            0)) THEN
               N3ICOR = N3ICOR + 1
            END IF
 200     CONTINUE
C
         KDIM = 64
         CALL DZERO(EQUMTX,KDIM)
C
         ICIJK = 0
         DO 300 IK = N3STR, N3END, N3STP
         DO 300 IJ = N2STR, N2END, N2STP
         DO 300 II = N1STR, N1END, N1STP
            ICIJK = ICIJK + 1
C
            ICRWMX(ICIJK,1) = IC3 - ICRIRP(IC3,2)*N3ICOR + (IK-1)*N3ICOR
            ICRWMX(ICIJK,2) = IC2 - ICRIRP(IC2,2)*N2ICOR + (IJ-1)*N2ICOR
            ICRWMX(ICIJK,3) = IC1 - ICRIRP(IC1,2)*N1ICOR + (II-1)*N1ICOR
C
            ICTUV = 0
            DO 400 IV = N3STR, N3END, N3STP
            DO 400 IU = N2STR, N2END, N2STP
            DO 400 IT = N1STR, N1END, N1STP
               ICTUV = ICTUV + 1
C
               IV1 = ISTII + 2*(II-1) + IT
               IV2 = ISTIJ + 2*(IJ-1) + IU
               IV3 = ISTIK + 2*(IK-1) + IV
               DO 500 IGORDR = 1, NGORDR
                  EQUMTX(ICTUV,ICIJK) = EQUMTX(ICTUV,ICIJK)
     &                                +(GRIREP(IGORDR,IV1)
     &                                 *GRIREP(IGORDR,IV2)
     &                                 *GRIREP(IGORDR,IV3))/DBLE(NGORDR)
 500           CONTINUE
 400        CONTINUE
 300     CONTINUE
C
         DPNDCY = .FALSE.
         IF (ICIJK.GT.1) THEN
            SUMMTX = D0
            DO 600 IIIJK = 1, ICIJK
            DO 600 IITUV = 1, ICTUV
               SUMMTX = SUMMTX + ABS(EQUMTX(IITUV,IIIJK))
 600        CONTINUE
            IF (SUMMTX .GT. DMTHR) DPNDCY = .TRUE.
         END IF
C
         IF (DPNDCY) THEN
            DO 700 II = 1, ICIJK
               EQUMTX(II,II) = EQUMTX(II,II)-D1
 700        CONTINUE
C
            IF (IPRINT .GT. 20) THEN
               WRITE (LUPRI,'(A)') 'Matrix to determine dependent' //
     &                             ' force constants'
               WRITE (LUPRI,'(2X,A,I4)') 'Number of components',
     &                   (N1STR+N1END-1)*(N2STR+N2END-1)*(N3STR+N3END-1)
               WRITE (LUPRI,'(2X,A,3I4)') 'Component', IC1, IC2, IC3
               DO IIIJK = 1, ICIJK
                  WRITE (LUPRI,'(8F8.4)')
     &                  (EQUMTX(IITUV,IIIJK),IITUV=1,ICTUV)
               END DO
               WRITE (LUPRI,'(A)') '                                   '
               WRITE (LUPRI,'(A)') '                                   '
            END IF
C
C           *** Sorting the components, such that the force-constants ***
C           *** calculated is situated last.                          ***
C
            KITMPE = 1
            KITMPR = KITMPE + 8
            CALL SRTNCF(EQUMTX,WORK(KITMPR),ICRWMX,WORK(KITMPE),8,3,
     &                  ICIJK,IPRINT)
            ICTUV = ICIJK
C
C           *** Solving the homogeneous linear equation system. ***
C
            KIRNDX = 1
            KZINDX = KIRNDX + 8
            KIDXTP = KZINDX + 8
            CALL DIAGUD(EQUMTX,TCOEFF,IDEP,WORK(KIRNDX),WORK(KZINDX),
     &                  WORK(KIDXTP),8,ICIJK,NROW,NIDEP,IPRINT)
C
C           *** Independent component.                                    ***
C           *** Checking if it is av valid component, and assigning the   ***
C           *** "coordinate-numbers"                                      ***
C
            DO 800 IDP = 1, NIDEP
               DO 900 ILNGTH = 1, ICIJK
                  IF (ILNGTH.EQ.IDEP(1,IDP+1)) THEN
                     IDXTST(1,IDP) = ICRWMX(IDEP(1,IDP+1),1)
                     IDXTST(2,IDP) = ICRWMX(IDEP(1,IDP+1),2)
                     IDXTST(3,IDP) = ICRWMX(IDEP(1,IDP+1),3)
C
C                    *** The other reason for discarding the batch   ***
C                    *** based on the independent component, is that ***
C                    *** the batch already EXIST                     ***
C
                     CALL CHKCMP(KDPMTX,IDXTST(1,IDP),3,LDPMTX,IFRSTD,
     &                           NLDPST,NLDPMX,NIDEP,IPRINT,EXIST)
                  END IF
 900           CONTINUE
 800        CONTINUE
C
C           *** Finding the dependent components ***
C
            NSTART = NLDPMX
            IF (.NOT.EXIST) THEN
               DO 1100 IROW = 1, NROW-1
                  DO 1200 ILNGTH = 1, ICIJK
                     IF (ILNGTH.EQ.IDEP(IROW,1)) THEN
C
C                       *** Checking whether there are several equal solutions for  ***
C                       *** one dependent force constant.                           ***
C
                        IDP3 = ICRWMX(IDEP(IROW,1),1)
                        IDP2 = ICRWMX(IDEP(IROW,1),2)
                        IDP1 = ICRWMX(IDEP(IROW,1),3)
                        IEXIST = .FALSE.
                        DO 1300 ICOUNT = NSTART+1, NLDPMX
                           IF (.NOT.IEXIST) THEN
                              IEXIST = (IDP3.EQ.KDPMTX(ICOUNT,1,1)).AND.
     &                                 (IDP2.EQ.KDPMTX(ICOUNT,2,1)).AND.
     &                                 (IDP1.EQ.KDPMTX(ICOUNT,3,1))
                           END IF
 1300                   CONTINUE
C
C                       *** Checking if this is a trivial solution of the dependent ***
C                       *** force constant of the type K_{aaa} = K_{aaa}            ***
C
                        IF (NIDEP.EQ.1) THEN
                           TRIVIA = (IDP3.EQ.IDXTST(1,1)).AND.
     &                              (IDP2.EQ.IDXTST(2,1)).AND.
     &                              (IDP1.EQ.IDXTST(3,1))
                        END IF
C
                        IF ((.NOT.IEXIST).AND.(.NOT.TRIVIA)) THEN
                           NLDPMX = NLDPMX + 1
                           KDPMTX(NLDPMX,1,1) = IDP3
                           KDPMTX(NLDPMX,2,1) = IDP2
                           KDPMTX(NLDPMX,3,1) = IDP1
                           NMIDPC(NLDPMX        ) = NIDEP
                           DO 1400 IDP = 1, NIDEP
                              KDPMTX(NLDPMX,1,IDP+1) = IDXTST(1,IDP)
                              KDPMTX(NLDPMX,2,IDP+1) = IDXTST(2,IDP)
                              KDPMTX(NLDPMX,3,IDP+1) = IDXTST(3,IDP)
                              DCOEFF(NLDPMX  ,IDP  ) = TCOEFF(IROW,IDP)
 1400                      CONTINUE
                        ELSE
                           NROW = NROW - 1
                        END IF
                     END IF
 1200             CONTINUE
 1100          CONTINUE
C
C
C
               KIDXTP = 1
               KNMIDP = KIDXTP + 3
               KMVTMP = KNMIDP + 3
               KICMPI = KMVTMP + 3
               CALL CHKLCP(DCOEFF,KDPMTX,ICRIRP,WORK(KIDXTP),
     &                     WORK(KNMIDP),WORK(KMVTMP),WORK(KICMPI),
     &                     LDPMTX,IFRSTD,3,NROW-1,NLDPMX-NROW+NIDEP+1,
     &                     NIDEP,IPRINT)
C
C              *** Test print ***
C
               CALL HEADER('Coordinate dependency in cartesians.',-1)
               WRITE (LUPRI,'(5X,A,I4)') 'Number of independent' //
     &                 'force constants.', NMIDPC(NLDPMX)
               WRITE (LUPRI,'(5X,A)') 'Dependent components      ' //
     &                          'Independent components     Coeffisient'
               DO II = NLDPMX-NROW+2, NLDPMX
                 WRITE (LUPRI,'(A,3I4,A,3I4,A,F8.4)') '          ',
     &                 (KDPMTX(II,IX,1),IX=1,3), '          ',
     &                 (KDPMTX(II,IX,2),IX=1,3), '              ',
     &                  DCOEFF(II,1)
                  DO IDP = 2, NIDEP
                     WRITE (LUPRI,'(A,3I4,A,F8.4)')
     &                      '                                ',
     &                  (KDPMTX(II,IX,IDP+1),IX=1,3), '               ',
     &                  DCOEFF(II,IDP)
                  END DO
               END DO
               WRITE (LUPRI,'(A)')'                                 '
               WRITE (LUPRI,'(A)')'                                 '
            ELSE
               IF (IPRINT.GT.20) THEN
                  WRITE (LUPRI,'(A)') '                                '
                  WRITE (LUPRI,'(A)') '                                '
                  WRITE (LUPRI,'(A)')
     &                             'Component not going through due to:'
                  WRITE (LUPRI,'(A,L1)')'Exist               : ', EXIST
                  WRITE (LUPRI,'(A)') '                                '
                  WRITE (LUPRI,'(A)') '                                '
               END IF
            END IF
         END IF
C
 100  CONTINUE
C
      RETURN
      END
C
C
C     /* Deck scnfcs */
      SUBROUTINE SCNFCS(DCOEFF,KDPMTX,ICRIRP,LDPMTX,IFRSTD,NLDPMX,
     &                  IPRINT)
C     ***************************************************
C     **** Figures out which second derivatives are  ****
C     **** symmetrical dependent. The dependent      ****
C     **** components are put in first row of KDPMTX ****
C     **** The independent are put in the second row ****
C     **** The dependency coeffecients are put in    ****
C     **** DCOEFF.                                   ****
C     ***************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D1 = 1.0D0)
#include <trkoor.h>
#include <numder.h>
#include <fcsym.h>
      DIMENSION DCOEFF(LDPMTX,IFRSTD), KDPMTX(LDPMTX,NSTRDR,IFRSTD),
     &          ICRIRP(NCOOR,2)
C
C     *** Loop over all the 2 dimensional ireps ***
C
      KDPSTR = NLDPMX+1
      DO 100 IIRP = N1DIME+1, NCVERT
C
C        *** Find how many functions there are of this irrep ***
C
         INUM = 0
         DO 200 IC = 1, NDCOOR
            IF ((ICRIRP(IC,1).EQ.IIRP).AND.(ICRIRP(IC,2).EQ.0)) THEN
               IF (INUM.EQ.0) ISTART = IC
               INUM = INUM + 1
            END IF
 200     CONTINUE
C
C        *** Assign dependencies ***
C
         DO 300 IC2 = ISTART, ISTART+INUM-1
         DO 300 IC1 = ISTART, IC2
            NLDPMX = NLDPMX + 1
            KDPMTX(NLDPMX,1,1) = IC2 + INUM
            KDPMTX(NLDPMX,2,1) = IC1 + INUM
            KDPMTX(NLDPMX,1,2) = IC2
            KDPMTX(NLDPMX,2,2) = IC1
            DCOEFF(NLDPMX,1  ) = D1
 300     CONTINUE
C
 100  CONTINUE
C
C     *** Test print ***
C
      CALL HEADER('Dependent and indepenent force coefficients',-1)
      WRITE (LUPRI,'(5X,A)') 'Dependent components   Independent ' //
     &                          'components     coefficients'
      DO II = KDPSTR, NLDPMX
         WRITE (LUPRI,'(I15,I4,I19,I4,15X,F6.2)')
     &        (KDPMTX(II,J,1),J=1,2), (KDPMTX(II,J,2),J=1,2),
     &        DCOEFF(II,1)
      END DO
C
      RETURN
      END
C
C
C     /* Deck frthfc*/
      SUBROUTINE FRTHFC(DCOEFF,GRIREP,EQUMTX,WORK,KDPMTX,NMIDPC,ICRIRP,
     &                  IDXTST,LDPMTX,IFRSTD,NLDPMX,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0=0.0D0, D1=1.0D0, DMTHR = 1.0D-6)
      PARAMETER (ITHDDM = 16)
#include <trkoor.h>
#include <numder.h>
#include <fcsym.h>
      LOGICAL DPNDCY, EXIST, IEXIST, TRIVIA
      DIMENSION GRIREP(NGORDR,NGVERT), DCOEFF(LDPMTX,IFRSTD),
     &          TCOEFF(ITHDDM,ITHDDM), EQUMTX(ITHDDM,ITHDDM),WORK(10000)
      DIMENSION KDPMTX(LDPMTX,NSTRDR,IFRSTD), NMIDPC(LDPMTX),
     &          ICRIRP(NCOOR,2),IDEP(ITHDDM,ITHDDM),
     &          IDXTST(NMORDR,ITHDDM), ICRWMX(ITHDDM,4), ITLRGR(ITHDDM)
C
c      iprint = 55
      NLDPST = NLDPMX
      DO 100 IC4 = 1, NDCOOR
      DO 100 IC3 = 1, IC4
      DO 100 IC2 = 1, IC3
      DO 100 IC1 = 1, IC2
C
         IF (ICRIRP(IC1,1).GT.N1DIME) THEN
            ISTII = N1DIME + 4*(ICRIRP(IC1,1)-N1DIME-1)
            IF (ICRIRP(IC1,2).EQ.0) THEN
               N1STR =  1
               N1END =  2
               N1STP =  1
            ELSE
               N1STR =  2
               N1END =  1
               N1STP = -1
            END IF
         ELSE
            ISTII = ICRIRP(IC1,1)-1
            N1STR = 1
            N1END = 1
            N1STP = 1
         END IF
C
         IF (ICRIRP(IC2,1).GT.N1DIME) THEN
            ISTIJ = N1DIME + 4*(ICRIRP(IC2,1)-N1DIME-1)
            IF (ICRIRP(IC2,2).EQ.0) THEN
               N2STR =  1
               N2END =  2
               N2STP =  1
            ELSE
               N2STR =  2
               N2END =  1
               N2STP = -1
            END IF
         ELSE
            ISTIJ = ICRIRP(IC2,1)-1
            N2STR = 1
            N2END = 1
            N2STP = 1
         END IF
C
         IF (ICRIRP(IC3,1).GT.N1DIME) THEN
            ISTIK = N1DIME + 4*(ICRIRP(IC3,1)-N1DIME-1)
            IF (ICRIRP(IC3,2).EQ.0) THEN
               N3STR =  1
               N3END =  2
               N3STP =  1
            ELSE
               N3STR =  2
               N3END =  1
               N3STP = -1
            END IF
         ELSE
            ISTIK = ICRIRP(IC3,1)-1
            N3STR = 1
            N3END = 1
            N3STP = 1
         END IF
C
         IF (ICRIRP(IC4,1).GT.N1DIME) THEN
            ISTIL = N1DIME + 4*(ICRIRP(IC4,1)-N1DIME-1)
            IF (ICRIRP(IC4,2).EQ.0) THEN
               N4STR =  1
               N4END =  2
               N4STP =  1
            ELSE
               N4STR =  2
               N4END =  1
               N4STP = -1
            END IF
         ELSE
            ISTIL = ICRIRP(IC4,1)-1
            N4STR = 1
            N4END = 1
            N4STP = 1
         END IF
C
C        *** Number of coordinates that transform as the  ***
C        *** same row of the same irep (which is equal to ***
C        *** the spacing between function 1 and 2 of the  ***
C        *** same mode in a 2-dimensional irep)           ***
C
         N1ICOR = 0
         N2ICOR = 0
         N3ICOR = 0
         N4ICOR = 0
         DO 200 IC = 1, NDCOOR
            IF ((ICRIRP(IC,1).EQ.ICRIRP(IC1,1)).AND.
     &          (ICRIRP(IC,2).EQ.            0)) THEN
               N1ICOR = N1ICOR + 1
            END IF
C
            IF ((ICRIRP(IC,1).EQ.ICRIRP(IC2,1)).AND.
     &          (ICRIRP(IC,2).EQ.            0)) THEN
               N2ICOR = N2ICOR + 1
            END IF
C
            IF ((ICRIRP(IC,1).EQ.ICRIRP(IC3,1)).AND.
     &          (ICRIRP(IC,2).EQ.            0)) THEN
               N3ICOR = N3ICOR + 1
            END IF
C
            IF ((ICRIRP(IC,1).EQ.ICRIRP(IC4,1)).AND.
     &          (ICRIRP(IC,2).EQ.            0)) THEN
               N4ICOR = N4ICOR + 1
            END IF
 200     CONTINUE
C
         KDIM = ITHDDM**2
         CALL DZERO(EQUMTX,KDIM)
C
         ICIJKL = 0
         DO 300 IL = N4STR, N4END, N4STP
         DO 300 IK = N3STR, N3END, N3STP
         DO 300 IJ = N2STR, N2END, N2STP
         DO 300 II = N1STR, N1END, N1STP
            ICIJKL = ICIJKL + 1
C
            ICRWMX(ICIJKL,1)= IC4 - ICRIRP(IC4,2)*N4ICOR + (IL-1)*N4ICOR
            ICRWMX(ICIJKL,2)= IC3 - ICRIRP(IC3,2)*N3ICOR + (IK-1)*N3ICOR
            ICRWMX(ICIJKL,3)= IC2 - ICRIRP(IC2,2)*N2ICOR + (IJ-1)*N2ICOR
            ICRWMX(ICIJKL,4)= IC1 - ICRIRP(IC1,2)*N1ICOR + (II-1)*N1ICOR

C
            ICTUVX = 0
            DO 400 IX = N4STR, N4END, N4STP
            DO 400 IV = N3STR, N3END, N3STP
            DO 400 IU = N2STR, N2END, N2STP
            DO 400 IT = N1STR, N1END, N1STP
               ICTUVX = ICTUVX + 1
C
               IV1 = ISTII + 2*(II-1) + IT
               IV2 = ISTIJ + 2*(IJ-1) + IU
               IV3 = ISTIK + 2*(IK-1) + IV
               IV4 = ISTIL + 2*(IL-1) + IX
               DO 500 IGORDR = 1, NGORDR
                  EQUMTX(ICTUVX,ICIJKL) = EQUMTX(ICTUVX,ICIJKL)
     &                                +(GRIREP(IGORDR,IV1)
     &                                 *GRIREP(IGORDR,IV2)
     &                                 *GRIREP(IGORDR,IV3)
     &                                 *GRIREP(IGORDR,IV4))/DBLE(NGORDR)
 500           CONTINUE
 400        CONTINUE
 300     CONTINUE
C
         DPNDCY = .FALSE.
         IF (ICIJKL.GT.1) THEN
            SUMMTX = D0
            DO 600 IIIJK = 1, ICIJKL
            DO 600 IITUV = 1, ICTUVX
               SUMMTX = SUMMTX + ABS(EQUMTX(IITUV,IIIJK))
 600        CONTINUE
            IF (SUMMTX .GT. DMTHR) DPNDCY = .TRUE.
         END IF
C
         IF (DPNDCY) THEN
            DO 700 II = 1, ICIJKL
               EQUMTX(II,II) = EQUMTX(II,II)-D1
 700        CONTINUE
C
            IF (IPRINT .GT. 20) THEN
               WRITE (LUPRI,'(A)') 'Matrix to determine dependent' //
     &                             ' force constants'
               WRITE (LUPRI,'(2X,A,I4)') 'Number of components',
     &                                 (N1STR+N1END-1)*(N2STR+N2END-1)
     &                                *(N3STR+N3END-1)*(N4STR+N4END-1)
               WRITE (LUPRI,'(2X,A,4I4)') 'Component', IC1,IC2,IC3,IC4
               DO IIIJKL = 1, ICIJKL
                  WRITE (LUPRI,'(16F8.4)')
     &                  (EQUMTX(IITUVX,IIIJKL),IITUVX=1,ICTUVX)
               END DO
               WRITE (LUPRI,'(A)') '                                   '
               WRITE (LUPRI,'(A)') '                                   '
            END IF
C
C           *** Sorting the components, such thet the highest number ***
C           *** component is first. In addition matrix elements for  ***
C           *** equal force-constants are added together             ***
C
            KITMPE = 1
            KITMPR = KITMPE + ITHDDM
            CALL SRTNCF(EQUMTX,WORK(KITMPR),ICRWMX,WORK(KITMPE),ITHDDM,
     &                  4,ICIJKL,IPRINT)
C
C           *** Solving the homogeneous linear equation system. ***
C
            KIRNDX = 1
            KZINDX = KIRNDX + ICIJKL
            KIDXTP = KZINDX + ICIJKL
            CALL DIAGUD(EQUMTX,TCOEFF,IDEP,WORK(KIRNDX),WORK(KZINDX),
     &                  WORK(KIDXTP),ITHDDM,ICIJKL,NROW,NIDEP,IPRINT)
C
C           *** Independent component. ***
C
            DO 800 IDP = 1, NIDEP
               DO 900 ILNGTH = 1, ICIJKL
                  IF (ILNGTH.EQ.IDEP(1,IDP+1)) THEN
                     IDXTST(1,IDP) = ICRWMX(IDEP(1,IDP+1),1)
                     IDXTST(2,IDP) = ICRWMX(IDEP(1,IDP+1),2)
                     IDXTST(3,IDP) = ICRWMX(IDEP(1,IDP+1),3)
                     IDXTST(4,IDP) = ICRWMX(IDEP(1,IDP+1),4)
C
C                    *** One reason of reason for discarding the batch ***
C                    *** based on the independent component, is that   ***
C                    *** the batch already EXIST                       ***
C
                     CALL CHKCMP(KDPMTX,IDXTST(1,IDP),4,LDPMTX,IFRSTD,
     &                           NLDPST,NLDPMX,NIDEP,IPRINT,EXIST)
                  END IF
 900           CONTINUE
 800        CONTINUE
C
C           *** Finding the dependent components ***
C
            NSTART = NLDPMX
            IF (.NOT.EXIST) THEN
               DO 1100 IROW = 1, NROW-NIDEP
                  DO 1200 ILNGTH = 1, ICIJKL
                     IF (ILNGTH.EQ.IDEP(IROW,1)) THEN
C
C                       *** Checking whether there are several equal    ***
C                       *** solutions for one dependent force constant. ***
C
C
                        IDP4 = ICRWMX(IDEP(IROW,1),1)
                        IDP3 = ICRWMX(IDEP(IROW,1),2)
                        IDP2 = ICRWMX(IDEP(IROW,1),3)
                        IDP1 = ICRWMX(IDEP(IROW,1),4)
                        IEXIST = .FALSE.
                        DO 1300 ICOUNT = NSTART+1, NLDPMX
                           IF (.NOT.IEXIST) THEN
                              IEXIST = (IDP4.EQ.KDPMTX(ICOUNT,1,1)).AND.
     &                                 (IDP3.EQ.KDPMTX(ICOUNT,2,1)).AND.
     &                                 (IDP2.EQ.KDPMTX(ICOUNT,3,1)).AND.
     &                                 (IDP1.EQ.KDPMTX(ICOUNT,4,1))
                           END IF
 1300                   CONTINUE
C
C                       *** Checking if this is a trivial solution of the dependent ***
C                       *** force constant of the type K_{abcd} = K_{abcd}          ***
C
                        IF (NIDEP.EQ.1) THEN
                           TRIVIA = (IDP4.EQ.IDXTST(1,1)).AND.
     &                              (IDP3.EQ.IDXTST(2,1)).AND.
     &                              (IDP2.EQ.IDXTST(3,1)).AND.
     &                              (IDP1.EQ.IDXTST(4,1))
                        END IF
C
                        IF ((.NOT.IEXIST).AND.(.NOT.TRIVIA)) THEN
                           NLDPMX = NLDPMX + 1
                           KDPMTX(NLDPMX,1,1) = IDP4
                           KDPMTX(NLDPMX,2,1) = IDP3
                           KDPMTX(NLDPMX,3,1) = IDP2
                           KDPMTX(NLDPMX,4,1) = IDP1
                           NMIDPC(NLDPMX    ) = NIDEP
                           write (lupri,*) 'hallo', NLDPMX
                           DO 1400 IDP = 1, NIDEP
                              KDPMTX(NLDPMX,1,IDP+1) = IDXTST(1,IDP)
                              KDPMTX(NLDPMX,2,IDP+1) = IDXTST(2,IDP)
                              KDPMTX(NLDPMX,3,IDP+1) = IDXTST(3,IDP)
                              KDPMTX(NLDPMX,4,IDP+1) = IDXTST(4,IDP)
                              DCOEFF(NLDPMX  ,IDP  ) = TCOEFF(IROW,IDP)
 1400                      CONTINUE
                        ELSE
                           NROW = NROW - 1
                        END IF
                     END IF
 1200             CONTINUE
 1100          CONTINUE
C
C              *** Checking whether the calulated independent components  ***
C              *** are the the components with least energy calculations. ***
C
               KIDXTP = 1
               KNMIDP = KIDXTP + 4
               KMVTMP = KNMIDP + 4
               KICMPI = KMVTMP + 4
               CALL CHKLCP(DCOEFF,KDPMTX,ICRIRP,WORK(KIDXTP),
     &                     WORK(KNMIDP),WORK(KMVTMP),WORK(KICMPI),
     &                     LDPMTX,IFRSTD,4,NROW-1,NLDPMX-NROW+NIDEP+1,
     &                     NIDEP,IPRINT)
C
C              *** Test print ***
C
               IF (NROW.GT.1) THEN
                  CALL HEADER('Coordinate dependency in cartesians.',-1)
                  WRITE (LUPRI,'(5X,A,I4)') 'Number of independent' //
     &                 'force constants.', NMIDPC(NLDPMX)
                  WRITE (LUPRI,'(5X,A)') 'Dependent components      ' //
     &                          'Independent components     Coeffisient'
                  DO II = NLDPMX-NROW+NIDEP+1, NLDPMX
                     write (lupri,*) 'component', ii
                     WRITE (LUPRI,'(A,4I4,A,4I4,A,F8.4)')
     &                   '        ', (KDPMTX(II,IX,1),IX=1,4), '      ',
     &                    (KDPMTX(II,IX,2),IX=1,4), '             ',
     &                     DCOEFF(II,1)
                     DO IDP = 2, NIDEP
                        WRITE (LUPRI,'(A,4I4,A,F8.4)')
     &                      '                              ',
     &                    (KDPMTX(II,IX,IDP+1),IX=1,4), '             ',
     &                     DCOEFF(II,IDP)
                     END DO
                  END DO
                  WRITE (LUPRI,'(A)')'                                 '
                  WRITE (LUPRI,'(A)')'                                 '
               ELSE
                  WRITE (LUPRI,'(A)') '                                '
                  WRITE (LUPRI,'(A)') '                                '
                  WRITE (LUPRI,'(5X,A)') 'No valid component of the' //
     &                                   ' independent force constant'
                  WRITE (LUPRI,'(A)') '                                '
                  WRITE (LUPRI,'(A)') '                                '
               END IF
            ELSE
               IF (IPRINT.GT.20) THEN
                  WRITE (LUPRI,'(A)') '                                '
                  WRITE (LUPRI,'(A)') '                                '
                  WRITE (LUPRI,'(A)')
     &                             'Component not going through due to:'
                  WRITE (LUPRI,'(A,L1)')'Exist               : ', EXIST
                  WRITE (LUPRI,'(A)') '                                '
                  WRITE (LUPRI,'(A)') '                                '
               END IF
            END IF
         END IF
C
 100  CONTINUE
      iprint = 0
C
      RETURN
      END
C
C
C     /*Deck addpfc*/
      SUBROUTINE ADDPFC(DERIV,DCOEFF,KDPMTX,NMIDPC,LDPMTX,IFRSTD,NDERIV,
     &                  NLDPMX,IPRINT)
C     *******************************************************************
C     **** This assigns values to the dependent force-constants from ****
C     **** the independent.                                          ****
C     *******************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
#include <taymol.h>
#include <numder.h>
#include <fcsym.h>
      DIMENSION DERIV(NDERIV), DCOEFF(LDPMTX,IFRSTD)
      DIMENSION KDPMTX(LDPMTX,NSTRDR,IFRSTD), NMIDPC(LDPMTX)
C
      NSTART = 1
      DO 100 IRDR = 1, NMORDR
         IF (IRDR .EQ. 2) THEN
            ISTART = 0
            MINTST = MIN(NMORDR,3)
            DO 200 ILDP = 1, NLDPMX
               IF ((KDPMTX(ILDP,MINTST,1).EQ.0).OR.
     &                            (MINTST.EQ.2)) THEN
                  HESMOL(KDPMTX(ILDP,1,1),KDPMTX(ILDP,2,1)) =
     &                         HESMOL(KDPMTX(ILDP,1,2),KDPMTX(ILDP,2,2))
                  ISTART = ISTART + 1
               END IF
 200        CONTINUE
            NSTART = NSTART + ISTART
         ELSE IF (IRDR.EQ.3) THEN
            ISTART = 0
            MINTST = MIN(NMORDR,4)
            DO 300 ILDP = NSTART, NLDPMX
               IF (((KDPMTX(ILDP,MINTST,1).EQ.0).AND.(MINTST.GT.3))
     &            .OR.((MINTST.EQ.3).AND.(KDPMTX(ILDP,3,1).NE.0))) THEN
C
                  ISTART = ISTART + 1
C
                  IDPDRV = ((KDPMTX(ILDP,1,1)-1)*KDPMTX(ILDP,1,1)
     &                    * (KDPMTX(ILDP,1,1)+1))/6
     &                    +((KDPMTX(ILDP,2,1)-1)*KDPMTX(ILDP,2,1))/2
     &                    +  KDPMTX(ILDP,3,1)
                  DERIV(IDPDRV) = D0
                  DO 400 II = 1, NMIDPC(ILDP)
                     IIDPDR =
     &                   ((KDPMTX(ILDP,1,II+1)-1)*KDPMTX(ILDP,1,II+1)
     &                  * (KDPMTX(ILDP,1,II+1)+1))/6
     &                  +((KDPMTX(ILDP,2,II+1)-1)*KDPMTX(ILDP,2,II+1))/2
     &                  +  KDPMTX(ILDP,3,II+1)
                     DERIV(IDPDRV) = DERIV(IDPDRV)
     &                             + DCOEFF(ILDP,II)*DERIV(IIDPDR)
 400              CONTINUE
               END IF
 300        CONTINUE
            NSTART = NSTART + ISTART
         ELSE IF (IRDR.EQ.4) THEN
            MINTST = MIN(NMORDR,5)
            IOFFST =(NDCOOR*(NDCOOR+1)*(NDCOOR+2))/6
C
            DO 500 ILDP = NSTART, NLDPMX
               IF (((KDPMTX(ILDP,MINTST,1).EQ.0).AND.(MINTST.GT.4))
     &             .OR.((MINTST.EQ.4).AND.(KDPMTX(ILDP,4,1).NE.0))) THEN
                  IDPDRV = IOFFST +
     &                    ((KDPMTX(ILDP,1,1)-1)* KDPMTX(ILDP,1,1)
     &                    *(KDPMTX(ILDP,1,1)+1)*(KDPMTX(ILDP,1,1)+2))/24
     &                   +((KDPMTX(ILDP,2,1)-1)*KDPMTX(ILDP,2,1)
     &                   * (KDPMTX(ILDP,2,1)+1))/6
     &                   +((KDPMTX(ILDP,3,1)-1)*KDPMTX(ILDP,3,1))/2
     &                   +  KDPMTX(ILDP,4,1)
                  DERIV(IDPDRV) = D0
                  DO 600 II = 1, NMIDPC(ILDP)
                     IJ = II+1
                     IIDPDR = IOFFST +
     &                  ((KDPMTX(ILDP,1,IJ)-1)* KDPMTX(ILDP,1,IJ)
     &                  *(KDPMTX(ILDP,1,IJ)+1)*(KDPMTX(ILDP,1,IJ)+2))/24
     &                 +((KDPMTX(ILDP,2,IJ)-1)* KDPMTX(ILDP,2,IJ)
     &                  *(KDPMTX(ILDP,2,IJ)+1))/6
     &                 +((KDPMTX(ILDP,3,IJ)-1)*KDPMTX(ILDP,3,IJ))/2
     &                 +  KDPMTX(ILDP,4,IJ)
                     DERIV(IDPDRV) = DERIV(IDPDRV)
     &                             + DCOEFF(ILDP,II)*DERIV(IIDPDR)
 600              CONTINUE
               END IF
 500        CONTINUE
         END IF
 100  CONTINUE
C
      RETURN
      END
C
C
C     /* Deck diagud */
      SUBROUTINE DIAGUD(RMTRIX,COEFF,IDEP,IRINDX,IZINDX,IDXTMP,NRDIM,
     &                  NDIM,NROW,NIDEP,IPRINT)
C     ******************************************************************
C     *** Subroutine that solves a underdetermined set of homogenous ***
C     *** linear equations.                                          ***
C     *** In : RMTRIX - coeffisient matrix to the equation set.      ***
C     *** Out: NIDEP -> Number of independent variables              ***
C     ***      IDEP  -> First row: Dependendent variables            ***
C     ***               Scn.  row: Dependend on this variable        ***
C     ***      COEFF-> Dependence coeffisient.                      ***
C     ******************************************************************
#include <implicit.h>
#include <priunit.h>
      LOGICAL EXIT
      PARAMETER (D0 = 0.0D0, D1=1.0D0, DMTHR = 1.0D-6)
C
      DIMENSION RMTRIX(NRDIM,NRDIM), COEFF(NRDIM,NRDIM)
      DIMENSION IDEP(NRDIM,NRDIM), IRINDX(NRDIM), IZINDX(NRDIM),
     &          IDXTMP(NRDIM,2)
C
      CALL RMZROR(RMTRIX,IRINDX,IZINDX,NRDIM,NDIM,NROW,IPRINT)
C
C     *** The matrix is then "Gaussian-eliminated" to get it ***
C     *** almost on diagonal form.                           ***
C
C     *** Forward substitution. ***
C
      DO 100 ICOL = 1   , NROW
         IF (ABS(RMTRIX(ICOL,ICOL)) .GT. DMTHR) THEN
            BCOEFF = D1/RMTRIX(ICOL,ICOL)
C
            DO 200 IROW = ICOL+1, NROW
               TCOEFF = RMTRIX(IROW,ICOL)*BCOEFF
               DO 300 IICOL = ICOL, NROW
                  RMTRIX(IROW,IICOL) = RMTRIX(IROW,IICOL)
     &                        - TCOEFF*RMTRIX(ICOL,IICOL)
 300           CONTINUE
 200        CONTINUE
C
         END IF
 100  CONTINUE
C
C     *** Back substitution. ***
C
      DO 400 ICOL = NROW, 1, -1
         IF (ABS(RMTRIX(ICOL,ICOL)) .GT. DMTHR) THEN
            BCOEFF = D1/RMTRIX(ICOL,ICOL)
            DO 500 IROW = 1, ICOL-1
               TCOEFF = RMTRIX(IROW,ICOL)*BCOEFF
               DO 600 IICOL = ICOL, NROW
                  RMTRIX(IROW,IICOL) = RMTRIX(IROW,IICOL)
     &                        - TCOEFF*RMTRIX(ICOL,IICOL)
 600           CONTINUE
 500        CONTINUE
         END IF
 400  CONTINUE
C
C     *** Each row is then normalized **
C
      DO 700 IROW = 1, NROW
         IF (ABS(RMTRIX(IROW,IROW)).GT.DMTHR) THEN
            BCOEFF = D1/RMTRIX(IROW,IROW)
            DO 800 ICOL = IROW, NROW
               RMTRIX(IROW,ICOL) = BCOEFF*RMTRIX(IROW,ICOL)
 800        CONTINUE
         END IF
 700  CONTINUE
C
C     *** Finding rows in the matrix that are zero ***
C     *** These are the independent variables.     ***
C
      IZRO = 0
      IDP = 0
      EXIT = .FALSE.
      DO 900 IROW = 1, NROW
         SUMROW = D0
         DO 1100 ICOL = 1, NROW
            SUMROW = SUMROW + ABS(RMTRIX(IROW,ICOL))
 1100    CONTINUE
         IF (SUMROW.LT.DMTHR) THEN
            IZRO = IZRO + 1
            IDXTMP(IZRO,2) = IROW
         ELSE
            IDP = IDP + 1
            IDXTMP(IDP,1) = IROW
            IDEP(IDP,1)   = IRINDX(IROW)
         END IF
 900  CONTINUE
      NIDEP = IZRO
      NDDEP = IDP
C
C     *** Assigning dependencies. ***
C
      DO 1200 IDP  = 1, NDDEP
      DO 1200 IZRO = 1, NIDEP
         IDEP (IDP,IZRO+1) = IRINDX(              IDXTMP(IZRO,2))
         COEFF(IDP,IZRO  ) =-RMTRIX(IDXTMP(IDP,1),IDXTMP(IZRO,2))
 1200 CONTINUE
C
C     *** Test print ***
C
      IF (IPRINT .GT. 20) THEN
         WRITE (LUPRI,'(A,I4)') 'Number of nonzero rows', NROW-NIDEP
         WRITE (LUPRI,'(A,16I4)') 'The non-zero rows are:',
     &                                        (IRINDX(II),II=1,NROW)
         DO II = 1, NDDEP
            WRITE (LUPRI,'(16F8.4)') (RMTRIX(II,IJ),IJ=1,NROW)
         END DO
         WRITE (LUPRI,'(A)') '                                   '
         WRITE (LUPRI,'(A)') '                                   '
C
         WRITE (LUPRI,'(2X,A)') 'Dependent component    ' //
     &                       'independent component     coeffisients'
         DO II = 1, NDDEP
            WRITE (LUPRI,'(10X,I4,20X,I3,15X,F8.4)') IDEP(II,1),
     &                                   IDEP(II,2), COEFF(II,1)
            DO IDP = 2, NIDEP
               WRITE (LUPRI,'(34X,I3,15X,F8.4)') IDEP(II,IDP+1),
     &                                          COEFF(II,IDP  )
            END DO
         END DO
         WRITE (LUPRI,'(A)') '                                   '
         WRITE (LUPRI,'(A)') '                                   '
      END IF
C
      RETURN
      END
C
C
C     /*Deck rmzror*/
      SUBROUTINE RMZROR(RMTRIX,IRINDX,IZINDX,NRDIM,NDIM,NROW,IPRINT)
C     ****************************************************************
C     *** Subroutine that removes the rows of a set of homogenous ***
C     *** linear equations, that does not contribute -> type:      ***
C     *** a_i*x_i = 0                                              ***
C     ****************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DMTHR = 1.0D-6)
      DIMENSION RMTRIX(NRDIM,NRDIM)
      DIMENSION IRINDX(NRDIM), IZINDX(NRDIM)
C
C     *** Test print before removal of rows ***
C
      IF (IPRINT .GT. 20) THEN
         WRITE (LUPRI,'(A)') 'Matrix before removal of rows'
         DO II = 1, NDIM
            WRITE (LUPRI,'(16F8.4)') (RMTRIX(II,IJ),IJ=1,NDIM)
         END DO
      END IF
C
C     *** Indexing for removal of rows                     ***
C     *** IRINDX -> which rows are left after the removal. ***
C     *** IZINDX -> which rows are supposed to be removed. ***
C
      IRROW = 0
      IZROW = 0
      MXDIM = NDIM
      DO 100 IDIM1 = 1, NDIM
         INOZRO = 0
         DO 200 IDIM2= 1, NDIM
            IF (ABS(RMTRIX(IDIM1,IDIM2)).GT.DMTHR) INOZRO = INOZRO + 1
 200     CONTINUE
C
         IF (INOZRO.GT.1) THEN
            IRROW = IRROW + 1
            IRINDX(IRROW) = IDIM1
         ELSE
            IZROW = IZROW + 1
            IZINDX(IZROW) = IDIM1
         END IF
 100  CONTINUE
      NROW  = IRROW
      NZROW = IZROW
C
C     *** removal of the rows. ***
C
      MXDIM = NDIM
      DO 300 IROW = 1, NZROW
         DO 400 IDIM2 = 1, NDIM
         DO 400 IDIM1 = IZINDX(IROW)-IROW+1, MXDIM
            RMTRIX(IDIM1,IDIM2) = RMTRIX(IDIM1+1,IDIM2)
 400     CONTINUE
         MXDIM = MXDIM - 1
 300  CONTINUE
C
C     *** removal of the corresponding coloumns ***
C
      MXDIM = NDIM
      DO 500 ICOL = 1, NZROW
         DO 600 IDIM2 = IZINDX(ICOL)-ICOL+1, MXDIM
         DO 600 IDIM1 = 1, NDIM
            RMTRIX(IDIM1,IDIM2) = RMTRIX(IDIM1,IDIM2+1)
 600     CONTINUE
         MXDIM = MXDIM - 1
 500  CONTINUE
C
C     *** Test print after removal of rows ***
C
      IF (IPRINT .GT. 20) THEN
         WRITE (LUPRI,'(A)') 'Matrix after removal of rows'
         DO II = 1, NROW
            WRITE (LUPRI,'(16F8.4)') (RMTRIX(II,IJ),IJ=1,NROW)
         END DO
      END IF
C
      RETURN
      END
C
C
C     /* Deck chkcmp */
      SUBROUTINE CHKCMP(KDPMTX,IDXTST,KORDR,LDPMTX,IFRSTD,NSTART,NLDPMX,
     &                  NIDEP,IPRINT,EXIST)
C     *************************************************************************
C     *** This subroutine checks whether the force constant IDXTST exist in ***
C     *** KDPMTX. It returns exist = .true. if this is the case.            ***
C     *************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <numder.h>
#include <fcsym.h>
      LOGICAL EXIST
      DIMENSION KDPMTX(LDPMTX,NSTRDR,IFRSTD), IDXTST(NMORDR)
C
      DO 100 IDP = 1, NIDEP+1
      DO 100 II   = NSTART, NLDPMX
         EXIST = .TRUE.
         DO 200 IORDR = 1, KORDR
            EXIST = EXIST.AND.(KDPMTX(II,IORDR,IDP).EQ.IDXTST(IORDR))
 200     CONTINUE
         IF (EXIST) GOTO 300
 100  CONTINUE
C
C     *** if it gets here with exist = .true. (Through the goto), ***
C     *** the component exist.                                    ***
C
 300  CONTINUE
C
      IF (IPRINT.GT.20) THEN
         WRITE (LUPRI,'(A,I4,A,I4)') 'Search startst at',NSTART,
     &                               ' and ends at', NLDPMX
         WRITE (LUPRI,'(A,4I4)') 'The component',
     &                           (IDXTST(II),II=1,KORDR)
         WRITE (LUPRI,'(L1)') EXIST
      END IF
C
      RETURN
      END
C
C
C     /* Deck chklcp*/
      SUBROUTINE CHKLCP(DCOEFF,KDPMTX,ICRIRP,IDXTMP,INMIDP,IMVTMP,
     &                  ICMPID,LDPMTX,IFRSTD,NDMCMP,NUMCMP,NSTART,NIDEP,
     &                  IPRINT)
C     **********************************************************************
C     *** This subroutine checks if the independent component chosen     ***
C     *** is the component with the least energy calculations. If this   ***
C     *** is not the case, the components are rearranged such that it is ***
C     **********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D1 = 1.0D0, DTHR = 1.0D-6)
#include <trkoor.h>
#include <numder.h>
#include <fcsym.h>
      LOGICAL INCMNT, FOUND
      DIMENSION DCOEFF(LDPMTX,IFRSTD)
      DIMENSION KDPMTX(LDPMTX,NSTRDR,IFRSTD), ICRIRP(NCOOR,2),
     &          IDXTMP(NDMCMP), ICMPID(NDMCMP), INMIDP(NIDEP),
     &          IMVTMP(NDMCMP)
C
C     *** Test print ***
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Coordinate dependency before change.',-1)
         WRITE (LUPRI,'(5X,A)') 'Dependent components      ' //
     &                          'Independent components     Coeffisient'
         DO II = NSTART, NSTART+NUMCMP-NIDEP
            IF (NDMCMP.EQ.3) THEN
               WRITE (LUPRI,FMT=1003) (KDPMTX(II,IX,1),IX=1,NDMCMP),
     &                   (KDPMTX(II,IX,2),IX=1,NDMCMP), DCOEFF(II,1)
               DO IDP = 2, NIDEP
                  WRITE (LUPRI,FMT=1006)
     &              (KDPMTX(II,IX,IDP+1),IX=1,NDMCMP), DCOEFF(II,IDP)
               END DO
            ELSE IF (NDMCMP.EQ.4) THEN
               write (lupri,'(a,I4)') 'testkomp', II
               WRITE (LUPRI,FMT=1004) (KDPMTX(II,IX,1),IX=1,NDMCMP),
     &                   (KDPMTX(II,IX,2),IX=1,NDMCMP), DCOEFF(II,1)
               DO IDP = 2, NIDEP
                  WRITE (LUPRI,FMT=1008)
     &              (KDPMTX(II,IX,IDP+1),IX=1,NDMCMP), DCOEFF(II,IDP)
               END DO
            END IF
         END DO
         WRITE (LUPRI,'(A)')'                                 '
         WRITE (LUPRI,'(A)')'                                 '
      END IF
C
C     *** Finding the number of alike variables in the ***
C     *** independent variables.                       ***
C
      DO 100 IDP = 1, NIDEP
C
         CALL IZERO(IDXTMP,NDMCMP)
C
         INUM = 1
         IDXTMP(1) = 1
         ICMPID(1) = KDPMTX(NSTART,1,IDP+1)
         DO 200 II = 2, NDMCMP
            INCMNT = .TRUE.
            DO 300 IJ = 1, INUM
               IF (KDPMTX(NSTART,II,IDP+1).EQ.ICMPID(IJ)) THEN
                  IDXTMP(IJ) = IDXTMP(IJ) + 1
                  INCMNT = .FALSE.
               END IF
 300        CONTINUE
C
            IF (INCMNT) THEN
               INUM = INUM + 1
               IDXTMP(INUM) = 1
               ICMPID(INUM) = KDPMTX(NSTART,II,IDP+1)
            END IF
 200     CONTINUE
         NUMB = INUM
C
         INMIDP(IDP) = 0
         DO 400 II = 1, NUMB
           INMIDP(IDP) = MAX(INMIDP(IDP),IDXTMP(II))
 400     CONTINUE
C
 100  CONTINUE
C
C     *** Comparing this with the dependent variables. ***
C
      DO 500 ICMP = NSTART, NSTART+NUMCMP-NIDEP
         CALL IZERO(IDXTMP,NDMCMP)
C
         INUM = 1
         IDXTMP(1) = 1
         ICMPID(1) = KDPMTX(ICMP,1,1)
         DO 600 II = 2, NDMCMP
            INCMNT = .TRUE.
            DO 700 IJ = 1, INUM
               IF (KDPMTX(ICMP,II,1).EQ.ICMPID(IJ)) THEN
                  IDXTMP(IJ) = IDXTMP(IJ) + 1
                  INCMNT = .FALSE.
               END IF
 700        CONTINUE
C
            IF (INCMNT) THEN
               INUM = INUM + 1
               IDXTMP(INUM) = 1
               ICMPID(INUM) = KDPMTX(ICMP,II,1)
            END IF
 600     CONTINUE
         NUMB = INUM
C
         MAXDP = 0
         DO 800 II = 1, NUMB
            MAXDP = MAX(MAXDP,IDXTMP(II))
 800     CONTINUE
C
         FOUND = .FALSE.
         DO 900 IDP = 1, NIDEP
            IF ((MAXDP.GT.INMIDP(IDP)).AND.(.NOT.FOUND).AND.
     &                           (ABS(DCOEFF(ICMP,IDP)).GT.1.0D-4)) THEN
               FOUND = .TRUE.
               IMIN = 1
               NMIN = INMIDP(1)
               DO 1100 IDP2 = 2, NIDEP
                  IF (INMIDP(IDP2).LT.NMIN) THEN
                     IMIN = IDP2
                     NMIN = INMIDP(IDP2)
                  END IF
 1100          CONTINUE
C
C              *** Swapping coordinates with respect to this criteria ***
C
               DO 1200 II = 1, NDMCMP
                  IMVTMP(II) = KDPMTX(ICMP,II,1)
                  KDPMTX(ICMP,II,1) = KDPMTX(NSTART,II,IMIN+1)
                  DO 1300 IST = NSTART, NSTART+NUMCMP-NIDEP
                     KDPMTX(IST,II,IMIN+1) = IMVTMP(II)
 1300             CONTINUE
 1200          CONTINUE
               INMIDP(IDP) = MAXDP
C
C              *** Arranging the coeffisients for the new and ***
C              *** better point. And placing the old as       ***
C              *** independent.                               ***
C
               BCOEFF = D1/DCOEFF(ICMP,IMIN)
               DO 1400 IDP2 = 1, NIDEP
                  IF (IDP2 .EQ. IMIN) THEN
                     DCOEFF(ICMP,IDP2) = BCOEFF
                  ELSE
                     DCOEFF(ICMP,IDP2) = -DCOEFF(ICMP,IDP2)*BCOEFF
                  END IF
 1400          CONTINUE
C
C              *** Changing the other coeffisients to fit the ***
C              *** new and better point                       ***
C
               DO 1500 ICMP2 = NSTART, NSTART+NUMCMP-NIDEP
                  IF (ICMP2.NE.ICMP) THEN
                     DO 1600 IDP2 = 1, NIDEP
                        IF (IDP2 .NE. IMIN) THEN
                           DCOEFF(ICMP2,IDP2) = DCOEFF(ICMP2,IDP2)
     &                                        + DCOEFF(ICMP2,IMIN)
     &                                         *DCOEFF(ICMP ,IDP2)
                        END IF
 1600                CONTINUE
                     DCOEFF(ICMP2,IMIN) = DCOEFF(ICMP2,IMIN)*BCOEFF
                  END IF
 1500          CONTINUE

               IF (IPRINT .GT. 20) THEN
                  CALL HEADER('Coordinate dependency after 1. change.',
     &                        -1)
                  WRITE (LUPRI,'(5X,A)') 'Dependent components      ' //
     &                 'Independent components     Coeffisient'
                  DO II = NSTART, NSTART+NUMCMP-NIDEP
                     IF (NDMCMP.EQ.3) THEN
                        WRITE (LUPRI,FMT=1003)
     &                       (KDPMTX(II,IX,1),IX=1,NDMCMP),
     &                       (KDPMTX(II,IX,2),IX=1,NDMCMP), DCOEFF(II,1)
                        DO IDP2 = 2, NIDEP
                           WRITE (LUPRI,FMT=1006)
     &                               (KDPMTX(II,IX,IDP2+1),IX=1,NDMCMP),
     &                               DCOEFF(II,IDP2)
                        END DO
                     ELSE IF (NDMCMP.EQ.4) THEN
                        WRITE (LUPRI,FMT=1004)
     &                       (KDPMTX(II,IX,1),IX=1,NDMCMP),
     &                       (KDPMTX(II,IX,2),IX=1,NDMCMP), DCOEFF(II,1)
                        DO IDP2 = 2, NIDEP
                           WRITE (LUPRI,FMT=1008)
     &                           (KDPMTX(II,IX,IDP2+1),IX=1,NDMCMP),
     &                            DCOEFF(II,IDP2)
                        END DO
                     END IF
                  END DO
                  WRITE (LUPRI,'(A)')'                                 '
                  WRITE (LUPRI,'(A)')'                                 '
               END IF
            END IF
 900     CONTINUE
C
C        *** Checking wether the component is better with respect to ***
C        *** which row it belongs to (row 1 of a 2 dimensional irrep ***
C        *** is preffered.                                           ***
C
         IF (.NOT.FOUND) THEN
C
C           *** Finds the worst independent variable available ***
C
            MXIDEP = 0
            MXIDX  = 0
            DO 1700 IDP = 1, NIDEP
               IIDEP = 0
               DO 1800 II = 1, NDMCMP
                  IIDEP = IIDEP + ICRIRP(KDPMTX(ICMP,II,IDP+1),2)
 1800          CONTINUE
               IF ((IIDEP.GT.MXIDEP) .AND.
     &                         (ABS(DCOEFF(ICMP,IDP)).GT.DTHR)) THEN
                  MXIDX = IDP
                  MXIDEP = IIDEP
               END IF
 1700       CONTINUE
C
            IDEP = 0
            DO 1900 II = 1, NDMCMP
               IDEP  = IDEP  + ICRIRP(KDPMTX(ICMP,II,1),2)
 1900       CONTINUE
C
C           *** If mxidep .gt. idep then the number coordinates that ***
C           *** transforms as row 2 i larger in iidep then in idep.  ***
C           *** We then change places.                               ***
C
            IF (MXIDEP .GT. IDEP) THEN
               FOUND = .TRUE.
C
               DO 2100 II = 1, NDMCMP
                  IMVTMP(II) = KDPMTX(ICMP,II,1)
                  KDPMTX(ICMP,II,1) = KDPMTX(NSTART,II,MXIDX+1)
                  DO 2200 IST = NSTART, NSTART+NUMCMP-NIDEP
                     KDPMTX(IST,II,MXIDX+1) = IMVTMP(II)
 2200             CONTINUE
 2100          CONTINUE
C
C              *** Arranging the coeffisients for the new and ***
C              *** better point. And placing the old as       ***
C              *** independent.                               ***
C
               BCOEFF = D1/DCOEFF(ICMP,MXIDX)
               DO 2300 IDP2 = 1, NIDEP
                  IF (IDP2 .EQ. MXIDX) THEN
                     DCOEFF(ICMP,IDP2) = BCOEFF
                  ELSE
                     DCOEFF(ICMP,IDP2) = -DCOEFF(ICMP,IDP2)*BCOEFF
                  END IF
 2300          CONTINUE
C
C              *** Changing the other coeffisients to fit the ***
C              *** new and better point                       ***
C
               DO 2400 ICMP2 = NSTART, NSTART+NUMCMP-NIDEP
                  IF (ICMP2.NE.ICMP) THEN
                     DO 2500 IDP2 = 1, NIDEP
                        IF (IDP2 .NE. MXIDX) THEN
                           DCOEFF(ICMP2,IDP2) = DCOEFF(ICMP2,IDP2 )
     &                                        + DCOEFF(ICMP2,MXIDX)
     &                                         *DCOEFF(ICMP ,IDP2 )
                        END IF
 2500                CONTINUE
                     DCOEFF(ICMP2,MXIDX)=DCOEFF(ICMP2,MXIDX)*BCOEFF
                  END IF
 2400          CONTINUE
            END IF
         END IF
         IF (IPRINT .GT. 20) THEN
            CALL HEADER('Coordinate dependency after 2. change. ',-1)
            WRITE (LUPRI,'(5X,A)') 'Dependent components      ' //
     &           'Independent components     Coeffisient'
            DO II = NSTART, NSTART+NUMCMP-NIDEP
               IF (NDMCMP.EQ.3) THEN
                  WRITE (LUPRI,FMT=1003) (KDPMTX(II,IX,1),IX=1,NDMCMP),
     &                       (KDPMTX(II,IX,2),IX=1,NDMCMP), DCOEFF(II,1)
                  DO IDP = 2, NIDEP
                     WRITE (LUPRI,FMT=1006)
     &                 (KDPMTX(II,IX,IDP+1),IX=1,NDMCMP), DCOEFF(II,IDP)
                  END DO
               ELSE IF (NDMCMP.EQ.4) THEN
                  WRITE (LUPRI,FMT=1004) (KDPMTX(II,IX,1),IX=1,NDMCMP),
     &                     (KDPMTX(II,IX,2),IX=1,NDMCMP), DCOEFF(II,1)
                  DO IDP = 2, NIDEP
                     WRITE (LUPRI,FMT=1008)
     &                 (KDPMTX(II,IX,IDP+1),IX=1,NDMCMP), DCOEFF(II,IDP)
                  END DO
               END IF
            END DO
            WRITE (LUPRI,'(A)')'                                 '
            WRITE (LUPRI,'(A)')'                                 '
         END IF
 500  CONTINUE
C
 1003 FORMAT (10X,3I4,10X,3I4,15X,F8.4)
 1006 FORMAT (32X,3I4,15X,F8.4)
 1004 FORMAT (6X,4I4,10X,4I4,11X,F8.4)
 1008 FORMAT (28X,4I4,15X,F8.4)
C
      RETURN
      END
C
C
C     /* Deck srtncf */
      SUBROUTINE SRTNCF(EQUMTX,ITMPRM,ICRWMX,ITMPEQ,NRDIM,KORDR,LDIM,
     &                  IPRINT)
C     ****************************************************
C     *** Sorting the components such that they become ***
C     *** calculated ones (X1 .GE. X2 .GE. X3 ....)    ***
C     ****************************************************
#include <implicit.h>
#include <priunit.h>
C
      LOGICAL FOUND, EXIST
      DIMENSION EQUMTX(NRDIM,NRDIM)
      DIMENSION ICRWMX(NRDIM,KORDR), ITMPRM(KORDR), ITMPEQ(NRDIM)
C
C     *** Test print before for comparison ***
C
      IF (IPRINT .GE. 20) THEN
         WRITE (LUPRI,'(5X,A)') 'Components before sorting'
C
         DO II =1 , LDIM
            WRITE (LUPRI,'(10X,20I4)') (ICRWMX(II,IJ),IJ=1,KORDR)
         END DO
      END IF
C
      DO 100 II = 1, LDIM
C
C        *** Sorting such that the indeces are aligned ***
C        ***         in descending order               ***
C
         IMAX = 0
         ITMPRM(1) = ICRWMX(II,1)
         DO 200 IJ = 2, KORDR
C
            FOUND = .FALSE.
            IMAX = IMAX + 1
            DO 300 IK = 1, IMAX
               IF ((ITMPRM(IK).LT.ICRWMX(II,IJ)).AND.(.NOT.FOUND)) THEN
                  FOUND = .TRUE.
                  DO 400 IL = IMAX, IK, -1
                     ITMPRM(IL+1) = ITMPRM(IL)
 400              CONTINUE
                  ITMPRM(IK) = ICRWMX(II,IJ)
               END IF
 300        CONTINUE
            IF (.NOT.FOUND) THEN
               ITMPRM(IMAX+1) = ICRWMX(II,IJ)
            END IF
 200     CONTINUE
C
C        *** Copying it back to the index array ***
C
         DO 500 IJ = 1, KORDR
            ICRWMX(II,IJ) = ITMPRM(IJ)
 500     CONTINUE
C
 100  CONTINUE
C
C     *** Removing redundant rows (with respect to equal force ***
C     *** constants) in the equation matrix                    ***
C
      IREM = 0
      LSTART = LDIM
      DO 600 II = 1, LSTART
C
C        *** Finding out if this force constant really already ***
C        *** exist, and which force constant this is.          ***
C
         EXIST = .FALSE.
         DO 700 IJ = 1, II-1
            IF (.NOT.EXIST) THEN
               EXIST = .TRUE.
               DO 800 IK = 1, KORDR
                  IF (ICRWMX(II,IK).NE.ICRWMX(IJ,IK)) EXIST = .FALSE.
 800           CONTINUE
               IF (EXIST) IEXIST = IJ
            END IF
 700     CONTINUE
C
C        *** If this component exist the matrix elements of  ***
C        ***       component will be added to the old.       ***
C
         IF (EXIST) THEN
            DO 900 IK = 1, LSTART
               EQUMTX(IK,IEXIST) = EQUMTX(IK,IEXIST) + EQUMTX(IK,II)
 900        CONTINUE
            IREM = IREM + 1
            ITMPEQ(IREM) = II
         END IF
 600  CONTINUE
C
C     *** Removing the appropriate rows and coloumns. ***
C
      DO 1100 IK = 1, IREM
         DO 1200 IJ = ITMPEQ(IK), LDIM-1
         DO 1200 II = 1         , LDIM
            EQUMTX(II,IJ) =  EQUMTX(II,IJ+1)
 1200    CONTINUE
C
         DO 1300 IJ = 1         , LDIM
         DO 1300 II = ITMPEQ(IK), LDIM-1
            EQUMTX(II,IJ) =  EQUMTX(II+1,IJ)
 1300    CONTINUE
C
         DO 1400 IJ = 1, KORDR
         DO 1400 II = ITMPEQ(IK), LDIM
            ICRWMX(II,IJ) = ICRWMX(II+1,IJ)
 1400    CONTINUE
C
         DO 1500 II = IK+1, IREM
            ITMPEQ(II) = ITMPEQ(II) - 1
 1500    CONTINUE
C
         LDIM = LDIM - 1
 1100 CONTINUE
C
C     *** Test print ***
C
      IF (IPRINT .GE. 20) THEN
C
         WRITE (LUPRI,'(5X,A,I4)') 'Size of matrix:', LDIM
         WRITE (LUPRI,'(5X,A)') 'Components after sorting'
C
         DO II =1 , LDIM
            WRITE (LUPRI,'(10X,20I4)') (ICRWMX(II,IJ),IJ=1,KORDR)
         END DO
C
         WRITE (LUPRI,'(A)') 'Matrix after sorting:'
         DO II = 1, LDIM
            WRITE (LUPRI,'(10X,20F8.4)') (EQUMTX(II,IJ),IJ=1,LDIM)
         END DO
      END IF
C
      RETURN
      END
C
C
C     /* Deck prtnrp */
      SUBROUTINE PRTNRP(GRIREP,WORK,ICRIRP,INDSTP,NPRTNR,MAXINR,IINNER,
     &                  IORDR,IRSRDR,NMPRTN,IDXTMP,INDTMP,LWORK,CLNRGY,
     &                  PRTNR,ALRCAL)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, DTHRSH = 1.0D-6, MXIDCR=5)
#include <fcsym.h>
#include <numder.h>
      LOGICAL PRTNR, CLNRGY, C2NRGY, TSM, ALRCAL
      DIMENSION GRIREP(NGORDR,NGVERT), WORK(LWORK)
      DIMENSION ICRIRP(NDCOOR,2), INDSTP(NTORDR), NPRTNR(MAXINR),
     &          INDTMP(MXIDCR), IDXTMP(NMORDR)
C
      PRTNR = .TRUE.
C
C     *** Checking wether this energy is already calculated ***
C     *** using another geometry as "partner-geometry.      ***
C
      ALRCAL = .FALSE.
      DO 100 INM = 1, NMPRTN
         IF ((NPRTNR(INM).EQ.IINNER).AND.(.NOT.ALRCAL)) THEN
            CLNRGY = .FALSE.
            PRTNR  = .FALSE.
            ALRCAL = .TRUE.
         END IF
 100  CONTINUE
C
C     *** Checking if all stepping coordinates are tot. sym. ***
C
      ITST = 1
      TSM  = .FALSE.
      DO 200 I = 1, IRSRDR+1
         ITST = ITST*ICRIRP(INDSTP(I),1)
 200  CONTINUE
      IF (ITST.EQ.1) TSM = .TRUE.
C
C     *** Checking if this geometry has a partner geometry ***
C
      IF ((IORDR.LE.2) .AND. (NMORDR.EQ.2) .AND.PRTNR
     &                                     .AND..NOT.TSM) THEN
C
C        *** Stepping in two directions. ***
C
         CALL SCNPRT(GRIREP,WORK,INDSTP,ICRIRP,INDTMP,IORDR,LWORK,PRTNR)
c      ELSE ((IORDR.LE.3) .AND. (PRTNR).AND.(.NOT.TSM)) THEN
c
      ELSE
         PRTNR = .FALSE.
      END IF
C
C     *** If there is a partner geometry, it's about ***
C     ***             time to find it.               ***
C
      IF (PRTNR) THEN
         NMPRTN = NMPRTN + 1
         ITINNR = IINNER - 1
         DO 400 I = IORDR, 1, -1
            IS = 2**(I-1)
            IDXTMP(I) = (ITINNR-MOD(ITINNR,IS))/IS
            ITINNR = ITINNR - IDXTMP(I)*IS
 400     CONTINUE
C
         DO 500 I = 1, IORDR
            IF (ICRIRP(INDSTP(I),1).NE.1) THEN
               NPRTNR(NMPRTN) = NPRTNR(NMPRTN)
     &                        + ABS(IDXTMP(I)-1)*2**(I-1)
            END IF
 500     CONTINUE
         NPRTNR(NMPRTN) = NPRTNR(NMPRTN) + 1
C
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE SCNPRT(GRIREP,WORK,INDSTP,ICRIRP,INDTMP,IORDR,LWORK,
     &                  PRTNR)
C     ********************************************************
C     *** Checking if this geometry has a partner geometry ***
C     *** Sufficient demand for stepping in 1-2 directions ***
C     *** is that all possible F_{x_1,x_2,x_3} in taylor   ***
C     *** expansion is zero.                               ***
C     ********************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (MXIDCR=5)
#include <fcsym.h>
#include <numder.h>
      LOGICAL PRTNR, C2NRGY
      DIMENSION GRIREP(NGORDR,NGVERT), WORK(LWORK)
      DIMENSION ICRIRP(NDCOOR,2), INDSTP(NTORDR), INDTMP(MXIDCR)
C
      DO 100 IIRDR3 = 1, IORDR
      DO 100 IIRDR2 = 1, IIRDR3
      DO 100 IIRDR1 = 1, IIRDR2
         IIREP3 = ICRIRP(INDSTP(IIRDR3),1)
         IIREP2 = ICRIRP(INDSTP(IIRDR2),1)
         IIREP1 = ICRIRP(INDSTP(IIRDR1),1)
C
         IF ((IIREP3.NE.1).OR.(IIREP2.NE.1).OR.(IIREP1.NE.1)) THEN
C
            INDTMP(1) = INDSTP(IIRDR1)
            INDTMP(2) = INDSTP(IIRDR2)
            INDTMP(3) = INDSTP(IIRDR3)
C
            KIRPST = 1
            KIRPDG = KIRPST + MXIDCR
            KKIRPD = KIRPDG + MXIDCR
            KIDDEG = KKIRPD + MXIDCR
            KIDEGI = KIDDEG + MXIDCR
            KLAST  = KIDEGI + MXIDCR
            IF (KLAST.GT.LWORK) CALL QUIT('Memory exceeded in' //
     &                                    ' PRTPNR.')
C
            C2NRGY = .FALSE.
            CALL SYMMLT(GRIREP,INDTMP,ICRIRP,WORK(KIRPST),WORK(KIRPDG),
     &                  WORK(KKIRPD),WORK(KIDDEG),WORK(KIDEGI),MXIDCR,
     &                  NTORDR,2,3,IPRINT,C2NRGY)
            PRTNR = PRTNR .AND. .NOT. C2NRGY
         END IF
 100  CONTINUE
C
      RETURN
      END
