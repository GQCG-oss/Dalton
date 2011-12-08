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
C  /* Deck sochk2 */
      SUBROUTINE SOCHK2(WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
C
      DIMENSION WORK(LWORK)
C
#include <sotabs.h>
C
#include <cotabs.h>
C
      CALL MKTAB
      JCOMID = ICOMID*(ICOMID + 1)/2
      JCOMAX = ICOMAX*(ICOMAX + 1)/2
      KCOMAX = JCOMAX*(JCOMAX + 1)/2
C
C     ***** Two-electron part *****
C
      KSOINT = 1
      K2INT  = KSOINT + 3*JCOMID*JCOMID
      KBUF   = K2INT  + KCOMAX
      KIBUF  = KBUF + 600
      KLAST  = KIBUF + 300
      IF (KLAST .GT. LWORK) CALL STOPIT('SOCHK2','TWO',KLAST,LWORK)
      CALL SCDRV2(WORK(KSOINT),JCOMID,WORK(K2INT),KCOMAX,
     *            WORK(KBUF),WORK(KIBUF),IPRINT)
C
C     ***** One-electron part *****
C
      KSOINT = 1
      K1INT  = KSOINT + 3*ICOMID*ICOMID
      KBUF   = K1INT  +   ICOMAX*ICOMAX
      KIBUF  = KBUF + 600
      KSCR   = KIBUF + 300
      KSCR2  = KSCR + JCOMID
      KLAST  = KSCR2 + JCOMAX
      IF (KLAST .GT. LWORK) CALL STOPIT('SOCHK2','ONE',KLAST,LWORK)
      CALL SCDRV1(WORK(KSOINT),ICOMID,WORK(K1INT),ICOMAX,
     *            WORK(KBUF),WORK(KIBUF),WORK(KSCR),JCOMID,
     *            WORK(KSCR2),JCOMAX,IPRINT)
      RETURN
      END
C  /* Deck scdrv1 */
      SUBROUTINE SCDRV1(SOINT,ISODIM,ONEINT,IONDIM,BUF,IBUF,
     *                  SCR,JCOMID,SCR2,JCOMAX,IPRALL)
C
C Test drive for spin-orbit integrals.
C
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <shells.h>
#include <nuclei.h>
#include <primit.h>
#include <inftap.h>
#include <sotabs.h>
#include <cotabs.h>
C
      DIMENSION SOINT(ISODIM,ISODIM,3), ONEINT(IONDIM,IONDIM),
     *          BUF(600), IBUF(600), SCR(JCOMID), SCR2(JCOMAX)
C
      REWIND LUPROP
      CALL MOLLAB('X1SPNORB',LUPROP,LUPRI)
      READ (LUPROP) SCR
      CALL DAPTGE(ICOMID,SCR,SOINT(1,1,1))
      CALL HEADER('Spin-orbit integrals (x component)',-1)
      CALL OUTPUT(SOINT(1,1,1),1,ICOMID,1,ICOMID,ICOMID,ICOMID,1,LUPRI)
C
      REWIND LUPROP
      CALL MOLLAB('Y1SPNORB',LUPROP,LUPRI)
      READ (LUPROP) SCR
      CALL DAPTGE(ICOMID,SCR,SOINT(1,1,2))
      CALL HEADER('Spin-orbit integrals (y component)',-1)
      CALL OUTPUT(SOINT(1,1,2),1,ICOMID,1,ICOMID,ICOMID,ICOMID,1,LUPRI)
C
      REWIND LUPROP
      CALL MOLLAB('Z1SPNORB',LUPROP,LUPRI)
      READ (LUPROP) SCR
      CALL DAPTGE(ICOMID,SCR,SOINT(1,1,3))
      CALL HEADER('Spin-orbit integrals (z component)',-1)
      CALL OUTPUT(SOINT(1,1,3),1,ICOMID,1,ICOMID,ICOMID,ICOMID,1,LUPRI)
C
      CALL DZERO(SCR2,JCOMAX)
      CALL RDINT1(SCR2,BUF,IBUF)
      CALL DSPTSI(ICOMAX,SCR2,ONEINT)
      CALL HEADER('Nuclear attraction integrals',-1)
      CALL OUTPUT(ONEINT,1,ICOMAX,1,ICOMAX,ICOMAX,ICOMAX,1,LUPRI)
C
C Compare the equivalent integrals
C
      CALL CMP1I(SOINT,ICOMID,ONEINT,ICOMAX)
      RETURN
      END
C  /* Deck rdint1 */
      SUBROUTINE RDINT1(SCR,BUF,IBUF)
#include <implicit.h>
#include <priunit.h>
#include <inftap.h>
C
      DIMENSION SCR(*)
C
      DIMENSION BUF(600)
      DIMENSION IBUF(600)
C
      REWIND LUONEL
      CALL MOLLAB('ONEHAMIL',LUONEL,LUERR)
  100 READ (LUONEL) BUF, IBUF, LENGTH
      IF(LENGTH .EQ. 0) GO TO 100
      IF(LENGTH .LT. 0) GO TO 110
         DO 120 I = 1, LENGTH
            SCR(IBUF(I)) = BUF(I)
  120    CONTINUE
         GO TO 100
  110 CONTINUE
      REWIND LUONEL
      CALL MOLLAB('KINETINT',LUONEL,LUERR)
  200 READ (LUONEL) BUF, IBUF, LENGTH
      IF(LENGTH .EQ. 0) GO TO 200
      IF(LENGTH .LT. 0) GO TO 210
         DO 220 I = 1, LENGTH
            SCR(IBUF(I)) = SCR(IBUF(I)) - BUF(I)
  220    CONTINUE
         GO TO 200
  210 CONTINUE
      RETURN
      END
C  /* Deck rdint2 */
      SUBROUTINE RDINT2(SCR,BUF,IBUF)
C
C     Reads in one-electron hamilton integrals. kr, aug. 1992
C
#include <implicit.h>
#include <priunit.h>
#include <inftap.h>
C
      DIMENSION SCR(*), BUF(600), IBUF(600)
C
      REWIND LUONEL
      CALL MOLLAB('ONEHAMIL',LUONEL,LUERR)
 100  READ (LUONEL) BUF,IBUF,LENGTH
      IF (LENGTH .EQ. 0) GOTO 100
      IF (LENGTH .LT. 0) GOTO 110
      DO 120 I = 1, LENGTH
         SCR(IBUF(I)) = BUF(I)
 120  CONTINUE
      GOTO 100
 110  CONTINUE
      RETURN
      END
C  /* Deck cmp1i */
      SUBROUTINE CMP1I(SOINT,ISODIM,ONEINT,IONDIM)
C
C Compares one-electron spin-orbital integrals with the
C equivalent ones constructed explicitely from nuclear
C attraction integrals.
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxmom.h>
C
      DIMENSION SOINT(ISODIM,ISODIM,3)
      DIMENSION ONEINT(IONDIM,IONDIM)
C
#include <sotabs.h>
#include <cotabs.h>
#include <shells.h>
#include <primit.h>
#include <xyzpow.h>
C
      DIMENSION JR(2),KR(2),JS(2),KS(2)
      DIMENSION JRSHL(2),KSSHL(2),KRSHL(2),JSSHL(2)
      DIMENSION JRCMP(2),KSCMP(2),KRCMP(2),JSCMP(2)
      LOGICAL FIRST,LAST,DIFFER
      DIMENSION RENR1(2), RENR2(2), RENS1(2), RENS2(2)
C
      DIMENSION ICOOR1(3),ICOOR2(3)
      CHARACTER COOR(3)
      DATA ICOOR1/2,3,1/,ICOOR2/3,1,2/
      DATA COOR/'x','y','z'/
      DATA PRTHRS/1D-13/
C
#include <ibtfun.h>

      SDPRE(I) = (1.0D0*I)

C
      DIFMAX   = 0D0
      RENR1(1) = 1D0
      RENR2(1) = 1D0
      RENS1(1) = 1D0
      RENS2(1) = 1D0
      FIRST = .TRUE.
      LAST = .FALSE.
      DIFFER = .FALSE.
C
C For each possible spin-orbit integral construct the equivalent
C from the electron repulsion integrals.
C
      DO 320  IR = 1, ISODIM
         DO 330  IS = 1, IR
            DO 340 ICOOR = 1, 3
C
C The shell and component of IR and IS:
C
               IRSHL = ISHLCO(IR)
               IRCMP = ICMPCO(IR)
               ISSHL = ISHLCO(IS)
               ISCMP = ICMPCO(IS)
C
C The ICOOR component of a cross product, i.e. I,J,K always
C forms an even permutation of 1,2,3
C
C
               JCOOR = ICOOR1(ICOOR)
               KCOOR = ICOOR2(ICOOR)
C
               NRANG = NHKT(IRSHL) - 1
               NSANG = NHKT(ISSHL) - 1
C
               IF (JCOOR .EQ. 1) THEN
                  RENR1(2) = SDPRE(NRANG - ISTEP(IRCMP))
                  RENS2(2) = SDPRE(NSANG - ISTEP(ISCMP))
               ELSE IF (JCOOR .EQ. 2) THEN
                  RENR1(2) = SDPRE(MVAL(IRCMP))
                  RENS2(2) = SDPRE(MVAL(ISCMP))
               ELSE
                  RENR1(2) = SDPRE(NVAL(IRCMP))
                  RENS2(2) = SDPRE(NVAL(ISCMP))
               END IF
C
               IF (KCOOR .EQ. 1) THEN
                  RENS1(2) = SDPRE(NSANG - ISTEP(ISCMP))
                  RENR2(2) = SDPRE(NRANG - ISTEP(IRCMP))
               ELSE IF (KCOOR .EQ. 2) THEN
                  RENS1(2) = SDPRE(MVAL(ISCMP))
                  RENR2(2) = SDPRE(MVAL(IRCMP))
               ELSE
                  RENS1(2) = SDPRE(NVAL(ISCMP))
                  RENR2(2) = SDPRE(NVAL(IRCMP))
               END IF
C
C Find the shells and components of the orbitals that result from
C differentiation.
C
               DO 10,II=1,2
                 JRSHL(II) = ISHLDF(JCOOR,II,IRCMP,IRSHL)
                 JRCMP(II) = ICMPDF(JCOOR,II,IRCMP,IRSHL)
C
                 KRSHL(II) = ISHLDF(KCOOR,II,IRCMP,IRSHL)
                 KRCMP(II) = ICMPDF(KCOOR,II,IRCMP,IRSHL)
C
                 JSSHL(II) = ISHLDF(JCOOR,II,ISCMP,ISSHL)
                 JSCMP(II) = ICMPDF(JCOOR,II,ISCMP,ISSHL)
C
                 KSSHL(II) = ISHLDF(KCOOR,II,ISCMP,ISSHL)
                 KSCMP(II) = ICMPDF(KCOOR,II,ISCMP,ISSHL)
C
C Find the orbital label for each of shell-component pair
C
                 JR(II) = ICSHCM(JRSHL(II),JRCMP(II))
                 KR(II) = ICSHCM(KRSHL(II),KRCMP(II))
                 JS(II) = ICSHCM(JSSHL(II),JSCMP(II))
                 KS(II) = ICSHCM(KSSHL(II),KSCMP(II))
10            CONTINUE
C
              SQTEXR = SQRT(PRIEXP(IRSHL))
              SQTEXS = SQRT(PRIEXP(ISSHL))
              FACEXP = SQTEXR*SQTEXS
              VAL = 0D0
              DO 20 NR = 1,2
                 SGNR = 1D0
                 IF (NR .EQ. 2) SGNR = - 1D0
                 DO 21 NS = 1,2
                    SGNRS = SGNR
                    IF (NS .EQ. 2) SGNRS = - SGNR
                    FACTOR = SDPRE(NR*NS)*FACEXP*SGNRS
                    NX = JR(NR)
                    NY = KS(NS)
                    VAL = VAL + RENR1(NR)*RENS1(NS)*FACTOR*ONEINT(NX,NY)
                    NY = KR(NR)
                    NX = JS(NS)
                    VAL = VAL - RENR2(NR)*RENS2(NS)*FACTOR*ONEINT(NX,NY)
21               CONTINUE
20            CONTINUE
C
C Print the integrals if they differ
C
              SOVAL = SOINT(IR,IS,ICOOR)
              ABSDIF = ABS(SOVAL - VAL)
              DIFMAX = MAX(DIFMAX,ABSDIF)
              IF (ABSDIF.GT.PRTHRS) THEN
                 IF (FIRST) THEN
                    WRITE(LUPRI,*)
                    WRITE(LUPRI,'(A,D8.1,A)')
     *  'Differences found: (threshold,',PRTHRS,')'
                    WRITE(LUPRI,'(A)')'-----------------'
                    WRITE(LUPRI,'(17X,A)')
     *'   Spin-orbit integral     Test integral       Difference'
                    WRITE(LUPRI,*)
                    DIFFER = .TRUE.
                    FIRST  = .FALSE.
                 ENDIF
                 WRITE(LUPRI,'(1X,A,2I3,5X,2D23.15,D13.3)')
     *             COOR(ICOOR),IR,IS,SOVAL,VAL,ABSDIF
              ENDIF
340        CONTINUE
330     CONTINUE
320   CONTINUE
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) ' Maximum difference ', DIFMAX
      RETURN
      END
C  /* Deck scdrv2 */
      SUBROUTINE SCDRV2(SOINT,JCOMID,ERINT,KCOMAX,BUF,IBUF,IPRALL)
C
C Test drive for spin-orbit integrals.
C
C
#include <implicit.h>
#include <dummy.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <mxcent.h>
C
#include <inftap.h>
#include <shells.h>
#include <nuclei.h>
#include <primit.h>
#include <sotabs.h>
#include <cotabs.h>
C
      DIMENSION SOINT(JCOMID,JCOMID,3), ERINT(KCOMAX),
     *          BUF(600), IBUF(600)
C
      CALL GPOPEN(LU1,'AO2SOINT','OLD',' ',' ',IDUMMY,.FALSE.)
C
C     Clear integral arrays
C
      CALL DZERO(SOINT,3*JCOMID*JCOMID)
      CALL DZERO(ERINT,KCOMAX)
C
C Read files to buffers.
C
      CALL RDINT(LU1,'AO2SOINT',SOINT,BUF,IBUF,JCOMID*JCOMID,JCOMID)
      CALL RDINT(LUINTA,'BASTWOEL',ERINT,BUF,IBUF,KCOMAX,ICOMAX)
C
C Compare the equivalent integrals
C
      CALL CMP2I(SOINT,ICOMID,JCOMID,ERINT,ICOMAX,KCOMAX)
      CALL GPCLOSE(LU1,'KEEP')
      RETURN
      END
C  /* Deck mktab */
      SUBROUTINE MKTAB
C
C Sets up table referring differentiated orbitals to
C equivalent orbitals of lower or higher order.
C Ghost orbitals are  assumed at end of MOLECULE input corresponding
C to, in order, the equivalent orbitals.
C
C Sets up another table relating all orbital indices to a shell
C index and a component index and vice versa
C
#include <implicit.h>
C
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
C
#include <aovec.h>
#include <shells.h>
#include <nuclei.h>
#include <sotabs.h>
#include <cotabs.h>
C
C Following implies max f:s, ie d derivatives.
C
      DIMENSION IDIF(3,2,10)
      DIMENSION JOFF(3)
      DATA IDIF/
C s
     *         1,2,3,0,0,0,
C p
     *         1,2,3,1,0,0,
     *         2,4,5,0,1,0,
     *         3,5,6,0,0,1,
C d
     *         1,2,3,1,0,0,
     *         2,4,5,2,1,0,
     *         3,5,6,3,0,1,
     *         4,7,8,0,2,0,
     *         5,8,9,0,3,2,
     *         6,9,10,0,0,3/
      DATA JOFF/0,1,4/
C
C Find index of last original shell
C
      DO 10 K = 1,KMAX
         Z = CHARGE(NCENT(K))
         IF (Z.EQ.0D0) GOTO 100
10    CONTINUE
100   CONTINUE
      KMAX1 = K - 1
C
C Let the derivative of each shell component point to two
C other shells and components - ISHLDF(2),ICMPDF(2)
C
      KF = KMAX1
      DO 20 K = 1,KMAX1
         KJ  = NHKT(K)
         KM  = KHKT(K)
C
C s shells only contibute to the step up derivaties
C All others contribute to both step up and stop down
C
         DO 21 KDU=2,1,-1
            IF(KDU.EQ.2.AND.KJ.EQ.1) GOTO 21
            KF = KF + 1
            DO 22 KK=1,KM
               DO 23 KI = 1,3
                  ISHLDF(KI,KDU,KK,K) = KF
                  ICMPDF(KI,KDU,KK,K) = IDIF(KI,KDU,JOFF(KJ)+KK)
23             CONTINUE
22          CONTINUE
21       CONTINUE
C
20    CONTINUE
      IF (KF.NE.KMAX) THEN
         CALL QUIT('WRONG MOLECULE INPUT FOR SPIN-ORBIT CHECK')
      END IF
C
C Orbital - shell,  component tables.
C
      ICONTR = 0
      DO 30, K = 1, KMAX
         KJ = NHKT(K)
         KM = KHKT(K)
         DO 31, KK = 1, KM
            ICONTR = ICONTR + 1
            ISHLCO(ICONTR) = K
            ICMPCO(ICONTR) = KK
            ICSHCM(K,KK) = ICONTR
31       CONTINUE
         IF (K.EQ.KMAX1) ICOMID = ICONTR
30    CONTINUE
      ICOMAX = ICONTR
      RETURN
      END
C  /* Deck rdint */
      SUBROUTINE RDINT(LUINP,KEY,DINT,BUF,IBUF,IDIM,IBAS)
#include <implicit.h>
#include <priunit.h>
C
C Reads two-electron integrals from file to buffer.
C LUINP: File to read
C LABEL: What section of the file to be read
C BGBUF: Points to start element in big integral buffer.
C IBGBUF: Points to start element in big label buffer.
C IDIM: Self-explanatory
C ILAST : Offset of last element returned to calling routine
C
      LOGICAL SPNORB
      CHARACTER*8 KEY
      DIMENSION DINT(IDIM,3)
C
      DIMENSION BUF(600)
      DIMENSION IBUF(600)
C
#include <inftap.h>
#include <ibtfun.h>
      REWIND LUINP
      CALL MOLLAB(KEY,LUINP,LUPRI)
      SPNORB = KEY .EQ. 'AO2SOINT'
 150  READ (LUINP, END = 300 ) BUF, IBUF, LENGTH
         IF (LENGTH .GT. 0) THEN
            DO 100 I = 1, LENGTH
               LABEL = IBUF(I)
               IF (SPNORB) THEN
                  IS     = IBTAND(LABEL,255)
                  IF (IS .EQ. 0) THEN
                     ICOOR = IBTAND(IBTSHR(LABEL,8),255)
                  ELSE
                     IP = IBTAND(IBTSHR(LABEL,24),255)
                     IQ = IBTAND(IBTSHR(LABEL,16),255)
                     IR = IBTAND(IBTSHR(LABEL, 8),255)
                     IS = IBTAND(LABEL,255)
                     IPQ = IP*(IP - 1)/2 + IQ
                     IRS = IR*(IR - 1)/2 + IS
                     IPQRS = (IPQ - 1)*IBAS + IRS
                     DINT(IPQRS,ICOOR) = BUF(I)
                  ENDIF
               ELSE
                  IP = IBTAND(IBTSHR(LABEL,24),255)
                  IQ = IBTAND(IBTSHR(LABEL,16),255)
                  IR = IBTAND(IBTSHR(LABEL, 8),255)
                  IS = IBTAND(LABEL,255)
                  IPQ = IP*(IP - 1)/2 + IQ
                  IRS = IR*(IR - 1)/2 + IS
                  IPQRS = IPQ*(IPQ - 1)/2 + IRS
                  DINT(IPQRS,1) = BUF(I)
               END IF
 100        CONTINUE
         ELSE IF (LENGTH .LT. 0 ) THEN
            GO TO 300
         END IF
         GO TO 150
C
 300  CONTINUE
      RETURN
      END
C  /* Deck cmp2i */
      SUBROUTINE CMP2I(SOINT,ISODIM,JCOMID,ERINT,IERDIM,KCOMAX)
C
C Compares two-electron spin-orbital integrals with the
C equivalent ones constructed explicitely from electron
C repulsion integrals.
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxmom.h>
C
      DIMENSION SOINT(JCOMID,JCOMID,3)
      DIMENSION ERINT(KCOMAX)
C
#include <sotabs.h>
#include <cotabs.h>
#include <shells.h>
#include <primit.h>
#include <xyzpow.h>
C
      CHARACTER COOR(3)
      DIMENSION ICOOR1(3),ICOOR2(3)
      DIMENSION JR(2),KR(2),JS(2),KS(2)
      DIMENSION JRSHL(2),KSSHL(2),KRSHL(2),JSSHL(2)
      DIMENSION JRCMP(2),KSCMP(2),KRCMP(2),JSCMP(2)
C
      LOGICAL FIRST,LAST,DIFFER
      DIMENSION RENR1(2), RENR2(2), RENS1(2), RENS2(2)
C
      DATA PRTHRS/1D-13/
      DATA COOR/'x','y','z'/
      DATA ICOOR1/2,3,1/,ICOOR2/3,1,2/
C
#include <ibtfun.h>

      SDPRE(I) = (1.0D0*I)

C
      DIFMAX   = 0D0
      RENR1(1) = 1D0
      RENR2(1) = 1D0
      RENS1(1) = 1D0
      RENS2(1) = 1D0
      FIRST = .TRUE.
      LAST = .FALSE.
      DIFFER = .FALSE.
C
C For each possible spin-orbit integral construct the equivalent
C from the electron repulsion integrals.
C
      DO 300  IP = 1, ISODIM
         DO 310  IQ = 1, IP
            DO 320  IR = 1, ISODIM
               DO 330  IS = 1, IR
                  DO 340 ICOOR = 1, 3
C
C The shell and component of IR and IS:
C
                     IRSHL = ISHLCO(IR)
                     IRCMP = ICMPCO(IR)
                     ISSHL = ISHLCO(IS)
                     ISCMP = ICMPCO(IS)
C
C The ICOOR component of a cross product, i.e. I,J,K always
C forms an even permutation of 1,2,3
C
C
                     JCOOR = ICOOR1(ICOOR)
                     KCOOR = ICOOR2(ICOOR)
C
                     NRANG = NHKT(IRSHL) - 1
                     NSANG = NHKT(ISSHL) - 1
C
                     IF (JCOOR .EQ. 1) THEN
                        RENR1(2) = SDPRE(NRANG - ISTEP(IRCMP))
                        RENS2(2) = SDPRE(NSANG - ISTEP(ISCMP))
                     ELSE IF (JCOOR .EQ. 2) THEN
                        RENR1(2) = SDPRE(MVAL(IRCMP))
                        RENS2(2) = SDPRE(MVAL(ISCMP))
                     ELSE
                        RENR1(2) = SDPRE(NVAL(IRCMP))
                        RENS2(2) = SDPRE(NVAL(ISCMP))
                     END IF
C
                     IF (KCOOR .EQ. 1) THEN
                        RENS1(2) = SDPRE(NSANG - ISTEP(ISCMP))
                        RENR2(2) = SDPRE(NRANG - ISTEP(IRCMP))
                     ELSE IF (KCOOR .EQ. 2) THEN
                        RENS1(2) = SDPRE(MVAL(ISCMP))
                        RENR2(2) = SDPRE(MVAL(IRCMP))
                     ELSE
                        RENS1(2) = SDPRE(NVAL(ISCMP))
                        RENR2(2) = SDPRE(NVAL(IRCMP))
                     END IF
C
C Find the shells and components of the orbitals that result from
C differentiation.
C
                     DO 10,II=1,2

                       JRSHL(II) = ISHLDF(JCOOR,II,IRCMP,IRSHL)
                       JRCMP(II) = ICMPDF(JCOOR,II,IRCMP,IRSHL)
C

                       KRSHL(II) = ISHLDF(KCOOR,II,IRCMP,IRSHL)
                       KRCMP(II) = ICMPDF(KCOOR,II,IRCMP,IRSHL)
C
                       JSSHL(II) = ISHLDF(JCOOR,II,ISCMP,ISSHL)
                       JSCMP(II) = ICMPDF(JCOOR,II,ISCMP,ISSHL)
C
                       KSSHL(II) = ISHLDF(KCOOR,II,ISCMP,ISSHL)
                       KSCMP(II) = ICMPDF(KCOOR,II,ISCMP,ISSHL)
C
C Find the orbital label for each of shell-component pair
C
                       JR(II) = ICSHCM(JRSHL(II),JRCMP(II))
                       KR(II) = ICSHCM(KRSHL(II),KRCMP(II))
                       JS(II) = ICSHCM(JSSHL(II),JSCMP(II))
                       KS(II) = ICSHCM(KSSHL(II),KSCMP(II))
10                  CONTINUE
C
C For all eight combinations check the canonical order.
C Search the ER buffer for the one unless we have a zero label..
C A zero label corresponds to a differentiation yielding zero.
C
C
C NOTE This renormalization only holds for non-contracted
C basis sets
C
                    SQTEXR = SQRT(PRIEXP(IRSHL))
                    SQTEXS = SQRT(PRIEXP(ISSHL))
                    FACEXP = SQTEXR*SQTEXS
                    VAL = 0D0
                    DO 20 NR = 1,2
                       SGNR = 1D0
                       IF (NR .EQ. 2) SGNR = - 1D0
                       DO 21 NS = 1,2
                          SGNRS = SGNR
                          IF (NS .EQ. 2) SGNRS = - SGNR
                          FACTOR = SDPRE(NR*NS)*FACEXP*SGNRS
C
C First term
C
                          NA = IP
                          NB = IQ
                          NX = JR(NR)
                          NY = KS(NS)
                          CALL CANON(NA,NB,NX,NY)
                          IF (NB.NE.0.AND.NY.NE.0) THEN
                             NAB = NA*(NA - 1)/2 + NB
                             NXY = NX*(NX - 1)/2 + NY
                             NABXY = NAB*(NAB - 1)/2 + NXY
                             VAL = VAL +
     *                          RENR1(NR)*RENS1(NS)*FACTOR*ERINT(NABXY)
                          ENDIF
C
C Second term
C
                          NA = IP
                          NB = IQ
                          NY = KR(NR)
                          NX = JS(NS)
C
C NOTE This renormalization only holds for non-contracted
C basis sets
C
                          CALL CANON(NA,NB,NY,NX)
                          IF (NB.NE.0.AND.NX.NE.0) THEN
                             NAB = NA*(NA - 1)/2 + NB
                             NXY = NY*(NY - 1)/2 + NX
                             NABXY = NAB*(NAB - 1)/2 + NXY
                             VAL = VAL -
     *                          RENR2(NR)*RENS2(NS)*FACTOR*ERINT(NABXY)
                          ENDIF
21                     CONTINUE
20                  CONTINUE
C
C Print the integrals if they differ
C
                    IPQ = IP*(IP - 1)/2 + IQ
                    IRS = IR*(IR - 1)/2 + IS
                    SOVAL = SOINT(IRS,IPQ,ICOOR)
                    ABSDIF = ABS(SOVAL - VAL)
                    DIFMAX = MAX(DIFMAX,ABSDIF)
                    IF (ABSDIF.GT.PRTHRS) THEN
                       IF (FIRST) THEN
                          WRITE(LUPRI,*)
                          WRITE(LUPRI,'(A,D8.1,A)')
     *  'Differences found: (threshold,',PRTHRS,')'
                          WRITE(LUPRI,'(A)')'-----------------'
                          WRITE(LUPRI,'(17X,A)')
     *'   Spin-orbit integral     Test integral       Difference'
                          WRITE(LUPRI,*)
                          DIFFER = .TRUE.
                          FIRST  = .FALSE.
                       ENDIF
                       WRITE(LUPRI,'(1X,A,4I3,5X,2D23.15,D13.3)')
     *                   COOR(ICOOR),IP,IQ,IR,IS,SOVAL,VAL,ABSDIF
                    ENDIF
340              CONTINUE
330           CONTINUE
320        CONTINUE
310      CONTINUE
300   CONTINUE
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) ' Maximum difference ', DIFMAX
      IF (.NOT.DIFFER) THEN
         WRITE(LUPRI,'(A)') 'No differences found'
         WRITE(LUPRI,*)
      ENDIF
      RETURN
      END
C  /* Deck mg1tst */
      SUBROUTINE MG1TST(WORK,LWORK,IPRINT,LABEL,LABELT,ORIGIN,WHATIN)
#include <implicit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <inforb.h>
C
      DIMENSION WORK(LWORK), ORIGIN(3)
      CHARACTER*8 LABEL, LABELT(3)
      CHARACTER*4 WHATIN
C
      KIQM   = 1
      KJCO   = KIQM   + (NUCIND + 1)/IRAT
      KAUXVE = KJCO   + (4*NUCIND + 1)/IRAT
      KTSTVE = KAUXVE + NBASIS*NBASIS
      KRESVE = KTSTVE + NBASIS*NBASIS
      KQMN   = KRESVE + NBASIS*NBASIS
      KSTVEC = KQMN + 9
      KRSVEC = KSTVEC + 3
      KMXDIF = KRSVEC + 3
      KNONTV = KMXDIF + 3
      KDIFVC = KNONTV + (NUCIND + 1)/IRAT
      KLAST  = KDIFVC + NBASIS*NBASIS
      IF (KLAST .GT. LWORK) CALL STOPIT('MG1TST',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL MG1TS1(WORK(KIQM),WORK(KJCO),WORK(KAUXVE),
     &            WORK(KTSTVE),WORK(KRESVE),WORK(KQMN),LABEL,LABELT,
     &            WORK(KSTVEC),WORK(KRSVEC),WORK(KDIFVC),
     &            WORK(KMXDIF),WORK(KNONTV),ORIGIN,WHATIN,WORK(KLAST),
     &            LWRK,IPRINT)
      RETURN
      END
C  /* Deck mg1ts1 */
      SUBROUTINE MG1TS1(IQM,JCO,AUXVEC,TSTVEC,RESVEC,QMN,LABEL,
     &                  LABELT,STVEC,RSVEC,DIFVEC,MAXDIF,NONTVC,
     &                  ORIGIN,WHATIN,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <primit.h>
#include <nuclei.h>
      PARAMETER (D0 = 0.D0)
      PARAMETER (PRTHRS = 1D-13, D2 = 2.D0, DP5 = 0.5D0)
      LOGICAL FNDLAB, DIFFER
      CHARACTER*8 LABEL, LABELT(3)
      CHARACTER*4 WHATIN
C
      DIMENSION WORK(LWORK), IQM(NUCIND), JCO(NUCIND,4),
     &          AUXVEC(NBASIS*NBASIS),
     &          TSTVEC(NBASIS*NBASIS/3,3),
     &          RESVEC(NBASIS*NBASIS/3,3),
     &          QMN(3,3),STVEC(3),RSVEC(3),IY(6,3),RMXDIF(3),
     &          NONTVC(NUCIND), DIFVEC(NBASIS*NBASIS/3,3),
     &          ORIGIN(3)
C
#include <inftap.h>
#include <shells.h>
      DATA ((IY(I,J), J = 1,3), I = 1,6) /0,1,2,0,2,3,0,2,3,0,3,4,
     &                                    0,3,4,0,3,4/
C
      CALL DZERO(TSTVEC,NBASIS*NBASIS)
      CALL TSTDAT(NONTYP,NONTVC,IQM,JCO,ICOUNT,NCOUNT,NDIFNU)
C
      IF (WHATIN .EQ. 'S1ML' .OR. WHATIN .EQ. 'S1MR' .OR.
     &    WHATIN .EQ. 'MQDP') THEN
         NCOUN2 = NCOUNT*ICOUNT
      ELSE
         NCOUN2 = NCOUNT*(NCOUNT + 1)/2
      END IF
      ICOUN2 = ICOUNT*(ICOUNT + 1)/2
C
C     Read integrals
C
      DO 10 I = 1, 3
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABELT(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABELT(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READT(LUPROP,NCOUN2,AUXVEC)
         IF (WHATIN .EQ. 'S1ML' .OR. WHATIN .EQ. 'S1MR' .OR.
     &       WHATIN .EQ. 'MQDP') THEN
            DO 11 I1 = 1, NCOUNT
               DO 11 I2 = 1, NCOUNT
                  IADR1 = (I1 - 1)*NCOUNT + I2
                  IADR2 = (I1 - 1)*ICOUNT + I2
                  RESVEC(IADR1,I) = AUXVEC(IADR2)
 11         CONTINUE
         ELSE
            CALL DCOPY(NCOUN2,AUXVEC,1,RESVEC(1,I),1)
         END IF
 10   CONTINUE
      IF (LABEL .EQ. 'ONEHAMIL') THEN
         KBUF = 1
         KIBUF = KBUF + 600
         KLAST = KIBUF + 600
         IF (KLAST .GT. LWORK) CALL STOPIT('MG1TS1',' ',KLAST,LWORK)
         LWRK = LWORK - KLAST + 1
         CALL DZERO(AUXVEC,ICOUN2*IRAT)
         CALL RDINT2(AUXVEC,WORK(KBUF),WORK(KIBUF))
      ELSE
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABEL,LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &           'Label ', LABEL, ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READT(LUPROP,ICOUN2,AUXVEC)
      END IF
C
C     Calculation of the new integrals
C
      DIFFER = .FALSE.
      QMN(1,1) = D0
      QMN(2,2) = D0
      QMN(3,3) = D0
      CALL DZERO(RMXDIF,3)
      ICOUNR = 0
      NUCNM1 = 0
      JPRIMA = 0
      DO 30 I = 1, NDIFNU
         DO 30 II = 1, NONTVC(I)
         NUCNM1 = NUCNM1 + 1
         DO 30 J = 1, IQM(I)
         DO 30 K = 1, JCO(I,J)
         JPRIMA = JPRIMA + 1
         DO 30 L = 1, J*(J + 1)/2
            ICOUNR = ICOUNR + 1
            ICOUNL = 0
            NUCNM2 = 0
            DO 40 I1 = 1, NDIFNU
               DO 40 II1 = 1, NONTVC(I1)
               NUCNM2 = NUCNM2 +1
               DO 40 J1 = 1, IQM(I1)
               DO 40 K1 = 1, JCO(I1,J1)
               DO 40 L1 = 1, J1*(J1 + 1)/2
                  ICOUNL = ICOUNL + 1
                  IF (WHATIN .EQ. 'S1ML') THEN
                     IADR = (ICOUNL - 1)*NCOUNT + ICOUNR
                     QMN(1,2) = ORIGIN(3) - CORD(3,NUCNM1)
                     QMN(1,3) = CORD(2,NUCNM1) - ORIGIN(2)
                     QMN(2,3) = ORIGIN(1) - CORD(1,NUCNM1)
                  ELSE IF (WHATIN .EQ. 'S1MR') THEN
                     IADR = (ICOUNL - 1)*NCOUNT + ICOUNR
                     QMN(1,2) = CORD(3,NUCNM2) - ORIGIN(3)
                     QMN(1,3) = ORIGIN(2) - CORD(2,NUCNM2)
                     QMN(2,3) = CORD(1,NUCNM2) - ORIGIN(1)
                  ELSE IF (WHATIN .EQ. 'MQDP') THEN
                     IADR = (ICOUNL - 1)*NCOUNT + ICOUNR
                  ELSE IF (ICOUNR .GE. ICOUNL) THEN
                     IADR = ICOUNR*(ICOUNR - 1)/2 + ICOUNL
                     QMN(1,2) = CORD(3,NUCNM2) - CORD(3,NUCNM1)
                     QMN(1,3) = CORD(2,NUCNM1) - CORD(2,NUCNM2)
                     QMN(2,3) = CORD(1,NUCNM2) - CORD(1,NUCNM1)
                  ELSE
                     GOTO 101
                  END IF
                  QMN(2,1) = - QMN(1,2)
                  QMN(3,1) = - QMN(1,3)
                  QMN(3,2) = - QMN(2,3)
                  IADR2 = NCOUNT*(NCOUNT + 1)/2
                  NROW2 = NCOUNT
                  NNUC = 0
                  DO 100 I2 = NDIFNU + 1, NONTYP
                     DO 100 II2 = 1, NONTVC(I2)
                     NNUC = NNUC + 1
                     DO 100 J2 = 1, IQM(I2)
                        DO 100 K2 = 1, JCO(I2,J2)
                           DO  100 L2 = 1, J2*(J2 + 1)/2
                              NROW2 = NROW2 + 1
                              SIGN = 1.0D0
                              IF (ICOUNR .GE. ICOUNL) THEN
                                 IADR5 = ICOUNR*(ICOUNR - 1)/2 + ICOUNL
                              ELSE
                                 IF (WHATIN .EQ. 'MQDP') THEN 
                                    IADR5 = ICOUNL*(ICOUNL - 1)/2+ICOUNR
                                    SIGN = -1.0D0
                                 ELSE
                                    IADR5 = ICOUNL*(ICOUNL - 1) + ICOUNR
                                 END IF
                              END IF
                              IF ((NUCNM1 .EQ. NNUC) .AND.
     &                             (J2 .EQ. (J + 1)) .AND. (K .EQ. K2)
     &                              .AND. (L .EQ. L2)) THEN
                                 IADR3 = IADR2
                                 IADR4 = IADR2
                                 DO 110 NN = 1, IY(L,2)
                                    IADR3 = IADR3 + NROW2 + NN - 1
 110                             CONTINUE
                                 DO 120 NN = 1, IY(L,3)
                                    IADR4 = IADR4 + NROW2 + NN - 1
 120                             CONTINUE
                                 IADRX = IADR2 + ICOUNL
                                 IADRY = IADR3 + ICOUNL
                                 IADRZ = IADR4 + ICOUNL
                              ELSE
                                 IADR2 = IADR2 + NROW2
                              END IF
 100                       CONTINUE
                  STVEC(1) = AUXVEC(IADRX)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + SIGN*CORD(1,NUCNM1)*AUXVEC(IADR5)
                  STVEC(2) = AUXVEC(IADRY)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + SIGN*CORD(2,NUCNM1)*AUXVEC(IADR5)
                  STVEC(3) = AUXVEC(IADRZ)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + SIGN*CORD(3,NUCNM1)*AUXVEC(IADR5)
                  IF (WHATIN .EQ. 'MQDP') THEN
                     CALL DCOPY(3,STVEC,1,RSVEC,1)
                  ELSE
                     CALL DGEMM('N','N',3,1,3,1.D0,
     &                          QMN,3,
     &                          STVEC,3,0.D0,
     &                          RSVEC,3)
                  END IF
                  DO 50 I2 = 1, 3
                     IF (WHATIN .EQ. 'MQDP') THEN
                        TSTVEC(IADR,I2) = 1.0D0/3.0D0*RSVEC(I2)
                     ELSE
                        TSTVEC(IADR,I2) = DP5*RSVEC(I2)
                     END IF
                     DIFRS = RESVEC(IADR,I2)
     &                     - TSTVEC(IADR,I2)
                     DIFFER = DIFFER .OR. (ABS(DIFRS) .GT. PRTHRS)
                     RMXDIF(I2) = MAX(RMXDIF(I2),ABS(DIFRS))
                     DIFVEC(IADR,I2) = DIFRS
 50               CONTINUE
 101           CONTINUE
 40         CONTINUE
 30   CONTINUE
      IF (.NOT. DIFFER) THEN
         WRITE(LUPRI,'(/A/)')
     &         'No differences found in magnetic properties test'
      ELSE
         WRITE(LUPRI,'(/A/)')
     &         'Differences found in magnetic properties test!'
         WRITE(LUPRI,'(/A/)') 'Output of difference matrices follows'
         DO 300 I = 1, 3
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUN2,DIFVEC(1,I),1,AUXVEC,1)
            IF (WHATIN .EQ. 'S1ML' .OR. WHATIN .EQ. 'S1MR' .OR.
     &          WHATIN .EQ. 'MQDP') THEN
               CALL OUTPUT(AUXVEC,1,NCOUNT,1,NCOUNT,NCOUNT,NCOUNT,1,
     &                     LUPRI)
            ELSE
               CALL OUTPAK(AUXVEC,NCOUNT,1,LUPRI)
            END IF
 300     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &         'Output of test program integrals is as follows:'
         DO 320 I = 1, 3
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUN2,TSTVEC(1,I),1,AUXVEC,1)
            IF (WHATIN .EQ. 'S1ML' .OR. WHATIN .EQ. 'S1MR' .OR.
     &          WHATIN .EQ. 'MQDP') THEN
               CALL OUTPUT(AUXVEC,1,NCOUNT,1,NCOUNT,NCOUNT,NCOUNT,1,
     &                     LUPRI)
            ELSE
               CALL OUTPAK(AUXVEC,NCOUNT,1,LUPRI)
            END IF
 320     CONTINUE
      ENDIF
      RETURN
      END
C  /* Deck mgtst2 */
      SUBROUTINE MGTST2(WORK,LWORK,IPRINT,LABEL,LABELT)
#include <implicit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <inforb.h>
C
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL, LABELT(3)
C
      KIQM   = 1
      KJCO   = KIQM   + (NUCIND + 1)/IRAT
      KAUXVE = KJCO   + (MXQN*NUCIND + 1)/IRAT
      KTSTVE = KAUXVE + NBASIS*(NBASIS + 1)/2
      KRESVE = KTSTVE + 3*NBASIS*(NBASIS + 1)
      KQMN   = KRESVE + 3*NBASIS*(NBASIS + 1)
      KSTVEC = KQMN + 9
      KRSVEC = KSTVEC + 9
      KMXDIF = KRSVEC + 9
      KNONTV = KMXDIF + 6
      KDIFVC = KNONTV + (NUCIND + 1)/IRAT
      KLAST  = KDIFVC + 3*NBASIS*(NBASIS + 1)
      IF (KLAST .GT. LWORK) CALL STOPIT('MG1TST',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL MG1TS2(WORK(KIQM),WORK(KJCO),WORK(KAUXVE),
     &            WORK(KTSTVE),WORK(KRESVE),WORK(KQMN),LABEL,LABELT,
     &            WORK(KSTVEC),WORK(KRSVEC),WORK(KDIFVC),
     &            WORK(KMXDIF),WORK(KNONTV),WORK(KLAST),LWRK,IPRINT)
      RETURN
      END
C  /* Deck mg1ts2 */
      SUBROUTINE MG1TS2(IQM,JCO,AUXVEC,TSTVEC,RESVEC,QMN,LABEL,
     &                  LABELT,STVEC,RSVEC,DIFVEC,MAXDIF,NONTVC,
     &                  WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <primit.h>
#include <nuclei.h>
      PARAMETER (D0 = 0.D0, DP25 = 0.25D0)
      PARAMETER (PRTHRS = 1D-13, D2 = 2.D0, D4 = 4.D0)
      CHARACTER*8 LABEL, LABELT(3)
      LOGICAL FNDLAB, DIFFER
C
      DIMENSION WORK(LWORK), IQM(NUCIND), JCO(NUCIND,MXQN),
     &          AUXVEC(NBASIS*(NBASIS+1)/2),
     &          TSTVEC(NBASIS*(NBASIS + 1)/2,6),IADR(6),
     &          RESVEC(NBASIS*(NBASIS + 1)/2,6),IZ(6,6),
     &          QMN(3,3),STVEC(3,3),RSVEC(3,3),IY(6,3),RMXDIF(6),
     &          NONTVC(NUCIND), DIFVEC(NBASIS*(NBASIS + 1)/2,6)
C
#include <inftap.h>
#include <shells.h>
      DATA ((IY(I,J), J = 1,3), I = 1,6) /0,1,2,0,2,3,0,2,3,0,3,4,
     &                                    0,3,4,0,3,4/
      DATA ((IZ(I,J), J = 1,6), I = 1,6) /0,1,2,3,4,5,0,2,3,5,6,7,
     &                                    0,2,3,5,6,7,0,3,4,7,8,9,
     &                                    0,3,4,7,8,9,0,3,4,7,8,9/
C
      CALL TSTDAT(NONTYP,NONTVC,IQM,JCO,ICOUNT,NCOUNT,NDIFNU)
C
      ICOUN2 = ICOUNT*(ICOUNT + 1)/2
      NCOUN2 = NCOUNT*(NCOUNT + 1)/2
C
C     Read integrals
C
      DO 10 I = 1, 6
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABELT(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABELT(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,NCOUN2*IRAT,AUXVEC)
         CALL DCOPY(NCOUN2,AUXVEC,1,RESVEC(1,I),1)
 10   CONTINUE
      REWIND LUPROP
      IF (.NOT. FNDLAB(LABEL,LUPROP)) THEN
         WRITE (LUPRI,'(1X,3A)')
     &         'Label ', LABEL, ' not found on file AOPROPER'
         CALL QUIT('Error in magnetic properties test')
      END IF
      CALL READI(LUPROP,ICOUN2*IRAT,AUXVEC)
C
C     Calculation of the new integrals
C
      DIFFER = .FALSE.
      QMN(1,1) = D0
      QMN(2,2) = D0
      QMN(3,3) = D0
      CALL DZERO(RMXDIF,6)
      ICOUNR = 0
      NUCNM1 = 0
      JPRIMA = 0
      DO 30 I = 1, NDIFNU
         DO 30 II = 1, NONTVC(I)
         NUCNM1 = NUCNM1 + 1
         DO 30 J = 1, IQM(I)
         DO 30 K = 1, JCO(I,J)
         JPRIMA = JPRIMA + 1
         DO 30 L = 1, J*(J + 1)/2
            ICOUNR = ICOUNR + 1
            ICOUNL = 0
            NUCNM2 = 0
            DO 40 I1 = 1, NDIFNU
               DO 40 II1 = 1, NONTVC(I1)
               NUCNM2 = NUCNM2 +1
               DO 40 J1 = 1, IQM(I1)
               DO 40 K1 = 1, JCO(I1,J1)
               DO 40 L1 = 1, J1*(J1 + 1)/2
                  ICOUNL = ICOUNL + 1
                  IF (ICOUNR .GE. ICOUNL) THEN
                     IADR0 = ICOUNR*(ICOUNR - 1)/2 + ICOUNL
                  QMN(1,2) = CORD(3,NUCNM2) - CORD(3,NUCNM1)
                  QMN(2,1) = - QMN(1,2)
                  QMN(1,3) = CORD(2,NUCNM1) - CORD(2,NUCNM2)
                  QMN(3,1) = - QMN(1,3)
                  QMN(2,3) = CORD(1,NUCNM2) - CORD(1,NUCNM1)
                  QMN(3,2) = - QMN(2,3)
                  IADR2 = NCOUN2
                  NROW2 = NCOUNT
                  NNUC = 0
                  DO 100 I2 = NDIFNU + 1, NONTYP
                     DO 100 II2 = 1, NONTVC(I2)
                     NNUC = NNUC + 1
                     DO 100 J2 = 1, IQM(I2)
                        DO 100 K2 = 1, JCO(I2,J2)
                           DO  100 L2 = 1, J2*(J2 + 1)/2
                              NROW2 = NROW2 + 1
                              IF ((NUCNM1 .EQ. NNUC) .AND. (J2
     &                             .EQ. (J + 1)) .AND. (K .EQ. K2)
     &                             .AND. (L .EQ. L2)) THEN
                                 IADR3 = IADR2
                                 IADR4 = IADR2
                                 DO 130 NN = 1, IY(L,2)
                                    IADR3 = IADR3 + NROW2 + NN - 1
 130                             CONTINUE
                                 DO 140 NN = 1, IY(L,3)
                                    IADR4 = IADR4 + NROW2 + NN - 1
 140                             CONTINUE
                                 IADRX = IADR2 + ICOUNL
                                 IADRY = IADR3 + ICOUNL
                                 IADRZ = IADR4 + ICOUNL
                              ELSE IF ((NUCNM1 .EQ. NNUC) .AND. (J2
     &                             .EQ. (J + 2)) .AND. (K .EQ. K2)
     &                             .AND. (L .EQ. L2)) THEN
                                 DO 120 NN = 2, 6
                                 IADR(NN) = IADR2
                                 DO 110 NN2 = 1, IZ(L,NN)
                                    IADR(NN) = IADR(NN) + NROW2 + NN2-1
 110                             CONTINUE
 120                             CONTINUE
                                 IADRXX = IADR2 + ICOUNL
                                 IADRXY = IADR(2) + ICOUNL
                                 IADRXZ = IADR(3) + ICOUNL
                                 IADRYY = IADR(4) + ICOUNL
                                 IADRYZ = IADR(5) + ICOUNL
                                 IADRZZ = IADR(6) + ICOUNL
                              END IF
                              IADR2 = IADR2 + NROW2
 100                       CONTINUE
                  STVEC(1,1) = AUXVEC(IADRXX)/(D4*PRIEXP(JPRIMA))
     &                       + D2*CORD(1,NUCNM1)*AUXVEC(IADRX)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(1,NUCNM1)*CORD(1,NUCNM1)*AUXVEC(IADR0)
                  STVEC(2,2) = AUXVEC(IADRYY)/(D4*PRIEXP(JPRIMA))
     &                       + D2*CORD(2,NUCNM1)*AUXVEC(IADRY)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(2,NUCNM1)*CORD(2,NUCNM1)*AUXVEC(IADR0)
                  STVEC(3,3) = AUXVEC(IADRZZ)/(D4*PRIEXP(JPRIMA))
     &                       + D2*CORD(3,NUCNM1)*AUXVEC(IADRZ)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(3,NUCNM1)*CORD(3,NUCNM1)*AUXVEC(IADR0)
                  STVEC(1,2) = AUXVEC(IADRXY)/(D4*PRIEXP(JPRIMA))
     &                       + (AUXVEC(IADRX)*CORD(2,NUCNM1)
     &                       +  AUXVEC(IADRY)*CORD(1,NUCNM1))/
     &                          (D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(1,NUCNM1)*CORD(2,NUCNM1)*AUXVEC(IADR0)
                  STVEC(1,3) = AUXVEC(IADRXZ)/(D4*PRIEXP(JPRIMA))
     &                       + (AUXVEC(IADRX)*CORD(3,NUCNM1)
     &                       +  AUXVEC(IADRZ)*CORD(1,NUCNM1))/
     &                          (D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(1,NUCNM1)*CORD(3,NUCNM1)*AUXVEC(IADR0)
                  STVEC(2,3) = AUXVEC(IADRYZ)/(D4*PRIEXP(JPRIMA))
     &                       + (AUXVEC(IADRY)*CORD(3,NUCNM1)
     &                       +  AUXVEC(IADRZ)*CORD(2,NUCNM1))/
     &                          (D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(2,NUCNM1)*CORD(3,NUCNM1)*AUXVEC(IADR0)
                  STVEC(2,1) = STVEC(1,2)
                  STVEC(3,2) = STVEC(2,3)
                  STVEC(3,1) = STVEC(1,3)
                  CALL DGEMM('N','N',3,3,3,1.D0,
     &                       QMN,3,
     &                       STVEC,3,0.D0,
     &                       RSVEC,3)
                  CALL DGEMM('N','N',3,3,3,1.D0,
     &                       RSVEC,3,
     &                       QMN,3,0.D0,
     &                       STVEC,3)
                  INTNUM = 0
                  DO 50 I2 = 1, 3
                     DO 50 I3 = I2, 3
                        INTNUM = INTNUM + 1
                        TSTVEC(IADR0,INTNUM) = DP25*STVEC(I2,I3)
                        DIFRS = RESVEC(IADR0,INTNUM)
     &                        - TSTVEC(IADR0,INTNUM)
                        DIFFER = DIFFER .OR. (ABS(DIFRS) .GT. PRTHRS)
                        RMXDIF(INTNUM) =
     &                         MAX(RMXDIF(INTNUM),ABS(DIFRS))
                        DIFVEC(IADR0,INTNUM) = DIFRS
 50               CONTINUE
            END IF
 40         CONTINUE
 30   CONTINUE
      IF (.NOT. DIFFER) THEN
         WRITE(LUPRI,'(/A/)')
     &         'No differences found in magnetic properties test'
      ELSE
         WRITE(LUPRI,'(/A/)')
     &         'Differences found in magnetic properties test!'
         WRITE(LUPRI,'(/A/)') 'Output of difference matrices follows'
         DO 300 I = 1, 6
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUN2,DIFVEC(1,I),1,AUXVEC,1)
            CALL OUTPAK(AUXVEC,NCOUNT,1,LUPRI)
 300     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &         'Output of test program integrals is as follows:'
         DO 320 I = 1, 6
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUN2,TSTVEC(1,I),1,AUXVEC,1)
            CALL OUTPAK(AUXVEC,NCOUNT,1,LUPRI)
 320     CONTINUE
      ENDIF
      RETURN
      END
C  /* Deck hdbtst */
      SUBROUTINE HDBTST(WORK,LWORK,IPRINT,LABELD,NATOM,ORIGIN)
#include <implicit.h>
#include <iratdef.h>
#include <dummy.h>
#include <maxmom.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <inforb.h>
C
      DIMENSION WORK(LWORK), INTREP(9*MXCENT), INTRE2(9*MXCENT),
     &          ORIGIN(3)
      CHARACTER*8 LABEL(9*MXCENT), LABELD(3)
C
      KATOM = 1
      CALL SETATM(WORK(KATOM),NATOM,10)
      CALL HDBTYP(NOPTYP,INTRE2,IDUMMY,LABEL,WORK(KATOM),NATOM)
      KIQM   = KATOM + NUCDEP
      KJCO   = KIQM  + (NUCIND + 1)/IRAT
      KAUXVE = KJCO + (4*NUCIND + 1)/IRAT
      KTSTVE = KAUXVE + NBASIS*NBASIS
      KRESVE = KTSTVE + 9*NBASIS*NBASIS
      KDIPVL = KRESVE + 9*NBASIS*NBASIS
      KNONTV = KDIPVL + 3*NBASIS*(NBASIS + 1)/2
      KDIFVC = KNONTV + (NUCIND + 1)/IRAT
      KLAST  = KDIFVC + 9*NBASIS*NBASIS
      IF (KLAST .GT. LWORK) CALL STOPIT('NSTST2',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL HDBTS1(WORK(KIQM),WORK(KJCO),WORK(KAUXVE),
     &            WORK(KTSTVE),WORK(KRESVE),LABEL,LABELD,WORK(KDIFVC),
     &            WORK(KNONTV),WORK(KLAST),LWRK,IPRINT,WORK(KDIPVL),
     &            NOPTYP,ORIGIN)
      RETURN
      END
C  /* Deck hdbts1 */
      SUBROUTINE HDBTS1(IQM,JCO,AUXVEC,TSTVEC,RESVEC,LABEL,LABELD,
     &                  DIFVEC,NONTVC,WORK,LWORK,IPRINT,DIPVEL,NOPTYP,
     &                  ORIGIN)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <primit.h>
#include <nuclei.h>
      PARAMETER (D0 = 0.D0)
      PARAMETER (PRTHRS = 1D-13, D2 = 2.D0, DP5 = 0.5D0)
      LOGICAL FNDLAB, DIFFER
      CHARACTER*8 LABEL(NOPTYP), LABELD(3)
C
      DIMENSION WORK(LWORK), IQM(NUCIND), JCO(NUCIND,4),
     &          AUXVEC(NBASIS*NBASIS),
     &          TSTVEC(NBASIS*NBASIS/3,9),
     &          RESVEC(NBASIS*NBASIS/3,9),
     &          DIPVEL(NBASIS*(NBASIS + 1)/2,3),
     &          QMN(3,3),STVEC(3),RSVEC(3),IY(6,3),
     &          NONTVC(NUCIND), DIFVEC(NBASIS*NBASIS/3,9),
     &          ORIGIN(3)
C
#include <inftap.h>
#include <shells.h>
      DATA ((IY(I,J), J = 1,3), I = 1,6) /0,1,2,0,2,3,0,2,3,0,3,4,
     &                                    0,3,4,0,3,4/
C
      CALL TSTDAT(NONTYP,NONTVC,IQM,JCO,ICOUNT,NCOUNT,NDIFNU)
C
      NCOUN2 = NCOUNT*ICOUNT
      ICOUN2 = ICOUNT*(ICOUNT + 1)/2
C
C     Read integrals
C
      CALL DZERO(RESVEC,NBASIS*NBASIS*3)
      DO 10 I = 1, NOPTYP
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABEL(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABEL(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,NCOUN2*IRAT,AUXVEC)
         ICOMP = MOD(I,9)
         IF (ICOMP .EQ. 0) ICOMP = 9
         DO 11 I1 = 1, NCOUNT
            DO 11 I2 = 1, NCOUNT
               IADR1 = (I1 - 1)*NCOUNT + I2
               IADR2 = (I1 - 1)*ICOUNT + I2
               RESVEC(IADR1,ICOMP) = RESVEC(IADR1,ICOMP) + AUXVEC(IADR2)
 11      CONTINUE
 10   CONTINUE
      DO 15 I = 1, 3
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABELD(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &           'Label ', LABELD(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,ICOUN2*IRAT,AUXVEC)
         CALL DCOPY(ICOUN2,AUXVEC,1,DIPVEL(1,I),1)
 15   CONTINUE
C
C     Calculation of the new integrals
C
      DIFFER = .FALSE.
      QMN(1,1) = D0
      QMN(2,2) = D0
      QMN(3,3) = D0
      ICOUNR = 0
      NUCNM1 = 0
      DO 30 I = 1, NDIFNU
         DO 30 II = 1, NONTVC(I)
         NUCNM1 = NUCNM1 + 1
         DO 30 J = 1, IQM(I)
         DO 30 K = 1, JCO(I,J)
         DO 30 L = 1, J*(J + 1)/2
            ICOUNR = ICOUNR + 1
            ICOUNL = 0
            NUCNM2 = 0
            JPRIMB = 0
            DO 40 I1 = 1, NDIFNU
               DO 40 II1 = 1, NONTVC(I1)
               NUCNM2 = NUCNM2 +1
               DO 40 J1 = 1, IQM(I1)
               DO 40 K1 = 1, JCO(I1,J1)
                  JPRIMB = JPRIMB + 1
               DO 40 L1 = 1, J1*(J1 + 1)/2
                  ICOUNL = ICOUNL + 1
                  IADR = (ICOUNL - 1)*NCOUNT + ICOUNR
                  QMN(1,2) = CORD(3,NUCNM2) - ORIGIN(3)
                  QMN(1,3) = ORIGIN(2) - CORD(2,NUCNM2)
                  QMN(2,3) = CORD(1,NUCNM2) - ORIGIN(1)
                  QMN(2,1) = - QMN(1,2)
                  QMN(3,1) = - QMN(1,3)
                  QMN(3,2) = - QMN(2,3)
                  IF (ICOUNR .GE. ICOUNL) THEN
                     IADR5 = ICOUNR*(ICOUNR - 1)/2 + ICOUNL
                  ELSE
                     IADR5 = ICOUNL*(ICOUNL - 1)/2 + ICOUNR
                  END IF
                  NNUC = 0
                  NROW = NCOUNT
                  DO 100 I2 = NDIFNU + 1, NONTYP
                     DO 100 II2 = 1, NONTVC(I2)
                     NNUC = NNUC + 1
                     DO 100 J2 = 1, IQM(I2)
                        DO 100 K2 = 1, JCO(I2,J2)
                           DO  100 L2 = 1, J2*(J2 + 1)/2
                              NROW = NROW + 1
                              IF ((NUCNM2 .EQ. NNUC) .AND.
     &                             (J2 .EQ. (J1 + 1)) .AND. (K1 .EQ. K2)
     &                              .AND. (L1 .EQ. L2)) THEN
                                 IADR2 = (NROW - 1)*NROW/2
                                 IADR3 = IADR2
                                 IADR4 = IADR2
                                 DO 110 NN = 1, IY(L2,2)
                                    IADR3 = IADR3 + NROW + NN - 1
 110                             CONTINUE
                                 DO 120 NN = 1, IY(L2,3)
                                    IADR4 = IADR4 + NROW + NN - 1
 120                             CONTINUE
                                 IADRX = IADR2 + ICOUNR
                                 IADRY = IADR3 + ICOUNR
                                 IADRZ = IADR4 + ICOUNR
                              END IF
 100              CONTINUE
               DO 130 NN = 1, 3
                     STVEC(1) = DIPVEL(IADRX,NN)/
     &                    (D2*SQRT(PRIEXP(JPRIMB)))
     &                     + CORD(1,NUCNM2)*DIPVEL(IADR5,NN)
                  STVEC(2) = DIPVEL(IADRY,NN)/
     &                       (D2*SQRT(PRIEXP(JPRIMB)))
     &                     + CORD(2,NUCNM2)*DIPVEL(IADR5,NN)
                  STVEC(3) = DIPVEL(IADRZ,NN)/
     &                       (D2*SQRT(PRIEXP(JPRIMB)))
     &                     + CORD(3,NUCNM2)*DIPVEL(IADR5,NN)
                  CALL DGEMM('N','N',3,1,3,1.D0,
     &                       QMN,3,
     &                       STVEC,3,0.D0,
     &                       RSVEC,3)
                  DO 50 I2 = 1, 3
                     TSTVEC(IADR,3*(NN - 1) + I2) = DP5*RSVEC(I2)
                     DIFRS = RESVEC(IADR,3*(NN - 1) + I2)
     &                     - TSTVEC(IADR,3*(NN - 1) + I2)
                     DIFFER = DIFFER .OR. (ABS(DIFRS) .GT. PRTHRS)
                     DIFVEC(IADR,3*(NN - 1) + I2) = DIFRS
 50               CONTINUE
 130           CONTINUE
 40         CONTINUE
 30   CONTINUE
      IF (.NOT. DIFFER) THEN
         WRITE(LUPRI,'(/A/)')
     &         'No differences found in magnetic properties test'
      ELSE
         WRITE(LUPRI,'(/A/)')
     &         'Differences found in magnetic properties test!'
         WRITE(LUPRI,'(/A/)') 'Output of difference matrices follows'
         DO 300 I = 1, 9
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUNT*NCOUNT,DIFVEC(1,I),1,AUXVEC,1)
            CALL OUTPUT(AUXVEC,1,NCOUNT,1,NCOUNT,NCOUNT,NCOUNT,1,LUPRI)
 300     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &         'Output of test program integrals is as follows:'
         DO 320 I = 1, 9
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUNT*NCOUNT,TSTVEC(1,I),1,AUXVEC,1)
            CALL OUTPUT(AUXVEC,1,NCOUNT,1,NCOUNT,NCOUNT,NCOUNT,1,LUPRI)
 320     CONTINUE
      ENDIF
      RETURN
      END
C  /* Deck sustst */
      SUBROUTINE SUSTST(WORK,LWORK,IPRINT,LABEL,LABELT)
#include <implicit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <inforb.h>
C
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL(3), LABELT(6)
C
      KIQM   = 1
      KJCO   = KIQM   + (NUCIND + 1)/IRAT
      KAUXVE = KJCO   + (4*NUCIND + 1)/IRAT
      KTSTVE = KAUXVE + NBASIS*(NBASIS + 1)/2
      KRESVE = KTSTVE + 3*NBASIS*(NBASIS + 1)/2
      KSINTV = KRESVE + 6*NBASIS*(NBASIS + 1)/2
      KQMN   = KSINTV + 3*NBASIS*(NBASIS + 1)/2
      KSTVEC = KQMN + 9
      KRSVEC = KSTVEC + 9
      KMXDIF = KRSVEC + 9
      KNONTV = KMXDIF + 6
      KDIFVC = KNONTV + (NUCIND + 1)/IRAT
      KLAST  = KDIFVC + 3*NBASIS*(NBASIS + 1)/2
      IF (KLAST .GT. LWORK) CALL STOPIT('MG1TST',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL SUSTS1(WORK(KIQM),WORK(KJCO),WORK(KAUXVE),
     &            WORK(KTSTVE),WORK(KRESVE),WORK(KSINTV),WORK(KQMN),
     &            LABEL,LABELT,WORK(KSTVEC),WORK(KRSVEC),WORK(KDIFVC),
     &            WORK(KMXDIF),WORK(KNONTV),WORK(KLAST),LWRK,IPRINT)
      RETURN
      END
C  /* Deck susts1 */
      SUBROUTINE SUSTS1(IQM,JCO,AUXVEC,TSTVEC,RESVEC,SINTVC,QMN,
     &                  LABEL,LABELT,STVEC,RSVEC,DIFVEC,MAXDIF,NONTVC,
     &                  WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <primit.h>
#include <nuclei.h>
      PARAMETER (D0 = 0.D0, DP25 = 0.25D0)
      PARAMETER (PRTHRS = 1D-13, D2 = 2.D0)
      CHARACTER*8 LABEL(3), LABELT(6)
      LOGICAL FNDLAB, DIFFER
C
      DIMENSION WORK(LWORK), IQM(NUCIND), JCO(NUCIND,4),
     &          AUXVEC(NBASIS*(NBASIS+1)/2),
     &          TSTVEC(NBASIS*(NBASIS + 1)/6,3),
     &          RESVEC(NBASIS*(NBASIS + 1)/6,6),
     &          SINTVC(NBASIS*(NBASIS + 1)/2,3),
     &          QMN(3,3),STVEC(3,3),RSVEC(3,3),IY(6,3),RMXDIF(3),
     &          NONTVC(NUCIND), DIFVEC(NBASIS*(NBASIS + 1)/2,3)
C
#include <inftap.h>
#include <shells.h>
      DATA ((IY(I,J), J = 1,3), I = 1,6) /0,1,2,0,2,3,0,2,3,0,3,4,
     &                                    0,3,4,0,3,4/
C
      CALL TSTDAT(NONTYP,NONTVC,IQM,JCO,ICOUNT,NCOUNT,NDIFNU)
C
      ICOUN2 = ICOUNT*(ICOUNT + 1)/2
      NCOUN2 = NCOUNT*(NCOUNT + 1)/2
C
C     Read integrals
C
      DO 10 I = 1, 6
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABELT(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABELT(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,NCOUN2*IRAT,AUXVEC)
         CALL DCOPY(NCOUN2,AUXVEC,1,RESVEC(1,I),1)
 10   CONTINUE
      DO 15 I = 1, 3
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABEL(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABEL(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,ICOUN2*IRAT,AUXVEC)
         CALL DCOPY(ICOUN2,AUXVEC,1,SINTVC(1,I),1)
 15   CONTINUE
C
C     Calculation of the new integrals
C
      DIFFER = .FALSE.
      QMN(1,1) = D0
      QMN(2,2) = D0
      QMN(3,3) = D0
      CALL DZERO(RMXDIF,3)
      ICOUNR = 0
      NUCNM1 = 0
      JPRIMA = 0
      DO 30 I = 1, NDIFNU
         DO 30 II = 1, NONTVC(I)
         NUCNM1 = NUCNM1 + 1
         DO 30 J = 1, IQM(I)
         DO 30 K = 1, JCO(I,J)
         JPRIMA = JPRIMA + 1
         DO 30 L = 1, J*(J + 1)/2
            ICOUNR = ICOUNR + 1
            ICOUNL = 0
            NUCNM2 = 0
            DO 40 I1 = 1, NDIFNU
               DO 40 II1 = 1, NONTVC(I1)
               NUCNM2 = NUCNM2 +1
               DO 40 J1 = 1, IQM(I1)
               DO 40 K1 = 1, JCO(I1,J1)
               DO 40 L1 = 1, J1*(J1 + 1)/2
                  ICOUNL = ICOUNL + 1
                  IF (ICOUNR .GE. ICOUNL) THEN
                     IADR = ICOUNR*(ICOUNR - 1)/2 + ICOUNL
                  QMN(1,2) = CORD(3,NUCNM2) - CORD(3,NUCNM1)
                  QMN(2,1) = - QMN(1,2)
                  QMN(1,3) = CORD(2,NUCNM1) - CORD(2,NUCNM2)
                  QMN(3,1) = - QMN(1,3)
                  QMN(2,3) = CORD(1,NUCNM2) - CORD(1,NUCNM1)
                  QMN(3,2) = - QMN(2,3)
                  IADR2 = NCOUN2
                  NROW2 = NCOUNT
                  NNUC = 0
                  DO 100 I2 = NDIFNU + 1, NONTYP
                     DO 100 II2 = 1, NONTVC(I2)
                     NNUC = NNUC + 1
                     DO 100 J2 = 1, IQM(I2)
                        DO 100 K2 = 1, JCO(I2,J2)
                           DO  100 L2 = 1, J2*(J2 + 1)/2
                              NROW2 = NROW2 + 1
                              IF ((NUCNM1 .EQ. NNUC) .AND. (J2
     &                             .EQ. (J + 1)) .AND. (K .EQ. K2)
     &                             .AND. (L .EQ. L2)) THEN
                                 IADR3 = IADR2
                                 IADR4 = IADR2
                                 DO 110 NN = 1, IY(L,2)
                                    IADR3 = IADR3 + NROW2 + NN - 1
 110                             CONTINUE
                                 DO 120 NN = 1, IY(L,3)
                                    IADR4 = IADR4 + NROW2 + NN - 1
 120                             CONTINUE
                                 IADRX = IADR2 + ICOUNL
                                 IADRY = IADR3 + ICOUNL
                                 IADRZ = IADR4 + ICOUNL
                              ELSE
                                 IADR2 = IADR2 + NROW2
                              END IF
 100                       CONTINUE
                  STVEC(1,1) = SINTVC(IADRX,1)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(1,NUCNM1)*SINTVC(IADR,1)
                  STVEC(1,2) = SINTVC(IADRX,2)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(1,NUCNM1)*SINTVC(IADR,2)
                  STVEC(1,3) = SINTVC(IADRX,3)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(1,NUCNM1)*SINTVC(IADR,3)
                  STVEC(2,1) = SINTVC(IADRY,1)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(2,NUCNM1)*SINTVC(IADR,1)
                  STVEC(2,2) = SINTVC(IADRY,2)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(2,NUCNM1)*SINTVC(IADR,2)
                  STVEC(2,3) = SINTVC(IADRY,3)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(2,NUCNM1)*SINTVC(IADR,3)
                  STVEC(3,1) = SINTVC(IADRZ,1)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(3,NUCNM1)*SINTVC(IADR,1)
                  STVEC(3,2) = SINTVC(IADRZ,2)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(3,NUCNM1)*SINTVC(IADR,2)
                  STVEC(3,3) = SINTVC(IADRZ,3)/(D2*SQRT(PRIEXP(JPRIMA)))
     &                     + CORD(3,NUCNM1)*SINTVC(IADR,3)
                  CALL DGEMM('N','N',3,3,3,1.D0,
     &                       QMN,3,
     &                       STVEC,3,0.D0,
     &                       RSVEC,3)
                  INTNUM = 0
                  DO 50 I2 = 1, 3
                     DO 50 I3 = I2, 3
                        INTNUM = INTNUM + 1
                        TSTVEC(IADR,INTNUM) =
     &                       DP25*(RSVEC(I2,I3) +RSVEC(I3,I2))
                     DIFRS = RESVEC(IADR,INTNUM)
     &                     - TSTVEC(IADR,INTNUM)
                     DIFFER = DIFFER .OR. (ABS(DIFRS) .GT. PRTHRS)
                     RMXDIF(INTNUM) = MAX(RMXDIF(INTNUM),ABS(DIFRS))
                     DIFVEC(IADR,INTNUM) = DIFRS
 50               CONTINUE
            END IF
 40         CONTINUE
 30   CONTINUE
      IF (.NOT. DIFFER) THEN
         WRITE(LUPRI,'(/A/)')
     &         'No differences found in magnetic properties test'
      ELSE
         WRITE(LUPRI,'(/A/)')
     &         'Differences found in magnetic properties test!'
         WRITE(LUPRI,'(/A/)') 'Output of difference matrices follows'
         DO 300 I = 1, 6
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUN2,DIFVEC(1,I),1,AUXVEC,1)
            CALL OUTPAK(AUXVEC,NCOUNT,1,LUPRI)
 300     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &         'Output of test program integrals is as follows:'
         DO 320 I = 1, 6
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUN2,TSTVEC(1,I),1,AUXVEC,1)
            CALL OUTPAK(AUXVEC,NCOUNT,1,LUPRI)
 320     CONTINUE
      ENDIF
      RETURN
      END
C  /* Deck ns1tst */
      SUBROUTINE NS1TST(WORK,LWORK,IPRINT,DOATOM,NATOM)
#include <implicit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxmom.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <inforb.h>
C
      DIMENSION WORK(LWORK), INTREP(9*MXCENT), INTRE2(9*MXCENT)
      CHARACTER*8 LABEL(9*MXCENT), LABELT(9*MXCENT)
      LOGICAL DOATOM(NUCIND)
C
      KATOM  = 1
      CALL SETATM(WORK(KATOM),NATOM,10)
      INTTYP = 10
      CALL PSOTYP(NOPTYP,INTREP,LABEL,DOATOM,WORK(KATOM),NATOM,INTTYP)
      CALL NSTTYP(NOPTY2,INTRE2,LABELT,WORK(KATOM),NATOM,'NSLO',INTREP)
      NOPTY2 = NOPTY2/2
      NOPTYP = NOPTYP/2
      KIQM   = KATOM + NUCDEP
      KJCO   = KIQM   + (NUCIND + 1)/IRAT
      KAUXVE = KJCO   + (4*NUCIND + 1)/IRAT
      KTSTVE = KAUXVE + NBASIS*(NBASIS + 1)/2
      KRESVE = KTSTVE + NOPTY2*NBASIS*(NBASIS + 1)/2
      KPSO   = KRESVE + NOPTY2*NBASIS*(NBASIS + 1)/2
      KQMN   = KPSO + NOPTYP*NBASIS*(NBASIS + 1)/2
      KSTVEC = KQMN + 9
      KRSVEC = KSTVEC + 9
      KMXDIF = KRSVEC + 9
      KNONTV = KMXDIF + 3
      KDIFVC = KNONTV + (NUCIND + 1)/IRAT
      KLAST  = KDIFVC + 3*NBASIS*(NBASIS + 1)/2
      IF (KLAST .GT. LWORK) CALL STOPIT('MG1TST',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL NS1TS1(WORK(KIQM),WORK(KJCO),WORK(KAUXVE),
     &            WORK(KTSTVE),WORK(KRESVE),WORK(KQMN),LABEL,LABELT,
     &            WORK(KSTVEC),WORK(KRSVEC),WORK(KDIFVC),
     &            WORK(KMXDIF),WORK(KNONTV),WORK(KLAST),LWRK,IPRINT,
     &            WORK(KPSO),NOPTYP,NOPTY2)
      RETURN
      END
C  /* Deck ns1ts1 */
      SUBROUTINE NS1TS1(IQM,JCO,AUXVEC,TSTVEC,RESVEC,QMN,LABEL,
     &                  LABELT,STVEC,RSVEC,DIFVEC,MAXDIF,NONTVC,
     &                  WORK,LWORK,IPRINT,PSOVEC,NOPTYP,NOPTY2)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <primit.h>
#include <nuclei.h>
      PARAMETER (D0 = 0.D0, DP5 = 0.5D0)
      PARAMETER (PRTHRS = 1D-13, D2 = 2.D0)
      CHARACTER*8 LABEL(NOPTYP), LABELT(NOPTY2)
      LOGICAL FNDLAB, DIFFER
C
      DIMENSION WORK(LWORK), IQM(NUCIND), JCO(NUCIND,4),
     &          AUXVEC(NBASIS*(NBASIS+1)/2),
     &          TSTVEC(NBASIS*(NBASIS + 1)/6,NOPTY2),
     &          RESVEC(NBASIS*(NBASIS + 1)/6,NOPTY2),
     &          PSOVEC(NBASIS*(NBASIS + 1)/2,NOPTYP),
     &          QMN(3,3),STVEC(3,3),RSVEC(3,3),IY(6,3),RMXDIF(3),
     &          NONTVC(NUCIND), DIFVEC(NBASIS*(NBASIS + 1)/2,NOPTY2)
C
#include <inftap.h>
#include <shells.h>
      DATA ((IY(I,J), J = 1,3), I = 1,6) /0,1,2,0,2,3,0,2,3,0,3,4,
     &                                    0,3,4,0,3,4/
C
      CALL TSTDAT(NONTYP,NONTVC,IQM,JCO,ICOUNT,NCOUNT,NDIFNU)
C
      ICOUN2 = ICOUNT*(ICOUNT + 1)/2
      NCOUN2 = NCOUNT*(NCOUNT + 1)/2
C
C     Read integrals
C
      DO 10 I = 1, NOPTY2
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABELT(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABELT(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,NCOUN2*IRAT,AUXVEC)
         CALL DCOPY(NCOUN2,AUXVEC,1,RESVEC(1,I),1)
 10   CONTINUE
      DO 15 I = 1, NOPTYP
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABEL(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABEL(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,ICOUN2*IRAT,AUXVEC)
         CALL DCOPY(ICOUN2,AUXVEC,1,PSOVEC(1,I),1)
 15   CONTINUE
C
C     Calculation of the new integrals
C
      DIFFER = .FALSE.
      QMN(1,1) = D0
      QMN(2,2) = D0
      QMN(3,3) = D0
      CALL DZERO(RMXDIF,3)
      ICOUNR = 0
      NUCNM1 = 0
      JPRIMA = 0
      DO 30 I = 1, NDIFNU
         DO 30 II = 1, NONTVC(I)
         NUCNM1 = NUCNM1 + 1
         DO 30 J = 1, IQM(I)
         DO 30 K = 1, JCO(I,J)
         JPRIMA = JPRIMA + 1
         DO 30 L = 1, J*(J + 1)/2
            ICOUNR = ICOUNR + 1
            ICOUNL = 0
            NUCNM2 = 0
            DO 40 I1 = 1, NDIFNU
               DO 40 II1 = 1, NONTVC(I1)
               NUCNM2 = NUCNM2 +1
               DO 40 J1 = 1, IQM(I1)
               DO 40 K1 = 1, JCO(I1,J1)
               DO 40 L1 = 1, J1*(J1 + 1)/2
                  ICOUNL = ICOUNL + 1
                  IF (ICOUNR .GE. ICOUNL) THEN
                     IADR = ICOUNR*(ICOUNR - 1)/2 + ICOUNL
                  QMN(1,2) = CORD(3,NUCNM2) - CORD(3,NUCNM1)
                  QMN(2,1) = - QMN(1,2)
                  QMN(1,3) = CORD(2,NUCNM1) - CORD(2,NUCNM2)
                  QMN(3,1) = - QMN(1,3)
                  QMN(2,3) = CORD(1,NUCNM2) - CORD(1,NUCNM1)
                  QMN(3,2) = - QMN(2,3)
                  IADR2 = NCOUN2
                  NROW2 = NCOUNT
                  NNUC = 0
                  DO 100 I2 = NDIFNU + 1, NONTYP
                     DO 100 II2 = 1, NONTVC(I2)
                     NNUC = NNUC + 1
                     DO 100 J2 = 1, IQM(I2)
                        DO 100 K2 = 1, JCO(I2,J2)
                           DO  100 L2 = 1, J2*(J2 + 1)/2
                              NROW2 = NROW2 + 1
                              IF ((NUCNM1 .EQ. NNUC) .AND. (J2
     &                             .EQ. (J + 1)) .AND. (K .EQ. K2)
     &                             .AND. (L .EQ. L2)) THEN
                                 IADR3 = IADR2
                                 IADR4 = IADR2
                                 DO 110 NN = 1, IY(L,2)
                                    IADR3 = IADR3 + NROW2 + NN - 1
 110                             CONTINUE
                                 DO 120 NN = 1, IY(L,3)
                                    IADR4 = IADR4 + NROW2 + NN - 1
 120                             CONTINUE
                                 IADRX = IADR2 + ICOUNL
                                 IADRY = IADR3 + ICOUNL
                                 IADRZ = IADR4 + ICOUNL
                              ELSE
                                 IADR2 = IADR2 + NROW2
                              END IF
 100                       CONTINUE
            NUCLEU = 0
            DO 45 N = 1, NDIFNU
               DO 45 N2 = 1, NONTVC(N)
                  NUCLEU = NUCLEU + 1
                  ICOMPX = 3*(NUCLEU - 1) + 1
                  ICOMPY = 3*(NUCLEU - 1) + 2
                  ICOMPZ = 3*(NUCLEU - 1) + 3
                  STVEC(1,1) = PSOVEC(IADRX,ICOMPX)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                       + CORD(1,NUCNM1)*PSOVEC(IADR,ICOMPX)
                  STVEC(1,2) = PSOVEC(IADRX,ICOMPY)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                       + CORD(1,NUCNM1)*PSOVEC(IADR,ICOMPY)
                  STVEC(1,3) = PSOVEC(IADRX,ICOMPZ)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                       + CORD(1,NUCNM1)*PSOVEC(IADR,ICOMPZ)
                  STVEC(2,1) = PSOVEC(IADRY,ICOMPX)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                       + CORD(2,NUCNM1)*PSOVEC(IADR,ICOMPX)
                  STVEC(2,2) = PSOVEC(IADRY,ICOMPY)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                       + CORD(2,NUCNM1)*PSOVEC(IADR,ICOMPY)
                  STVEC(2,3) = PSOVEC(IADRY,ICOMPZ)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                       + CORD(2,NUCNM1)*PSOVEC(IADR,ICOMPZ)
                  STVEC(3,1) = PSOVEC(IADRZ,ICOMPX)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                       + CORD(3,NUCNM1)*PSOVEC(IADR,ICOMPX)
                  STVEC(3,2) = PSOVEC(IADRZ,ICOMPY)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                       + CORD(3,NUCNM1)*PSOVEC(IADR,ICOMPY)
                  STVEC(3,3) = PSOVEC(IADRZ,ICOMPZ)/
     &                         (D2*SQRT(PRIEXP(JPRIMA)))
     &                       + CORD(3,NUCNM1)*PSOVEC(IADR,ICOMPZ)
                  CALL DGEMM('N','N',3,3,3,1.D0,
     &                       QMN,3,
     &                       STVEC,3,0.D0,
     &                       RSVEC,3)
                  DO 50 I2 = 1, 3
                  DO 50 I3 = 1, 3
                     ICOMP = 9*(NUCLEU - 1) + 3*(I3 - 1) + I2
                     TSTVEC(IADR,ICOMP) = DP5*RSVEC(I2,I3)
                     DIFRS = RESVEC(IADR,ICOMP)
     &                     - TSTVEC(IADR,ICOMP)
                     DIFFER = DIFFER .OR. (ABS(DIFRS) .GT. PRTHRS)
                     RMXDIF(I2) = MAX(RMXDIF(I2),ABS(DIFRS))
                     DIFVEC(IADR,ICOMP) = DIFRS
 50               CONTINUE
 45         CONTINUE
            END IF
 40         CONTINUE
 30   CONTINUE
      IF (.NOT. DIFFER) THEN
         WRITE(LUPRI,'(/A/)')
     &         'No differences found in magnetic properties test'
      ELSE
         WRITE(LUPRI,'(/A/)')
     &         'Differences found in magnetic properties test!'
         WRITE(LUPRI,'(/A/)') 'Output of difference matrices follows'
         DO 300 I = 1, NOPTY2
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUN2,DIFVEC(1,I),1,AUXVEC,1)
            CALL OUTPAK(AUXVEC,NCOUNT,1,LUPRI)
 300     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &         'Output of test program integrals is as follows:'
         DO 320 I = 1, NOPTY2
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUN2,TSTVEC(1,I),1,AUXVEC,1)
            CALL OUTPAK(AUXVEC,NCOUNT,1,LUPRI)
 320     CONTINUE
      ENDIF
      RETURN
      END
C  /* Deck nstst2 */
      SUBROUTINE NSTST2(WORK,LWORK,IPRINT,DOATOM,NATOM)
#include <implicit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxmom.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <inforb.h>
C
      DIMENSION WORK(LWORK), INTREP(9*MXCENT), INTRE2(9*MXCENT)
      CHARACTER*8 LABEL(9*MXCENT), LABELT(9*MXCENT)
      LOGICAL DOATOM(NUCIND)
C
      KATOM  = 1
      CALL SETATM(WORK(KATOM),NATOM,10)
      CALL EFNTYP(NOPTYP,INTREP,LABEL,DOATOM,WORK(KATOM),NATOM,29)
      CALL NSTTYP(NOPTY2,INTRE2,LABELT,WORK(KATOM),NATOM,'NSNL',INTREP)
      NOPTY2 = NOPTY2/2
      NOPTYP = NOPTYP/2
      KIQM   = KATOM + NUCDEP
      KJCO   = KIQM  + (NUCIND + 1)/IRAT
      KAUXVE = KJCO   + (4*NUCIND + 1)/IRAT
      KTSTVE = KAUXVE + NBASIS*(NBASIS + 1)/2
      KRESVE = KTSTVE + NOPTY2*NBASIS*(NBASIS + 1)/2
      KNEF   = KRESVE + NOPTY2*NBASIS*(NBASIS + 1)/2
      KNONTV = KNEF + NOPTYP*NBASIS*(NBASIS + 1)/2
      KDIFVC = KNONTV + (NUCIND + 1)/IRAT
      KLAST  = KDIFVC + NOPTY2*NBASIS*(NBASIS + 1)/2
      IF (KLAST .GT. LWORK) CALL STOPIT('MG1TST',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL NS1TS2(WORK(KIQM),WORK(KJCO),WORK(KAUXVE),
     &            WORK(KTSTVE),WORK(KRESVE),LABEL,LABELT,WORK(KDIFVC),
     &            WORK(KNONTV),WORK(KLAST),LWRK,IPRINT,
     &            WORK(KNEF),NOPTYP,NOPTY2)
      RETURN
      END
C  /* Deck ns1ts2 */
      SUBROUTINE NS1TS2(IQM,JCO,AUXVEC,TSTVEC,RESVEC,LABEL,
     &                  LABELT,DIFVEC,NONTVC,
     &                  WORK,LWORK,IPRINT,EFVEC,NOPTYP,NOPTY2)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <primit.h>
#include <nuclei.h>
      PARAMETER (DP5 = 0.5D0)
      PARAMETER (PRTHRS = 1D-13, D2 = 2.D0)
      CHARACTER*8 LABEL(NOPTYP), LABELT(NOPTY2)
      LOGICAL FNDLAB, DIFFER
C
      DIMENSION WORK(LWORK), IQM(NUCIND), JCO(NUCIND,4),
     &          AUXVEC(NBASIS*(NBASIS+1)/2),
     &          TSTVEC(NBASIS*(NBASIS + 1)/6,NOPTY2),
     &          RESVEC(NBASIS*(NBASIS + 1)/6,NOPTY2),
     &          EFVEC(NBASIS*(NBASIS + 1)/2,NOPTYP),
     &          IY(6,3),NONTVC(NUCIND),
     &          DIFVEC(NBASIS*(NBASIS + 1)/2,NOPTY2)
C
#include <inftap.h>
#include <shells.h>
      DATA ((IY(I,J), J = 1,3), I = 1,6) /0,1,2,0,2,3,0,2,3,0,3,4,
     &                                    0,3,4,0,3,4/
C
      CALL TSTDAT(NONTYP,NONTVC,IQM,JCO,ICOUNT,NCOUNT,NDIFNU)
C
      ICOUN2 = ICOUNT*(ICOUNT + 1)/2
      NCOUN2 = NCOUNT*(NCOUNT + 1)/2
C
C     Read integrals
C
      DO 10 I = 1, NOPTY2
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABELT(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABELT(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,NCOUN2*IRAT,AUXVEC)
         CALL DCOPY(NCOUN2,AUXVEC,1,RESVEC(1,I),1)
 10   CONTINUE
      DO 15 I = 1, NOPTYP
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABEL(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABEL(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,ICOUN2*IRAT,AUXVEC)
         CALL DCOPY(ICOUN2,AUXVEC,1,EFVEC(1,I),1)
 15   CONTINUE
C
C     Calculation of the new integrals
C
      DIFFER = .FALSE.
      ICOUNR = 0
      NUCNM1 = 0
      JPRIMA = 0
      DO 30 I = 1, NDIFNU
         DO 30 II = 1, NONTVC(I)
         NUCNM1 = NUCNM1 + 1
         DO 30 J = 1, IQM(I)
         DO 30 K = 1, JCO(I,J)
         JPRIMA = JPRIMA + 1
         DO 30 L = 1, J*(J + 1)/2
            ICOUNR = ICOUNR + 1
            ICOUNL = 0
            NUCNM2 = 0
            DO 40 I1 = 1, NDIFNU
               DO 40 II1 = 1, NONTVC(I1)
               NUCNM2 = NUCNM2 +1
               DO 40 J1 = 1, IQM(I1)
               DO 40 K1 = 1, JCO(I1,J1)
               DO 40 L1 = 1, J1*(J1 + 1)/2
                  ICOUNL = ICOUNL + 1
                  IF (ICOUNR .GE. ICOUNL) THEN
                     IADR = ICOUNR*(ICOUNR - 1)/2 + ICOUNL
                  XCORD = CORD(1,NUCNM1) - CORD(1,NUCNM2)
                  YCORD = CORD(2,NUCNM1) - CORD(2,NUCNM2)
                  ZCORD = CORD(3,NUCNM1) - CORD(3,NUCNM2)
                  IADR2 = NCOUN2
                  NROW2 = NCOUNT
                  NNUC = 0
                  DO 100 I2 = NDIFNU + 1, NONTYP
                     DO 100 II2 = 1, NONTVC(I2)
                     NNUC = NNUC + 1
                     DO 100 J2 = 1, IQM(I2)
                        DO 100 K2 = 1, JCO(I2,J2)
                           DO  100 L2 = 1, J2*(J2 + 1)/2
                              NROW2 = NROW2 + 1
                              IF ((NUCNM1 .EQ. NNUC) .AND. (J2
     &                             .EQ. (J + 1)) .AND. (K .EQ. K2)
     &                             .AND. (L .EQ. L2)) THEN
                                 IADR3 = IADR2
                                 IADR4 = IADR2
                                 DO 110 NN = 1, IY(L,2)
                                    IADR3 = IADR3 + NROW2 + NN - 1
 110                             CONTINUE
                                 DO 120 NN = 1, IY(L,3)
                                    IADR4 = IADR4 + NROW2 + NN - 1
 120                             CONTINUE
                                 IADRX = IADR2 + ICOUNL
                                 IADRY = IADR3 + ICOUNL
                                 IADRZ = IADR4 + ICOUNL
                              ELSE
                                 IADR2 = IADR2 + NROW2
                              END IF
 100                       CONTINUE
            NUCLEU = 0
            DO 45 N = 1, NDIFNU
               DO 45 N2 = 1, NONTVC(N)
                  NUCLEU = NUCLEU + 1
                  ICOMPX = 3*(NUCLEU - 1) + 1
                  ICOMPY = 3*(NUCLEU - 1) + 2
                  ICOMPZ = 3*(NUCLEU - 1) + 3
                  ICOMP  = 9*(NUCLEU - 1)
                  TSTVEC(IADR,ICOMP + 1) = DP5*((EFVEC(IADRY,ICOMPY)
     &                                   +       EFVEC(IADRZ,ICOMPZ))/
     &                                     (D2*SQRT(PRIEXP(JPRIMA)))
     &                                   + YCORD*EFVEC(IADR,ICOMPY)
     &                                   + ZCORD*EFVEC(IADR,ICOMPZ))
                  TSTVEC(IADR,ICOMP + 2) = -DP5*(EFVEC(IADRX,ICOMPY)/
     &                                     (D2*SQRT(PRIEXP(JPRIMA)))
     &                                   + XCORD*EFVEC(IADR,ICOMPY))
                  TSTVEC(IADR,ICOMP + 3) = -DP5*(EFVEC(IADRX,ICOMPZ)/
     &                                     (D2*SQRT(PRIEXP(JPRIMA)))
     &                                   + XCORD*EFVEC(IADR,ICOMPZ))
                  TSTVEC(IADR,ICOMP + 4) = -DP5*(EFVEC(IADRY,ICOMPX)/
     &                                     (D2*SQRT(PRIEXP(JPRIMA)))
     &                                   + YCORD*EFVEC(IADR,ICOMPX))
                  TSTVEC(IADR,ICOMP + 5) = DP5*((EFVEC(IADRX,ICOMPX)
     &                                   +       EFVEC(IADRZ,ICOMPZ))/
     &                                     (D2*SQRT(PRIEXP(JPRIMA)))
     &                                   + XCORD*EFVEC(IADR,ICOMPX)
     &                                   + ZCORD*EFVEC(IADR,ICOMPZ))
                  TSTVEC(IADR,ICOMP + 6) = -DP5*(EFVEC(IADRY,ICOMPZ)/
     &                                     (D2*SQRT(PRIEXP(JPRIMA)))
     &                                   + YCORD*EFVEC(IADR,ICOMPZ))
                  TSTVEC(IADR,ICOMP + 7) = -DP5*(EFVEC(IADRZ,ICOMPX)/
     &                                     (D2*SQRT(PRIEXP(JPRIMA)))
     &                                   + ZCORD*EFVEC(IADR,ICOMPX))
                  TSTVEC(IADR,ICOMP + 8) = -DP5*(EFVEC(IADRZ,ICOMPY)/
     &                                     (D2*SQRT(PRIEXP(JPRIMA)))
     &                                   + ZCORD*EFVEC(IADR,ICOMPY))
                  TSTVEC(IADR,ICOMP + 9) = DP5*((EFVEC(IADRX,ICOMPX)
     &                                   +       EFVEC(IADRY,ICOMPY))/
     &                                     (D2*SQRT(PRIEXP(JPRIMA)))
     &                                   + XCORD*EFVEC(IADR,ICOMPX)
     &                                   + YCORD*EFVEC(IADR,ICOMPY))
                  DO 50 ICOMP2 = 1, 9
                     DIFRS = RESVEC(IADR,ICOMP + ICOMP2)
     &                     - TSTVEC(IADR,ICOMP + ICOMP2)
                     DIFFER = DIFFER .OR. (ABS(DIFRS) .GT. PRTHRS)
                     DIFVEC(IADR,ICOMP + ICOMP2) = DIFRS
 50               CONTINUE
 45         CONTINUE
            END IF
 40         CONTINUE
 30   CONTINUE
      IF (.NOT. DIFFER) THEN
         WRITE(LUPRI,'(/A/)')
     &         'No differences found in magnetic properties test'
      ELSE
         WRITE(LUPRI,'(/A/)')
     &         'Differences found in magnetic properties test!'
         WRITE(LUPRI,'(/A/)') 'Output of difference matrices follows'
         DO 300 I = 1, NOPTY2
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUN2,DIFVEC(1,I),1,AUXVEC,1)
            CALL OUTPAK(AUXVEC,NCOUNT,1,LUPRI)
 300     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &         'Output of test program integrals is as follows:'
         DO 320 I = 1, NOPTY2
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(NCOUN2,TSTVEC(1,I),1,AUXVEC,1)
            CALL OUTPAK(AUXVEC,NCOUNT,1,LUPRI)
 320     CONTINUE
      ENDIF
      RETURN
      END
C  /* Deck nstst3 */
      SUBROUTINE NSTST3(WORK,LWORK,IPRINT,DOATOM,NATOM)
#include <implicit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxmom.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <inforb.h>
C
      DIMENSION WORK(LWORK), INTREP(9*MXCENT), INTRE2(9*MXCENT)
      CHARACTER*8 LABEL(9*MXCENT), LABELT(9*MXCENT), LABELS(9*MXCENT)
      LOGICAL DOATOM(NUCIND)
C
      KATOM  = 1
      CALL SETATM(WORK(KATOM),NATOM,10)
      CALL NSTTYP(NOPTYP,INTREP,LABEL,WORK(KATOM),NATOM,'NSNL',INTREP)
      CALL NSTTYP(NOPTYP,INTREP,LABELS,WORK(KATOM),NATOM,'NSLO',INTREP)
      CALL NSTTYP(NOPTYP,INTRE2,LABELT,WORK(KATOM),NATOM,' NST',INTREP)
      KAUXVE = KATOM  + NUCDEP
      KTSTVE = KAUXVE + NBASIS*(NBASIS + 1)/2
      KRESVE = KTSTVE + NOPTYP*NBASIS*(NBASIS + 1)/2
      KNSNL  = KRESVE + NOPTYP*NBASIS*(NBASIS + 1)/2
      KNSLO  = KNSNL + NOPTYP*NBASIS*(NBASIS + 1)/2
      KDIFVC = KNSLO + NOPTYP*NBASIS*(NBASIS + 1)/2
      KLAST  = KDIFVC + NOPTYP*NBASIS*(NBASIS + 1)/2
      IF (KLAST .GT. LWORK) CALL STOPIT('NS1TS3',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL NS1TS3(WORK(KAUXVE),
     &            WORK(KTSTVE),WORK(KRESVE),LABEL,LABELS,LABELT,
     &            WORK(KDIFVC),WORK(KLAST),LWRK,IPRINT,
     &            WORK(KNSNL),WORK(KNSLO),NOPTYP)
      RETURN
      END
C  /* Deck ns1ts3 */
      SUBROUTINE NS1TS3(AUXVEC,TSTVEC,RESVEC,LABEL,LABELS,
     &                  LABELT,DIFVEC,WORK,LWORK,IPRINT,VCNSNL,
     &                  VCNSLO,NOPTYP)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <primit.h>
#include <nuclei.h>
      PARAMETER (D0 = 0.D0, DP5 = 0.5D0)
      PARAMETER (PRTHRS = 1D-13)
      CHARACTER*8 LABEL(NOPTYP), LABELT(NOPTYP), LABELS(NOPTYP)
      LOGICAL FNDLAB, DIFFER
C
      DIMENSION WORK(LWORK),
     &          AUXVEC(NBASIS*(NBASIS+1)/2),
     &          TSTVEC(NBASIS*(NBASIS + 1)/2,NOPTYP),
     &          RESVEC(NBASIS*(NBASIS + 1)/2,NOPTYP),
     &          VCNSNL(NBASIS*(NBASIS + 1)/2,NOPTYP),
     &          VCNSLO(NBASIS*(NBASIS + 1)/2,NOPTYP),
     &          DIFVEC(NBASIS*(NBASIS + 1)/2,NOPTYP)
C
#include <inftap.h>
#include <shells.h>
C
      ICOUN2 = NBASIS*(NBASIS + 1)/2
C
C     Read integrals
C
      DO 10 I = 1, NOPTYP
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABELT(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABELT(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,ICOUN2*IRAT,AUXVEC)
         CALL DCOPY(ICOUN2,AUXVEC,1,RESVEC(1,I),1)
 10   CONTINUE
      DO 15 I = 1, NOPTYP
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABEL(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABEL(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,ICOUN2*IRAT,AUXVEC)
         CALL DCOPY(ICOUN2,AUXVEC,1,VCNSNL(1,I),1)
 15   CONTINUE
      DO 17 I = 1, NOPTYP
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABELS(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &            'Label ', LABELS(I), ' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,ICOUN2*IRAT,AUXVEC)
         CALL DCOPY(ICOUN2,AUXVEC,1,VCNSLO(1,I),1)
 17   CONTINUE
C
C     Calculation of the new integrals
C
      DIFFER = .FALSE.
      DO 30 I = 1, ICOUN2
         DO 40 J = 1, NOPTYP
           TSTVEC(I,J) = VCNSNL(I,J) + VCNSLO(I,J)
           DIFRS = RESVEC(I,J)
     &           - TSTVEC(I,J)
           DIFFER = DIFFER .OR. (ABS(DIFRS) .GT. PRTHRS)
           DIFVEC(I,J) = DIFRS
 40      CONTINUE
 30   CONTINUE
      IF (.NOT. DIFFER) THEN
         WRITE(LUPRI,'(/A/)')
     &         'No differences found in magnetic properties test'
      ELSE
         WRITE(LUPRI,'(/A/)')
     &         'Differences found in magnetic properties test!'
         WRITE(LUPRI,'(/A/)') 'Output of difference matrices follows'
         DO 300 I = 1, NOPTYP
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(ICOUN2,DIFVEC(1,I),1,AUXVEC,1)
            CALL OUTPAK(AUXVEC,ICOUNT,1,LUPRI)
 300     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &         'Output of test program integrals is as follows:'
         DO 320 I = 1, NOPTYP
            WRITE(LUPRI,'(1X,A,I3)') 'Component:  ', I
            CALL DCOPY(ICOUN2,TSTVEC(1,I),1,AUXVEC,1)
            CALL OUTPAK(AUXVEC,ICOUNT,1,LUPRI)
 320     CONTINUE
      ENDIF
      RETURN
      END
C  /* Deck nptst */
      SUBROUTINE NPTST(WORK,LWORK,NATOM)
C
C     Test nuclear potential integrals, KR, Aug -92
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxmom.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
      PARAMETER (PRTHRS = 1D-13)
#include <inftap.h>
#include <nuclei.h>
#include <inforb.h>
C
      LOGICAL DIFFER, FNDLAB
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABINT(9*MXCENT)
C
      DIFFER  = .FALSE.
      KNPTRU = 1
      KNPTST = KNPTRU + NBASIS*(NBASIS + 1)/2
      KBUF   = KNPTST + NBASIS*(NBASIS + 1)/2
      KIBUF  = KBUF   + 600
      KATOM  = KIBUF  + 600
      KINTRE = KATOM  + NUCDEP
      KINTAD = KINTRE + (9*MXCENT + 1)/IRAT
      KAUXVE = KINTAD + (9*MXCENT + 1)/IRAT
      KDIFVE = KAUXVE + NBASIS*(NBASIS + 1)/2
      KLAST  = KDIFVE + NBASIS*(NBASIS + 1)/2
      IF (KLAST .GT. LWORK) CALL STOPIT('NPTST ',' ',KLAST,LWORK)
C
C     Read in correct integrals
C
      CALL RDINT1(WORK(KNPTRU),WORK(KBUF),WORK(KIBUF))
C
C     Read in test integrals
C
      CALL SETATM(WORK(KATOM),NATOM,35)
      CALL NPETYP(NOPTYP,WORK(KINTRE),WORK(KINTAD),LABINT,WORK(KATOM),
     &            NATOM)
      CALL DZERO(WORK(KNPTST),NBASIS*(NBASIS + 1)/2)
      DO 10 I = 1, NOPTYP
         REWIND LUPROP
         IF (.NOT. FNDLAB(LABINT(I),LUPROP)) THEN
            WRITE (LUPRI,'(1X,3A)')
     &           'Label ',LABINT(I),' not found on file AOPROPER'
            CALL QUIT('Error in magnetic properties test')
         END IF
         CALL READI(LUPROP,NBASIS*(NBASIS + 1)/2*IRAT,WORK(KAUXVE))
         DO 20 J = 0, NBASIS*(NBASIS + 1)/2 - 1
            WORK(KNPTST + J) = WORK(KNPTST + J) + WORK(KAUXVE + J)
 20      CONTINUE
 10   CONTINUE
C
      DO 30 J = 0, NBASIS*(NBASIS + 1)/2 - 1
         DIFRS = WORK(KNPTRU) - WORK(KNPTST)
         DIFFER = DIFFER .OR. (ABS(DIFRS) .GT. PRTHRS)
         WORK(KDIFVE + J) = DIFRS
 30   CONTINUE
      IF (.NOT. DIFFER) THEN
         WRITE (LUPRI,'(/A/)')
     &        ' No differences found in nuclear pot. energy test'
      ELSE
         WRITE (LUPRI,'(/A/)')
     &        ' Differences found in nuclear pot. energy test'
         WRITE (LUPRI,'(/A/)') 'Output of difference matrix follows'
         CALL OUTPAK(WORK(KDIFVE),NBASIS,1,LUPRI)
         WRITE (LUPRI,'(/A/)')
     &        ' Output of test program integrals follows'
         CALL OUTPAK(WORK(KNPTST),NBASIS,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck mm2tst */
      SUBROUTINE MM2TST(WORK,LWORK,PRTHRS,IPRINT)
C
C     K.Ruud, September 1992
C
#include <implicit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <inforb.h>
#include <cotabs.h>
C
      DIMENSION WORK(LWORK)
C
C     Allocate memory from work
C
      KIQM   = 1
      KJCO   = KIQM   + (NUCIND + 1)/IRAT
      KNONTV = KJCO   + (4*NUCIND + 1)/IRAT
      KBUF   = KNONTV + (NUCIND + 1)/IRAT
      IF (KBUF .GT. LWORK) CALL STOPIT('DS2BUF',' ',KBUF,LWORK)
C
C     Read necessary data from LUONEL
C
      CALL TSTDAT(NONTYP,WORK(KNONTV),WORK(KIQM),WORK(KJCO),ICOUNT,
     &            NCOUNT,NDIFNU)
C
C     More memory allocations
C
      KIBUF  = KBUF   + 600
      KQMN   = KIBUF  + 600
      KTSTVE = KQMN   + 9
      KTWOEL = KTSTVE + NCOUNT*NCOUNT*NCOUNT*NCOUNT*3
      KTSVE1 = KTWOEL + ICOUNT*ICOUNT*ICOUNT*ICOUNT
      KTSVE2 = KTSVE1 + 3
      KRESV1 = KTSVE2 + 3
      KRESV2 = KRESV1 + 3
      KDONE  = KRESV2 + 3
      KLAST  = KDONE  + (NCOUNT*NCOUNT*NCOUNT*NCOUNT*3 + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('DS2TS1',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL MM2TS1(WORK(KIQM),WORK(KJCO),WORK(KBUF),WORK(KIBUF),
     &            WORK(KNONTV),WORK(KQMN),WORK(KTSTVE),WORK(KTSVE1),
     &            WORK(KTSVE2),WORK(KRESV1),WORK(KRESV2),
     &            WORK(KTWOEL),WORK(KLAST),NONTYP,ICOUNT,NCOUNT,
     &            WORK(KDONE),LWRK,NDIFNU,PRTHRS,IPRINT)
      RETURN
      END
C  /* Deck mm2ts1 */
      SUBROUTINE MM2TS1(IQM,JCO,BUF,IBUF,NONTVC,QMN,TSTVEC,TSVEC1,
     &                  TSVEC2,RESVE1,RESVE2,TWOEL,WORK,NONTYP,ICOUNT,
     &                  NCOUNT,DONE,LWORK,NDIFNU,PRTHRS,IPRINT)
C
C     K.Ruud, September 1992
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxmom.h>
#include <primit.h>
#include <nuclei.h>
      PARAMETER (D0 = 0.D0, D2 = 2.D0, DP5 = 0.5D0)
      LOGICAL FNDLAB, DIFFER, FIRST, DONE(NCOUNT,NCOUNT,NCOUNT,NCOUNT,3)
      CHARACTER LABEL*8, COOR(3)*1
C
      DIMENSION IQM(NUCIND), JCO(NUCIND,4), BUF(600), IBUF(600),
     &          NONTVC(NUCIND), QMN(3,3), STVEC1(3),STVEC2(3),RSVEC1(3),
     &          RSVEC2(3),TSTVEC(NCOUNT,NCOUNT,NCOUNT,NCOUNT,3),
     &          TWOEL(ICOUNT,ICOUNT,ICOUNT,ICOUNT),WORK(LWORK),
     &          IY(6,3)
C
#include <inftap.h>
#include <shells.h>
      DATA COOR /'X','Y','Z'/
      DATA ((IY(I,J), J = 1, 3), I = 1, 6) /0,1,2,0,2,3,0,2,3,0,3,4,
     &                                      0,3,4,0,3,4/
C

C
C     Clear arrays and read in integrals
C
      CALL DZERO(TSTVEC,NCOUNT*NCOUNT*NCOUNT*NCOUNT*3)
      CALL DZERO(TWOEL ,ICOUNT*ICOUNT*ICOUNT*ICOUNT)
C
      CALL GPOPEN(LU1,'AO2MGINT','OLD',' ',' ',IDUMMY,.FALSE.)
C
C     Ordinary two-electron integrals
C
      CALL RDINTM(LUINTA,'BASTWOEL',TWOEL,DONE,BUF,IBUF,ICOUNT,1)
C
C     First-order differentiated two-el. integrals with respect to B
C
      CALL RDINTM(LU1,'AO2MGINT',TSTVEC,DONE,BUF,IBUF,NCOUNT,3)
C
C     Calculate new integrals
C
      DIFFER   = .FALSE.
      FIRST    = .TRUE.
      QMN(1,1) = D0
      QMN(2,2) = D0
      QMN(3,3) = D0
      ICOUR1   = 0
      NUCNM1   = 0
      JPRIM1   = 0
      DIFMAX   = D0
C
C     Loop over electron 1 indeces, row(30) and column(35)
C
      DO 30 I = 1, NDIFNU
      DO 30 II = 1, NONTVC(I)
      NUCNM1 = NUCNM1 + 1
      DO 30 J = 1, IQM(I)
      DO 30 K = 1, JCO(I,J)
      JPRIM1 = JPRIM1 + 1
      DO 30 L = 1, J*(J + 1)/2
         ICOUR1 = ICOUR1 + 1
         ICOUL1 = 0
         NUCNM2 = 0
         DO 35 I1 = 1, NDIFNU
         DO 35 II1 = 1, NONTVC(I1)
         NUCNM2 = NUCNM2 + 1
         DO 35 J1 = 1, IQM(I1)
         DO 35 K1 = 1, JCO(I1,J1)
         DO 35 L1 = 1, J1*(J1 + 1)/2
            ICOUL1 = ICOUL1 + 1
            ICOUR2 = 0
            NUCNM3 = 0
            JPRIM2 = 0
C
C     Loop over electron 2 indeces, row(40) and column(45)
C
            DO 40 I2 = 1, NDIFNU
            DO 40 II2 = 1, NONTVC(I2)
            NUCNM3 = NUCNM3 + 1
            DO 40 J2 = 1, IQM(I2)
            DO 40 K2 = 1, JCO(I2,J2)
            JPRIM2 = JPRIM2 + 1
            DO 40 L2 = 1, J2*(J2 + 1)/2
               ICOUR2 = ICOUR2 + 1
               ICOUL2 = 0
               NUCNM4 = 0
               DO 45 I3 = 1, NDIFNU
               DO 45 II3 = 1, NONTVC(I3)
               NUCNM4 = NUCNM4 + 1
               DO 45 J3 = 1, IQM(I3)
               DO 45 K3 = 1, JCO(I3,J3)
               DO 45 L3 = 1, J3*(J3 + 1)/2
                  ICOUL2 = ICOUL2 + 1
C
C     We first do the electron number 1 part
C
                  QMN(1,2) = CORD(3,NUCNM2) - CORD(3,NUCNM1)
                  QMN(1,3) = CORD(2,NUCNM1) - CORD(2,NUCNM2)
                  QMN(2,3) = CORD(1,NUCNM2) - CORD(1,NUCNM1)
                  QMN(2,1) = - QMN(1,2)
                  QMN(3,1) = - QMN(1,3)
                  QMN(3,2) = - QMN(2,3)
                  NROW2 = NCOUNT
                  DO 200 I5 = NDIFNU + 1, NONTYP
                     DO 200 II5 = 1, NONTVC(I5)
                     DO 200 J5 = 1, IQM(I5)
                        DO 200 K5 = 1, JCO(I5,J5)
                           DO  200 L5 = 1, J5*(J5 + 1)/2
                              NROW2 = NROW2 + 1
                              IF (((I+NDIFNU) .EQ. I5) .AND.
     &                             (J5 .EQ. (J + 1)) .AND. (K .EQ. K5)
     &                              .AND. (L .EQ. L5) .AND.
     &                             (II5 .EQ. II)) THEN
                  IADRX = NROW2
                  IADRY = NROW2 + IY(L,2)
                  IADRZ = NROW2 + IY(L,3)
                  STVEC1(1) = TWOEL(IADRX,ICOUL1,ICOUR2,ICOUL2)/
     &                        (D2*SQRT(PRIEXP(JPRIM1)))
     &                      + CORD(1,NUCNM1)*
     &                        TWOEL(ICOUR1,ICOUL1,ICOUR2,ICOUL2)
                  STVEC1(2) = TWOEL(IADRY,ICOUL1,ICOUR2,ICOUL2)/
     &                        (D2*SQRT(PRIEXP(JPRIM1)))
     &                      + CORD(2,NUCNM1)*
     &                        TWOEL(ICOUR1,ICOUL1,ICOUR2,ICOUL2)
                  STVEC1(3) = TWOEL(IADRZ,ICOUL1,ICOUR2,ICOUL2)/
     &                        (D2*SQRT(PRIEXP(JPRIM1)))
     &                      + CORD(3,NUCNM1)*
     &                        TWOEL(ICOUR1,ICOUL1,ICOUR2,ICOUL2)
                  CALL DGEMM('N','N',3,1,3,1.D0,
     &                       QMN,3,
     &                       STVEC1,3,0.D0,
     &                       RSVEC1,3)
               END IF
 200        CONTINUE
C
C     Second electron part
C
                  QMN(1,2) = CORD(3,NUCNM4) - CORD(3,NUCNM3)
                  QMN(1,3) = CORD(2,NUCNM3) - CORD(2,NUCNM4)
                  QMN(2,3) = CORD(1,NUCNM4) - CORD(1,NUCNM3)
                  QMN(2,1) = - QMN(1,2)
                  QMN(3,1) = - QMN(1,3)
                  QMN(3,2) = - QMN(2,3)
                  NROW2 = NCOUNT
                  DO 250 I6 = NDIFNU + 1, NONTYP
                     DO 250 II6 = 1, NONTVC(I6)
                     DO 250 J6 = 1, IQM(I6)
                        DO 250 K6 = 1, JCO(I6,J6)
                           DO  250 L6 = 1, J6*(J6 + 1)/2
                              NROW2 = NROW2 + 1
                              IF (((I2+NDIFNU) .EQ. I6) .AND.
     &                             (J6 .EQ. (J2 + 1)) .AND. (K6 .EQ. K2)
     &                              .AND. (L6 .EQ. L2) .AND.
     &                             (II6 .EQ. II2)) THEN
                  IADRX = NROW2
                  IADRY = NROW2 + IY(L2,2)
                  IADRZ = NROW2 + IY(L2,3)
                  STVEC2(1) = TWOEL(ICOUR1,ICOUL1,IADRX,ICOUL2)/
     &                        (D2*SQRT(PRIEXP(JPRIM2)))
     &                      + CORD(1,NUCNM3)*
     &                        TWOEL(ICOUR1,ICOUL1,ICOUR2,ICOUL2)
                  STVEC2(2) = TWOEL(ICOUR1,ICOUL1,IADRY,ICOUL2)/
     &                        (D2*SQRT(PRIEXP(JPRIM2)))
     &                      + CORD(2,NUCNM3)*
     &                        TWOEL(ICOUR1,ICOUL1,ICOUR2,ICOUL2)
                  STVEC2(3) = TWOEL(ICOUR1,ICOUL1,IADRZ,ICOUL2)/
     &                        (D2*SQRT(PRIEXP(JPRIM2)))
     &                      + CORD(3,NUCNM3)*
     &                        TWOEL(ICOUR1,ICOUL1,ICOUR2,ICOUL2)
                  CALL DGEMM('N','N',3,1,3,1.D0,
     &                       QMN,3,
     &                       STVEC2,3,0.D0,
     &                       RSVEC2,3)
               END IF
 250        CONTINUE
C
C     Construct test integrals, and compare them
C
                  DO 50 ICOMP = 1, 3
                     TEST = DP5*(RSVEC1(ICOMP) + RSVEC2(ICOMP))
                     DINT = TSTVEC(ICOUR1,ICOUL1,ICOUR2,ICOUL2,ICOMP)
                     DIFRS = DINT - TEST
                     ABSDF = ABS(DIFRS)
                     DIFFER = DIFFER .OR. (ABSDF .GT. PRTHRS)
                     IF (ABSDF .GT. DIFMAX) THEN
                        DIFMAX = ABSDF
                        DINTM  = DINT
                        TESTM  = TEST
                        IMAX   = ICOUR1
                        JMAX   = ICOUL1
                        KMAX   = ICOUR2
                        LMAX   = ICOUL2
                     END IF
                     IF (ABSDF .GT. PRTHRS) THEN
                        IF (FIRST) THEN
                           WRITE (LUPRI,*)
                           WRITE (LUPRI,'(A,D8.1,A)')
     &                          ' Differences found: (threshold,',
     &                          PRTHRS,')'
                           WRITE (LUPRI,'(A)')
     &                          '-----------------'
                           WRITE (LUPRI,'(17X,A)')
     &                          '  Magnetic moment integral    '//
     &                          'Test integral        Difference'
                           WRITE (LUPRI,*)
                           FIRST = .FALSE.
                        END IF
                        WRITE (LUPRI,'(1X,A,4I3,5X,2D23.15,D13.3)')
     &                       COOR(ICOMP),ICOUR1,ICOUL1,ICOUR2,ICOUL2,
     &                       DINT, TEST, DIFRS
                     END IF
 50               CONTINUE
 45            CONTINUE
 40         CONTINUE
 35      CONTINUE
 30   CONTINUE
      IF (.NOT. DIFFER) THEN
         WRITE (LUPRI,'(A)')
     &         ' No differences larger than threshold found.'
         WRITE (LUPRI,*)
      END IF
      WRITE (LUPRI,'(/,A,4I3,A,//,5X,D12.6,2(A,D12.6),/)')
     &         ' Largest difference found for element',
     &         IMAX,JMAX,KMAX,LMAX,': ',
     &         DINTM - TESTM, ' = ', DINTM ,' - ',TESTM
      CALL GPCLOSE(LU1,'KEEP')
      RETURN
      END
C  /* Deck rdintm */
      SUBROUTINE RDINTM(LUFILE,KEY,ARRINT,DONE,BUF,IBUF,IDIM,NCOMP)
#include <implicit.h>
#include <priunit.h>
      INTEGER P,Q,R,S
      LOGICAL DONE(IDIM,IDIM,IDIM,IDIM,NCOMP), NODUP
      CHARACTER*8 KEY
      DIMENSION ARRINT(IDIM,IDIM,IDIM,IDIM,*), BUF(600), IBUF(600)
#include <ibtfun.h>
C
      NODUP = .TRUE.
      IF (KEY .EQ. 'AO2MGINT') THEN
         DO 10 I1 = 1, IDIM
            DO 10 I2 = 1, IDIM
               DO 10 I3 = 1, IDIM
                  DO 10 I4 = 1, IDIM
                     DO 10 I5 = 1, NCOMP
                        DONE(I1,I2,I3,I4,I5) = .FALSE.
 10      CONTINUE
      END IF
      REWIND LUFILE
      CALL MOLLAB(KEY,LUFILE,LUPRI)
 150  READ (LUFILE, END=300) BUF,IBUF,LENGTH
      IF (LENGTH .GT. 0) THEN
         DO 100 I = 1, LENGTH
            LABEL = IBUF(I)
            VALUE = BUF(I)
            IF (KEY .EQ. 'BASTWOEL') THEN
               P = IBTAND(IBTSHR(LABEL,24),255)
               Q = IBTAND(IBTSHR(LABEL,16),255)
               R = IBTAND(IBTSHR(LABEL, 8),255)
               S = IBTAND(LABEL,255)
               ARRINT(P,Q,R,S,1) = VALUE
               ARRINT(P,Q,S,R,1) = VALUE
               ARRINT(Q,P,R,S,1) = VALUE
               ARRINT(Q,P,S,R,1) = VALUE
               ARRINT(R,S,P,Q,1) = VALUE
               ARRINT(R,S,Q,P,1) = VALUE
               ARRINT(S,R,P,Q,1) = VALUE
               ARRINT(S,R,Q,P,1) = VALUE
            ELSE
               S = IBTAND(LABEL,255)
               IF (S .EQ. 0) THEN
                  ICOMP = IBTAND(IBTSHR(LABEL,8),255)
               ELSE
                  P = IBTAND(IBTSHR(LABEL,24),255)
                  Q = IBTAND(IBTSHR(LABEL,16),255)
                  R = IBTAND(IBTSHR(LABEL, 8),255)
                  IF (DONE(P,Q,R,S,ICOMP)) THEN
                     WRITE (LUPRI,'(A,4I4)')
     &                 ' The following integral is duplicated: ',P,Q,R,S
                     NODUP = .FALSE.
                  END IF
                  IF (DONE(Q,P,S,R,ICOMP)) THEN
                     WRITE (LUPRI,'(A,4I4)')
     &                 ' The following integral is duplicated: ',Q,P,S,R
                     NODUP = .FALSE.
                  END IF
                  IF (DONE(R,S,P,Q,ICOMP)) THEN
                     WRITE (LUPRI,'(A,4I4)')
     &                 ' The following integral is duplicated: ',R,S,P,Q
                     NODUP = .FALSE.
                  END IF
                  IF (DONE(S,R,Q,P,ICOMP)) THEN
                     WRITE (LUPRI,'(A,4I4)')
     &                 ' The following integral is duplicated: ',S,R,Q,P
                     NODUP = .FALSE.
                  END IF
                  ARRINT(P,Q,R,S,ICOMP) = VALUE
                  ARRINT(R,S,P,Q,ICOMP) = VALUE
                  ARRINT(Q,P,S,R,ICOMP) = -VALUE
                  ARRINT(S,R,Q,P,ICOMP) = -VALUE
                  DONE(P,Q,R,S,ICOMP) = .TRUE.
                  DONE(Q,P,S,R,ICOMP) = .TRUE.
                  DONE(R,S,P,Q,ICOMP) = .TRUE.
                  DONE(S,R,Q,P,ICOMP) = .TRUE.
               END IF
            END IF
 100     CONTINUE
         ELSE IF (LENGTH .LT. 0) THEN
            GOTO 300
         END IF
         GOTO 150
 300  CONTINUE
      IF (KEY .EQ. 'AO2MGINT') THEN
         IF (NODUP) THEN
            WRITE (LUPRI,'(/,A,/)') ' No duplicates have been found.'
         ELSE
            WRITE (LUPRI,'(/,A,/)') ' Duplicates have been found.'
         END IF
      END IF
      RETURN
      END
C  /* Deck tstdat */
      SUBROUTINE TSTDAT(NONTYP,NONTVC,IQM,JCO,ICOUNT,NCOUNT,NDIFNU)
C
C     Resided previously inside magnetic properties test routines.
C     K.Ruud, sept 1992
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
      LOGICAL FNDLAB
      DIMENSION NONTVC(NUCIND), IQM(NUCIND), JCO(NUCIND,4)
#include <inftap.h>
C
      REWIND (LUONEL)
      IF (.NOT. FNDLAB('TESTDATA',LUONEL)) THEN
         WRITE (LUPRI,'(/A/)')
     &        ' Label TESTDATA not found on file AOONEINT'
         CALL QUIT('Error in magnetic properties test')
      END IF
      READ (LUONEL) NONTYP, (NONTVC(I), I = 1, NONTYP),
     &              (IQM(I), (JCO(I,J), J = 1, IQM(I)), I = 1, NONTYP)
C
C     Number of basis functions
C
      ICOUNT = 0
      NCOUNT = 0
      NDIFNU = NONTYP/2
      DO 5 I = 1, NONTYP
         DO 5 II = 1, NONTVC(I)
            DO 5 J = 1, IQM(I)
               DO 5 K = 1, JCO(I,J)
                  DO 5 L = 1, J*(J + 1)/2
                     ICOUNT = ICOUNT + 1
                     IF (I .LE. NDIFNU) NCOUNT = NCOUNT + 1
 5    CONTINUE
      RETURN
      END
