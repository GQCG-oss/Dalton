C  /* Deck cc_chopt */
      SUBROUTINE CC_CHOPT(FOCKD,T1AM,WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C     CCSD(T) energy correction with Cholesky decomposition.
C
#include <implicit.h>
#include <priunit.h>
#include <cc_cho.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CHARACTER*7 FILNAM, FILVI1, FILVI2
      CHARACTER*10 MODFIL
C
      DIMENSION FOCKD(*),T1AM(*)
C
      DIMENSION WORK(LWORK)
      DIMENSION CHOELE(MAXCHO),IOFSET(8),IOFF2(8),IOFF3(8),IOFF4(8)
      DIMENSION IOFF5(8)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <dummy.h>
C
      LOGICAL BCONV,ACONV,ICONV
C
      LOGICAL FBATCH
      PARAMETER (FBATCH = .FALSE.)
      PARAMETER (NONA = 2, NONB = 3)
C
      PARAMETER (FILVI1 = 'CHO_VI1', FILVI2 = 'CHO_VI2')
C
      PARAMETER (INFO = 3)
C
      CALL QENTER('CC_CHOPT')
C
      TITOT = SECOND()
      CALL GETTIM(TSTART,WSTART)
C 
C De-commented by Domenico
      CALL DZERO(SCNDSA,MAXCHO)
      CALL DZERO(SCNDSB,MAXCHO)
      CALL DZERO(SCNDSC,MAXCHO)
      CALL DZERO(SCNDSD,MAXCHO)
      CALL DZERO(SCNDSE,MAXCHO)
      CALL DZERO(SCNDSF,MAXCHO)
      CALL DZERO(SCNDSG,MAXCHO)
      CALL DZERO(SCNDSH,MAXCHO)
      CALL DZERO(SCNDSI,MAXCHO)
      CALL DZERO(SCNDSJ,MAXCHO)
      CALL DZERO(SCNDSX,MAXCHO)

      TLAST = SECOND()
C 
C De-commented by Domenico 
      CALL DZERO(ENERGA,MAXCHO)
      CALL DZERO(ENERGB,MAXCHO)
      CALL DZERO(ENERGC,MAXCHO)
      CALL DZERO(ENERGD,MAXCHO)
      CALL DZERO(ENERGE,MAXCHO)
      CALL DZERO(ENERGF,MAXCHO)
      CALL DZERO(ENERGG,MAXCHO)
      CALL DZERO(ENERGH,MAXCHO)
      CALL DZERO(ENERGI,MAXCHO)
      CALL DZERO(ENERGJ,MAXCHO)
      CALL DZERO(ENERGX,MAXCHO)
C
      CALL AROUND('Output from CC_CHOPT')
C
      LUCHOU = -1
      CALL GPOPEN(LUCHOU,'CHOPT_RST','UNKNOWN',' ','FORMATTED',
     &                     IDUM,LDUM)
      REWIND(LUCHOU)
C
      WRITE(LUCHOU,'(/,10X,A,/,10X,A,//)')
     &             'Cholesky CCSD(T) restart information',
     &             '------------------------------------'
      CALL FLSHFO(LUCHOU)
C
      IF (FBATCH) THEN
         WRITE(LUPRI,'(A)')
     &   '***NOTICE*** FBATCH=.TRUE.: Forcing batching for debug!'
         WRITE(LUPRI,'(A,I2,1X,I2,/)')
     &   'Enforced number of virtuals (B,A): ',NONB,NONA
      ENDIF
C
      CALL FLSHFO(LUPRI)
C
      ISYMT2 = 1
C
      TIMEX = ZERO
      TIMEZ = ZERO
      TIMA1 = ZERO
      TIMA2 = ZERO
      TIMEB = ZERO
      TIME1 = ZERO
      TIME2 = ZERO
      TIMEJ = ZERO
      TIMEI = ZERO
      TIMG1 = ZERO
      TIMEG = ZERO
      TIMED = ZERO
      TIMEC = ZERO
C
      E4AT = ZERO
      E4BT = ZERO
      E4C1 = ZERO
      E4C2 = ZERO
      E4D  = ZERO
      E4E1 = ZERO
      E4E2 = ZERO
      E4F1 = ZERO
      E4F2 = ZERO
      E5G  = ZERO
      E5G1 = ZERO
      E5H  = ZERO
      E5H1 = ZERO
      E4I  = ZERO
      E4J  = ZERO
C
      EYS  = ZERO
      EXR  = ZERO
      EWQ  = ZERO
      EVP  = ZERO
C
C--------------------
C     Open new files.
C--------------------
C
      LUVI1 = 0
      LUVI2 = 0
      CALL WOPEN2(LUVI1,FILVI1,64,0)
      CALL WOPEN2(LUVI2,FILVI2,64,0)
C
C-----------------------------------------------------------
C     Allocation for (jb|ki) integrals (sorted as (jk,i,b)).
C-----------------------------------------------------------
C
      KOINT = 1
      KEND0 = KOINT + NTRAOC(1)
      LWRK0 = LWORK - KEND0 + 1
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Insufficient memory in CC_CHOPT'
         WRITE(LUPRI,*) 'Need (*much* more than): ',KEND0-1
         WRITE(LUPRI,*) 'Available              : ',LWORK
         CALL QUIT(' Insufficient memory in CC_CHOPT ')
      ENDIF
C
C------------------------------------------------------------
C     Get virtual integrals (ai|bc) on disk, (jb|ki) in core.
C------------------------------------------------------------
C
      CALL CCHO_TRAINT(WORK(KEND0),LWRK0,WORK(KOINT),LUVI1,FILVI1,
     &                 LUVI2,FILVI2)
C
c     xnorm = DSQRT(DDOT(NTRAOC(1),WORK(KOINT),1,WORK(KOINT),1))
c      write(LUPRI,*) 'After CCHO_TRAINT: Norm of (ia|jk): ',xnorm
c      call flshfo(lupri)
C
C----------------------------------------------
C     Initial calculation of Cholesky ordering.
C----------------------------------------------
C
      CALL CCHO_INIT(FOCKD,CHOELE,MAXCHO,NUMCHO,WORK(KEND0),LWRK0)
C
      IF (NUMCHO .GT. MXCHVE) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
     &   'OBS!! ',NUMCHO,' vectors required for convergence'
         WRITE(LUPRI,*)
     &   'But only ',MXCHVE,' included as demanded in input'
         WRITE(LUPRI,*)
         NUMCHO = MXCHVE
      END IF
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KT2VO = KEND0
      KIAJB = KT2VO + NT2SQ(ISYMT2)
      KEND1 = KIAJB + MAX(NT2AM(ISYMT2),NT2AM(ISYMOP))
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in CC_CHOPT')
      ENDIF

C
C-------------------------------
C     Read amplitudes from disk.
C-------------------------------
C
      IOPT = 2
      CALL CC_RDRSP('R0',1,1,IOPT,MODFIL,DUMMY,WORK(KIAJB))
C
c     write(lupri,*) 'isymt2,NT2am(ISYMT2)',isymt2,NT2am(ISYMT2)
c     tnorm = dnorm2(NT2am(ISYMT2),WORK(KIAJB),1)
c     write(LUPRI,*) 'Before CCHO_T2VO: Norm of IAJB: ',tnorm
c     do i = 1,nt2am(isymt2)
c        write(lupri,*) i,work(kiajb+i-1)
c     end do
C---------------------------------
C     Sort amplitudes as t(ab,ij).
C---------------------------------
C
      CALL CCHO_T2A(WORK(KT2VO),WORK(KIAJB),ISYMT2) 

c     t1nrm = dnorm2(NT1AMX,T1AM,1)
c     tnorm = dnorm2(NT2SQ(ISYMT2),WORK(KT2VO),1)
c     write(LUPRI,*) 'After CCHO_T2VO: Norm of T1AM: ',t1nrm
c     write(LUPRI,*) 'After CCHO_T2VO: Norm of T2VO: ',tnorm
C
C------------------------------------------------
C     Read (ia|jb) integrals - sorted as (ai|bj).
C     Set up 2 Coulomb minus exchange in place.
C     NB: the integrals are packed as ai .le. bj
C         and unit 52 must be open!
C------------------------------------------------
C
      CALL CCG_RDIAJB(WORK(KIAJB),NT2AM(ISYMOP))
      IOPT = 1
      CALL CCSD_TCMEPK(WORK(KIAJB),ONE,ISYMOP,IOPT)

c     znorm = DSQRT(DDOT(NT2AM(ISYMOP),WORK(KIAJB),1,WORK(KIAJB),1))
c     write(LUPRI,*) 'After CCHO_RDIAJB: Norm of L(iajb): ',znorm,
c    &           ' (packed)'
C
C---------------------------------------
C     Start loop over virtual integrals.
C---------------------------------------
C
      MXCKA = 0
      DO ISYM = 1,NSYM
         IF (NCKATR(ISYM) .GT. MXCKA) MXCKA = NCKATR(ISYM)
      END DO
C
      LWORKB = LWRK1/2
C
      IF (IPRINT .GE. INFO) TIMEV = SECOND()
C
C     If restarting, go to the appropriate section
C
      C4ORO = OLD4O
      C5ORO = OLD5O
C
      IF (RSTH.OR.RSTH1.OR.RSTF1.OR.RSTC1.OR.RSTC2.OR.RSTVIR) THEN
         WRITE(LUPRI,'(//,A,/,A,//)') 
     &        '    Restarted CHO(T) calculation',
     &        '    ----------------------------'
         IF (UKNE4V) THEN
            WRITE(LUPRI,'(A,A)') 'Previous 4-V correction unknown. ',
     &                           'Be careful !!!'
         ELSE
            WRITE(LUPRI,'(A,D20.10)') '4-V correction : ',OLD4V
         END IF
         IF (UKNE5V) THEN
            WRITE(LUPRI,'(A,A)') 'Previous 5-V correction unknown. ',
     &                           'Be careful !!!'
         ELSE
            WRITE(LUPRI,'(A,D20.10)') '5-V correction : ',OLD5V
         END IF
      END IF
C
      IF (RSTH1 .OR. RSTF1 .OR. RSTC1 .OR. RSTC2) THEN
         IF (UKNE5O) THEN
            WRITE(LUPRI,'(A,A)') 'Previous 5-O correction unknown. ',
     &                           'Be careful !!!'
         ELSE 
            WRITE(LUPRI,'(A,D20.10)') '5-O correction : ',OLD5O
         END IF
      END IF
C
      IF (RSTC1 .OR. RSTC2) THEN
         IF (UKNE5O) THEN
            WRITE(LUPRI,'(A,A)') 'Previous 4-O correction unknown. ',
     &                           'Be careful !!!'
         ELSE 
            WRITE(LUPRI,'(A,D20.10)') '4-O correction : ',OLD4O
         END IF
      END IF
C
      IF (RSTH) THEN 
         WRITE(LUPRI,'(//,A,//)') 'Restarting from H term'
         GOTO 222
      ELSE IF (RSTH1) THEN
         WRITE(LUPRI,'(//,A,//)') 'Restarting from H1 term'
         GOTO 223
      ELSE IF (RSTF1) THEN
         WRITE(LUPRI,'(//,A,//)') 'Restarting from F1 term'
         GOTO 224
      ELSE IF (RSTC1) THEN
         WRITE(LUPRI,'(//,A,//)') 'Restarting from C1 term'
         GOTO 225
      ELSE IF (RSTC2) THEN
         WRITE(LUPRI,'(//,A,//)') 'Restarting from C2 term'
         GOTO 226
      ELSE IF (RSTVIR) THEN
         WRITE(LUPRI,'(//,A,/,A,I4,A,I2,/)')
     &       'Restarting virtual orbital loop ',
     &       'from orbital',IFVIOR,' of symmetry',IFVISY
      END IF
C
      DO 100 ISYMB = IFVISY,NSYM
C
         IF (RSTVIR) THEN
            NVIRB = NVIR(ISYMB) - IFVIOR + 1
         ELSE
            NVIRB = NVIR(ISYMB)
         END IF
C
         IF (NVIRB .EQ. 0) GOTO 100
C
         ISYEMD = MULD2H(ISYMB,ISYMOP)
         ISYIJC = MULD2H(ISYMB,ISYMT2)
         LWRK   = LWORKB - NCKATR(ISYEMD) - MXCKA
C
         LENB  = 4*NCKATR(ISYEMD) + NMAIJA(ISYIJC)
     &         + 2*NCKI(ISYMB)
         IF (FBATCH) THEN
            NEFB = MIN(NVIR(ISYMB),NONB)
            LEFF = NEFB*LENB + 1
            LWRK = MIN(LWRK,LEFF)
         ENDIF
C
         NUMB  = MIN(NVIRB,LWRK/LENB)
         IF (NUMB .EQ. 0) THEN
            WRITE(LUPRI,*) 'NUMB .EQ. 0 in CC_CHOPT'
            CALL QUIT('Not enough space in CC_CHOPT')
         END IF
C
         NBATB = (NVIRB-1)/NUMB + 1
C
         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(3X,A,I1,A,/,3X,A)')
     &      'Batch over B, symmetry ',ISYMB,':',
     &      '-------------------------'
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10)')
     &      'Minimum work space required   : ',LENB,
     &      'Work space available for batch: ',LWRK
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10,/)')
     &      'Number of virtual orbitals    : ',NVIR(ISYMB),
     &      'Required number of B-batches  : ',NBATB
         ENDIF
C
         IB2 = IFVIOR - 1
         DO IBATB = 1,NBATB
C
            IB1 = IB2 + 1
            IB2 = IB2 + NUMB
            IF (IB2 .GT. NVIR(ISYMB)) IB2 = NVIR(ISYMB)
            NUMIB = IB2 - IB1 + 1
C
            IF (IPRINT .GE. INFO) THEN
               WRITE(LUPRI,'(6X,A,I10,A,/,6X,A)')
     &         'B-batch number ',IBATB,':',
     &         '--------------------------'
               WRITE(LUPRI,'(6X,A,I10,1X,I10,/)')
     &         'First and last B: ',IB1,IB2
               CALL FLSHFO(LUPRI)
            ENDIF
C
C------------------------------
C           Dynamic allocation.
C------------------------------
C
            KIINT  = KEND1
            KJINT  = KIINT  + NCKATR(ISYEMD)*NUMIB
            KZMAT  = KJINT  + NCKATR(ISYEMD)*NUMIB
            KMXCKA = KZMAT  + NCKATR(ISYEMD)
            KTBIJC = KMXCKA + MXCKA
            KJTIL1 = KTBIJC + NMAIJA(ISYIJC)*NUMIB
            KJTIL2 = KJTIL1 + NCKATR(ISYEMD)*NUMIB
            KT2AM3 = KJTIL2 + NCKATR(ISYEMD)*NUMIB
            KT2AM4 = KT2AM3 + NCKI(ISYMB)*NUMIB
            KENDB1 = KT2AM4 + NCKI(ISYMB)*NUMIB
            LWRKB1 = LWORK  - KENDB1
            IF (LWRKB1 .LT. 0) THEN
               CALL QUIT('Insufficient memory in CC_CHOPT')
            ENDIF
C
C-------------------------------------------------------
C           Set work space available after batch over a.
C           (Only I and J integrals needed in memory).
C-------------------------------------------------------
C
            KENDB2 = KZMAT
            LWRKB2 = LWORK - KENDB2
C
C--------------------------
C           Read integrals.
C--------------------------
C
C           I(em,d,b) = (bm | ed)
C           ---------------------
C
            CALL CCHO_RDINT(LUVI2,FILVI2,WORK(KIINT),IB1,NUMIB,ISYMB)
C
C           J(em,d,b) = (em | bd)
C           ---------------------
C
            CALL CCHO_RDINT(LUVI1,FILVI1,WORK(KJINT),IB1,NUMIB,ISYMB)

c     len   = NCKATR(ISYMB)*NUMIB
c     xnorm = DSQRT(DDOT(len,WORK(KJINT),1,WORK(KJINT),1))
c     xinorm= DSQRT(DDOT(len,WORK(KIINT),1,WORK(KIINT),1))
c     write(LUPRI,*) 'After CCHO_RDINT: IB1, NUMIB, ISYMB: ',IB1,NUMIB,
c    &                ISYMB
c     write(LUPRI,*) 'After CCHO_RDINT: norm of J(em,d;#b)=(em|#bd): ',
c    &                xnorm
c     write(LUPRI,*) 'After CCHO_RDINT: norm of I(em,d;#b)=(#bm|ed): ',
c    &           xinorm
C
C-----------------------------------------
C           Resort J integrals for I-term.
C           Set up index array IOFF2.
C-----------------------------------------
C
            CALL CCHO_JRSRT(WORK(KJINT),WORK(KJTIL1),WORK(KJTIL2),
     &                      NUMIB,ISYMB,IOFF2)
C
C----------------------------------------------
C           Extract 2CME amplitudes for I-term.
C           Set up index array IOFF3.
C----------------------------------------------
C
            CALL CCHO_IX2CME(WORK(KT2VO),WORK(KT2AM3),WORK(KT2AM4),
     &                       IB1,NUMIB,ISYMB,IOFF3)
C
C-----------------------
C           Loop over a.
C-----------------------
C
            DO 200 ISYMA = 1,NSYM

c     write(LUPRI,*) 'entering symmetry A = ',ISYMA,
c    &           '; NVIR(ISYMA) = ',NVIR(ISYMA)
C
               IF (NVIR(ISYMA) .EQ. 0) GOTO 200
C
               ISYEMC = MULD2H(ISYMA,ISYMOP)
               ISYCKI = MULD2H(ISYMA,ISYMT2)
               ISYMAB = MULD2H(ISYMA,ISYMB)
               ISYMIJ = MULD2H(ISYMAB,ISYMT2)
               ISYMCD = ISYMAB
C
               IF (NMATAB(ISYMAB) .EQ. 0) GOTO 200
               IF (NMATIJ(ISYMIJ) .EQ. 0) GOTO 200
C
C---------------------------------------
C              Find largest ij subblock.
C---------------------------------------
C
               MAXIJ = -1
               DO ISYMJ = 1,NSYM
                  ISYMI = MULD2H(ISYMJ,ISYMIJ)
                  MAXIJ = MAX(MAXIJ,NRHF(ISYMI)*NRHF(ISYMJ))
               ENDDO
C
C-----------------------------------------------------
C              Allocate space for one Cholesky vector.
C-----------------------------------------------------
C
               LENCHO = NMATIJ(ISYMIJ)*NUMIB 
     &                + NCKATR(ISYEMC) + NMAIJA(ISYIJC) 
     &                + NCKI(ISYCKI)
C
               KCHOVE = KENDB1
               KENDCH = KCHOVE + LENCHO
               LWRKCH = LWORK  - KENDCH
               IF (LWRKCH .LT. 0)
     &         CALL QUIT('Not enough space in CC_CHOPT')
C
C---------------------------------
C              Batch over a index.
C---------------------------------
C
               LWORKA = LWRKCH
C
               LENA = 3*NCKATR(ISYEMC) + 3*NMATAB(ISYMCD)
     &              + NMATIJ(ISYMIJ)*NUMIB + 2*NUMIB
     &              + NMAIJA(ISYIJC) + 4*MAXIJ*NUMIB
     &              + 2*NCKI(ISYMA)
               IF (FBATCH) THEN
                  NEFA   = MIN(NVIR(ISYMA),NONA)
                  LEFF   = NEFA*LENA + 1
                  LWORKA = MIN(LWORKA,LEFF)
               ENDIF
               NUMA  = MIN(NVIR(ISYMA),LWORKA/LENA)
               IF (NUMA .EQ. 0) THEN
                  WRITE(LUPRI,*) 'NUMA .EQ. 0 in CC_CHOPT'
                  CALL QUIT('Not enough space in CC_CHOPT')
               END IF
C
               NBATA = (NVIR(ISYMA)-1)/NUMA + 1
C
               IF (IPRINT .GE. INFO) THEN
                  WRITE(LUPRI,'(9X,A,I1,A,/,9X,A)')
     &            'Batch over A, symmetry ',ISYMA,':',
     &            '-------------------------'
                  WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10)')
     &            'Minimum work space required   : ',LENA,
     &            'Work space available for batch: ',LWORKA
                  WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10,/)')
     &            'Number of virtual orbitals    : ',NVIR(ISYMA),
     &            'Required number of A-batches  : ',NBATA
               ENDIF
C
               IA2 = 0
               DO IBATA = 1,NBATA
C
                  IA1 = IA2 + 1
                  IA2 = IA2 + NUMA
                  IF (IA2 .GT. NVIR(ISYMA)) IA2 = NVIR(ISYMA)
                  NUMIA = IA2 - IA1 + 1
C
                  IF (IPRINT .GE. INFO) THEN
                     WRITE(LUPRI,'(12X,A,I10,A,/,12X,A)')
     &               'A-batch number ',IBATA,':',
     &               '--------------------------'
                     IF (ISYMA .NE. ISYMB) THEN
                        WRITE(LUPRI,'(12X,A,I10,1X,I10,/,12X,A,/)')
     &                  'First and last A: ',IA1,IA2,
     &              '(The A term will not be calculated for this batch)'
                     ELSE
                        WRITE(LUPRI,'(12X,A,I10,1X,I10,/)')
     &                  'First and last A: ',IA1,IA2
                     ENDIF
                     CALL FLSHFO(LUPRI)
                  ENDIF
C
C------------------------------------
C                 Dynamic allocation.
C------------------------------------
C
                  KTIJAB = KENDCH
                  KKINT  = KTIJAB + NMATIJ(ISYMIJ)*NUMIA*NUMIB
                  KLINT  = KKINT  + NCKATR(ISYEMC)*NUMIA
                  KXMAT  = KLINT  + NCKATR(ISYEMC)*NUMIA
                  KY1MAT = KXMAT  + NMATAB(ISYMCD)*NUMIA
                  KY2MAT = KY1MAT + NMATAB(ISYMCD)*NUMIA
                  KAMAT  = KY2MAT + NMATAB(ISYMCD)*NUMIA
                  KXAMA  = KAMAT  + NUMIA*NUMIB
                  KTIJCA = KXAMA  + NUMIA*NUMIB
                  KKTILD = KTIJCA + NMAIJA(ISYIJC)*NUMIA
                  KT2AM1 = KKTILD + NCKATR(ISYEMC)*NUMIA
                  KT2AM2 = KT2AM1 + NCKI(ISYMA)*NUMIA
                  KVMAT  = KT2AM2 + NCKI(ISYMA)*NUMIA
                  KWMAT  = KVMAT  + MAXIJ*NUMIB*NUMIA
                  KQMAT  = KWMAT  + MAXIJ*NUMIB*NUMIA
                  KRMAT  = KQMAT  + MAXIJ*NUMIB*NUMIA
                  KENDA1 = KRMAT  + MAXIJ*NUMIB*NUMIA
                  LWRKA1 = LWORK  - KENDA1
                  IF (LWRKA1 .LT. 0)
     &            CALL QUIT('Not enough space in CC_CHOPT')
C
C--------------------------------
C                 Read integrals.
C--------------------------------
C
C                 We read again from  CHO_VI1 and CHO_VI2 as this is CCSD(T)
C                 ----------------------------------------------------------
C
                  ICOUNT = 0
                  DO ISYMC = 1,NSYM
                     ISYMEM = MULD2H(ISYEMC,ISYMC)
                     IOFSET(ISYMC) = ICOUNT
                     ICOUNT = ICOUNT + NT1AM(ISYMEM)*NVIR(ISYMC)*NUMIA
                  END DO
C
C                 (ac | me) = K(em,c,a) ; (ae | mc) = L(em,c,a)
C                 packed according to ISYMC and then as (em,c,a)
C                 ----------------------------------------------
C
                  CALL CCHO_RDKIN(LUVI1,FILVI1,WORK(KKINT),WORK(KLINT),
     &                            WORK(KMXCKA),MXCKA,
     &                            IA1,NUMIA,ISYMA,IOFSET)
C
C-----------------------------------------
C                 Set up IOFF4 for I term.
C-----------------------------------------
C
                  ICOUN4 = 0
                  DO ISYMJ = 1,NSYM
                     ISYMDE = MULD2H(ISYMJ,ISYMA)
                     IOFF4(ISYMJ) = ICOUN4
                     ICOUN4 = ICOUN4 + NMATAB(ISYMDE)*NUMIA*NRHF(ISYMJ)
                  ENDDO
C
C---------------------------------------
C                 Extract T2-amplitudes.
C---------------------------------------
C
C                 TIJAB for B term.
C                 -----------------
C
                  NIJ = NMATIJ(ISYMIJ)
                  CALL CCHO_TTILD1(WORK(KT2VO),WORK(KTIJAB),
     &                 NUMIA,IA1,ISYMA,NUMIB,IB1,
     &                 ISYMB,NIJ,ISYMT2)
C
C                 Amplitudes for I term.
C                 Set up IOFF5.
C                 ----------------------
C
                  CALL CCHO_IXT2(WORK(KT2VO),WORK(KT2AM1),WORK(KT2AM2),
     &                           IA1,NUMIA,ISYMA,IOFF5)
C
C                 TBIJC for A term.
C                 -----------------
C
                  NIJC = NMAIJA(ISYIJC)
                  CALL CCHO_TBIJC(WORK(KT2VO),WORK(KTBIJC),
     &                 NUMIB,IB1,ISYMB,ISYIJC,NIJC)
C
C                 TIJCA for A term.
C                 -----------------
C
                  CALL CCHO_TIJCA(WORK(KT2VO),WORK(KTIJCA),
     &                 NUMIA,IA1,ISYMA,ISYIJC,NIJC)
C
C-----------------------------------------------
C                 Convergence control variables.
C-----------------------------------------------
C
                  BCONV = .FALSE.
                  ACONV = .FALSE.
                  ICONV = .FALSE.
C
C                 Start loop over Cholesky vectors.
C                 =================================
C
                  DO 999 ICHO = 1,NUMCHO
C
                  IF (ACONV .AND. BCONV .AND. ICONV) GOTO 999
C
                  BCOR  = ZERO
                  ACOR  = ZERO
                  CORI  = ZERO
C
C-------------------------------------------
C                 Construct Cholesky vector.
C-------------------------------------------
C
                  CALL CCHO_VEC(FOCKD,CHOELE,WORK(KCHOVE),ICHO,NUMCHO,
     &                          NUMIB,IB1,ISYMB,ISYMA,ISYEMC,
     &                          ISYMIJ,ISYIJC)
C
C------------------------------------------------------
C                 Construct occupied vector for I term.
C------------------------------------------------------
C
                  IF (.NOT. ICONV) THEN
                     KOFF3 = KCHOVE + NMATIJ(ISYMIJ)*NUMIB
     &                     + NCKATR(ISYEMC) + NMAIJA(ISYIJC)
                     CALL CCHO_DECHO3(FOCKD,CHOELE,NUMCHO,ICHO,
     &                                WORK(KOFF3),VDUMMY,
     &                                ISYCKI,ISYDUM,1,0)
                  ENDIF
C
C--------------------------------------------------------------------
C                 Scale integrals and amplitudes by Cholesky vectors.
C--------------------------------------------------------------------
C
                  KOFF3 = KCHOVE + NMATIJ(ISYMIJ)*NUMIB
                  DO ISYMC = 1,NSYM
C
                     KOFF1 = KKINT  + IOFSET(ISYMC)
                     KOFF2 = KLINT  + IOFSET(ISYMC)
C
                     ISYMEM = MULD2H(ISYMC,ISYEMC)
C
                     NEMC = NT1AM(ISYMEM)*NVIR(ISYMC) 
C
                     CALL CCHO_SCLINT(WORK(KOFF1),WORK(KOFF2),
     &                             WORK(KOFF3),NUMIA,NEMC)
C
                     KOFF3 = KOFF3 + NEMC
C
                  ENDDO
C
                  CALL CCHO_TTILD2(WORK(KTIJAB),WORK(KCHOVE),
     &                             NUMIA,NUMIB,NIJ)
C
                  IF (.NOT. ICONV) THEN
C
C--------------------------------------------------
C                    Resort scaled KINT for I term.
C--------------------------------------------------
C
                     CALL CCHO_IXKTIL(WORK(KKINT),WORK(KKTILD),
     &                                IOFSET,IOFF4,NUMIA,ISYMA)
C
C------------------------------------------------
C                    Scale amplitudes for I term.
C                    (Update vectors.)
C------------------------------------------------
C
                     KOFF3 = KCHOVE + NMATIJ(ISYMIJ)*NUMIB
     &                     + NCKATR(ISYEMC) + NMAIJA(ISYIJC)
                     CALL CCHO_ISCLT(WORK(KT2AM1),WORK(KT2AM2),
     &                               IOFF5,WORK(KOFF3),NUMIA,ISYMA)

c      ocnorm = dsqrt(ddot(NCKI(ISYIJC),WORK(KOFF3),1,WORK(KOFF3),1))
c      write(LUPRI,*) '   After ISCLT:'
c      write(LUPRI,*) '          ICHO,ISYIJC           : ',ICHO,ISYIJC
c      write(LUPRI,*) '          Norm of d(ki,c) update: ',ocnorm
C
                  ENDIF
C
C============================
C                 The B term.
C============================
C
                  IF (.NOT. BCONV) THEN
C
                  TB = SECOND()

C Decommented by Domenico                 
                 TNOW = SECOND()
                 DELTAT = TNOW - TLAST
                 TLAST = TNOW
                 SCNDSX(ICHO) = SCNDSX(ICHO)
     &                         + DELTAT
C Decommented by Domenico
C
C-----------------------------------------------------
C                 Loop over b to construct the B term.
C-----------------------------------------------------
C
                  DO B = 1,NUMIB
C
                     IB = B + IB1 - 1
C
C--------------------------------
C                    Construct X.
C--------------------------------
C
                     NCD = NMATAB(ISYMCD)
C
                     KOFF1 = KT2VO  + IT2VO(ISYMCD,ISYMIJ)
                     KOFF2 = KTIJAB + NMATIJ(ISYMIJ)*NUMIA*(B - 1)
C
                     T1 = SECOND()
                     CALL DGEMM('N','N',NCD,NUMIA,NIJ,ONE,WORK(KOFF1),
     &                   NCD,WORK(KOFF2),NIJ,ZERO,WORK(KXMAT),NCD)
                     T1 = SECOND() - T1
                     TIMEX = TIMEX + T1
C
                     DO 300 ISYMD = 1,NSYM
C
                        ISYMEM = MULD2H(ISYEMD,ISYMD)
                        ISYMC  = MULD2H(ISYMD,ISYMCD)
C
                        JSYEMC = MULD2H(ISYMEM,ISYMC)
                        IF (JSYEMC .NE. ISYEMC)
     &                  CALL QUIT('That is fishy')
C
                        NC   = NVIR(ISYMC)
                        ND   = NVIR(ISYMD)
                        NEM  = NT1AM(ISYMEM)
                        NCA  = NC*NUMIA
                        NEMD = NEM*ND
C
                        IF (NEM .EQ. 0) GOTO 300
                        IF (NC  .EQ. 0) GOTO 300
                        IF (ND  .EQ. 0) GOTO 300
C
C-----------------------------------
C                       Construct Y.
C-----------------------------------
C
                        CALL CCHO_Y1(WORK(KXMAT),WORK(KY1MAT),
     &                               WORK(KY2MAT),NUMIA,NCD,
     &                               NC,ISYMC,ND,ISYMD)
C
C------------------------------------------------------
C                       Construct Z used in first term.
C------------------------------------------------------
C
                        T1 = SECOND()
                        KOFF = KKINT + IOFSET(ISYMC)
                        CALL DGEMM('N','N',NEM,ND,NCA,ONE,
     &                             WORK(KOFF),NEM,WORK(KY1MAT),NCA,
     &                             ZERO,WORK(KZMAT),NEM)
                        T1 = SECOND() - T1
                        TIMEZ = TIMEZ + T1
C
C--------------------------------------------------
C                       Add contribution to energy.
C--------------------------------------------------
C
                        KOFF1 = KJINT + NCKATR(ISYEMD)*(B-1)
     &                                + ICKATR(ISYMEM,ISYMD)
                        KOFF2 = KIINT + NCKATR(ISYEMD)*(B-1)
     &                                + ICKATR(ISYMEM,ISYMD)
                        BCOR = BCOR 
     &                   - TWO*DDOT(NEMD,WORK(KZMAT),1,WORK(KOFF1),1)
     &                   +     DDOT(NEMD,WORK(KZMAT),1,WORK(KOFF2),1)
C
C-------------------------------------------------------
C                       Construct Z used in second term.
C-------------------------------------------------------
C
                        KOFF = KLINT + IOFSET(ISYMC)
                        CALL DGEMM('N','N',NEM,ND,NCA,ONE,
     &                             WORK(KOFF),NEM,WORK(KY2MAT),NCA,
     &                             ZERO,WORK(KZMAT),NEM)
C
C--------------------------------------------------
C                       Add contribution to energy.
C--------------------------------------------------
C
                        KOFF2 = KIINT + NCKATR(ISYEMD)*(B-1)
     &                                + ICKATR(ISYMEM,ISYMD)
                        BCOR = BCOR
     &                       + DDOT(NEMD,WORK(KZMAT),1,WORK(KOFF2),1)
C
C------------------------------------------------------
C                       Construct Z used in third term.
C------------------------------------------------------
C
                        KOFF = KLINT + IOFSET(ISYMC)
                        CALL DGEMM('N','N',NEM,ND,NCA,ONE,
     &                             WORK(KOFF),NEM,WORK(KY1MAT),NCA,
     &                             ZERO,WORK(KZMAT),NEM)
C
C--------------------------------------------------
C                       Add contribution to energy.
C--------------------------------------------------
C
                        KOFF1 = KJINT + NCKATR(ISYEMD)*(B-1)
     &                                + ICKATR(ISYMEM,ISYMD)
                        BCOR = BCOR
     &                       + DDOT(NEMD,WORK(KZMAT),1,WORK(KOFF1),1)
C
  300                CONTINUE
C
                  END DO
C
                  E4BT = E4BT + BCOR
C 
                 ENERGB(ICHO) = ENERGB(ICHO) + BCOR
C
C De-commented by Domenico
                 TNOW = SECOND()
                 DELTAT = TNOW - TLAST
                 TLAST = TNOW
                 SCNDSB(ICHO) = SCNDSB(ICHO) 
     &                         + DELTAT
 
c                 WRITE(LUPRI,'(A,5I4,2D15.6)')
c    &                 'ISB,IBB,ISA,IBA,IV,BCO,BT',
c    &                      ISYMB,IBATB,ISYMA,IBATA,ICHO,BCOR,E4BT
C
                  IF (ABS(BCOR) .LT. THRCHO) BCONV = .TRUE.
C
                  TB = SECOND() - TB
                  TIMEB = TIMEB + TB
C
                  END IF
C
C============================
C                 The A term.
C============================
C
                  IF (.NOT. ACONV .AND. (ISYMA .EQ. ISYMB)) THEN
C
                  KOFF = KCHOVE + NMATIJ(ISYMIJ)*NUMIB
     &                 + NCKATR(ISYEMC)
                  T1 = SECOND()
                  CALL CCHO_AT2SC(NIJC,NUMIA,WORK(KOFF),WORK(KTIJCA))
                  T1 = SECOND() - T1
                  TIMA1 = TIMA1 + T1
C
                  T1 = SECOND()
                  CALL CCHO_ATERM(WORK(KAMAT),WORK(KJINT),WORK(KKINT),
     &                            WORK(KLINT),WORK(KMXCKA),WORK(KXAMA),
     &                            WORK(KTBIJC),WORK(KTIJCA),NUMIA,
     &                            NUMIB,NIJC,ISYMA,ISYMB,ISYIJC,ACOR)
                  T1 = SECOND() - T1
                  TIMA2 = TIMA2 + T1
                  CALL FLSHFO(LUPRI)
C
                  E4AT = E4AT + ACOR
C 
                 ENERGA(ICHO) = ENERGA(ICHO) + ACOR
C
C 
C Decommented by Domenico
                 TNOW = SECOND()
                 DELTAT = TNOW - TLAST
                 TLAST = TNOW
                 SCNDSA(ICHO) = SCNDSA(ICHO) 
     &                         + DELTAT
 
C Decommented by Domenico
c                 WRITE(LUPRI,'(A,5I4,2D15.6)')
c    &            'ISB,IBB,ISA,IBA,IV,ACO,AT',
c    &                      ISYMB,IBATB,ISYMA,IBATA,ICHO,ACOR,E4AT
C
                  IF (ABS(ACOR) .LT. THRCHO) ACONV = .TRUE.
C
                  END IF
C
C============================
C                 The I term.
C============================
C
                  IF (.NOT. ICONV) THEN
c     write(LUPRI,*) '   Before I:'
c     write(LUPRI,'(4X,A,/,4X,I1,1X,I4,1X,I1,1X,I4,1X,I4,1X,2D15.6)')
c    &      'ISYMB,IBATB,ISYMA,IBATA,ICHO,CORI,E4I: ',
c    &       ISYMB,IBATB,ISYMA,IBATA,ICHO,CORI,E4I

c      write(LUPRI,*) '   Space for v,w,q,r: ',
c     &           MAXIJ*NUMIA*NUMIB
C
                     T1 = SECOND()
                     CALL CCHO_ITERM(WORK(KJTIL1),WORK(KJTIL2),
     &                               WORK(KT2AM3),WORK(KT2AM4),
     &                               WORK(KKTILD),WORK(KT2AM1),
     &                               WORK(KT2AM2),WORK(KVMAT),
     &                               WORK(KWMAT),WORK(KQMAT),
     &                               WORK(KRMAT),CORI,NUMIA,NUMIB,
     &                               ISYMA,ISYMB,IOFF2,IOFF3,IOFF4,
     &                               IOFF5)
                     T1 = SECOND() - T1
                     TIMEI = TIMEI + T1
                     CALL FLSHFO(LUPRI)
C
                     E4I = E4I + CORI
C 
                    ENERGI(ICHO) = ENERGI(ICHO) + CORI
C
C Decommented by Domenico
                    TNOW = SECOND()
                    DELTAT = TNOW - TLAST
                    TLAST = TNOW
                    SCNDSI(ICHO) = SCNDSI(ICHO) 
     &                            + DELTAT
C Decommented by Domenico
c     write(LUPRI,*) '   After  I:'
c     write(LUPRI,'(4X,A,/,4X,I1,1X,I4,1X,I1,1X,I4,1X,I4,1X,2D15.6)')
c    &      'ISYMB,IBATB,ISYMA,IBATA,ICHO,CORI,E4I: ',
c    &       ISYMB,IBATB,ISYMA,IBATA,ICHO,CORI,E4I
C
                      IF (ABS(CORI) .LT. THRCHO) ICONV = .TRUE.
C
                  ENDIF
C
c                 IF (IPRINT .GE. INFO) THEN
c                    WRITE(LUPRI,'(15X,A,I3,A,/,15X,A)')
c    &               'Status after Cholesky vector',ICHO,':',
c    &               '--------------------------------'
c                    IF (ISYMA .EQ. ISYMB) THEN
c                       IF (ACONV) THEN
c                          WRITE(LUPRI,'(15X,A)')
c    &                     'A term: converged'
c                       ELSE
c                          WRITE(LUPRI,'(15X,A)')
c    &                     'A term: not converged'
c                       ENDIF
c                    ENDIF
c                    IF (BCONV) THEN
c                       WRITE(LUPRI,'(15X,A)')
c    &                  'B term: converged'
c                    ELSE
c                       WRITE(LUPRI,'(15X,A)')
c    &                  'B term: not converged'
c                    ENDIF
c                    IF (ICONV) THEN
c                       WRITE(LUPRI,'(15X,A)')
c    &                  'I term: converged'
c                    ELSE
c                       WRITE(LUPRI,'(15X,A)')
c    &                  'I term: not converged'
c                    ENDIF
c                    IF (ISYMA .EQ. ISYMB) THEN
c                       WRITE(LUPRI,'(15X,A,F10.2,A)')
c    &                  'Accumulated A-time: ',TIMA1+TIMA2,' seconds'
c                    ENDIF
c                    WRITE(LUPRI,'(15X,A,F10.2,A)')
c    &               'Accumulated B-time: ',TIMEB,' seconds'
c                    WRITE(LUPRI,'(15X,A,F10.2,A,/)')
c    &               'Accumulated I-time: ',TIMEI,' seconds'
c                 ENDIF
C
  999             CONTINUE
C
               ENDDO
C
  200       CONTINUE
C
            IF (.TRUE.) THEN
C
C===================================
C              Calculate the D term.
C===================================
C
               T1 = SECOND()
               CALL CCHO_DTERM(WORK(KJINT),WORK(KOINT),WORK(KT2VO),
     &                         FOCKD,
     &                         NUMCHO,CHOELE,WORK(KENDB2),LWRKB2,E4D,
     &                         IB1,ISYMB,NUMIB,
     &                         FBATCH,IPRINT.GE.INFO,NONA)
               T1 = SECOND() - T1
               TIMED = TIMED + T1
C
C====================================
C              Calculate the E terms.
C====================================
C
c     write(LUPRI,*) 'Calling ETERM'
               T1 = SECOND()
               CALL CCHO_ETERM(WORK(KIINT),WORK(KJINT),WORK(KT2VO),
     &                         FOCKD,
     &                         NUMCHO,CHOELE,WORK(KENDB2),LWRKB2,E4E1,
     &                         IB1,ISYMB,NUMIB,
     &                         FBATCH,IPRINT.GE.INFO,NONA)
               T1 = SECOND() - T1
               TIME1 = TIME1 + T1
c     write(LUPRI,*) 'End of ETERM, E4E1 = ',E4E1
C
c     write(LUPRI,*) 'Calling ETERM2'
               T1 = SECOND()
               CALL CCHO_ETERM2(WORK(KIINT),WORK(KJINT),WORK(KOINT),
     &                          WORK(KT2VO),FOCKD,
     &                          NUMCHO,CHOELE,WORK(KENDB2),LWRKB2,E4E2,
     &                          IB1,ISYMB,NUMIB,
     &                          FBATCH,IPRINT.GE.INFO,NONA)
               T1 = SECOND() - T1
               TIME2 = TIME2 + T1
c     write(LUPRI,*) 'End of ETERM2, E4E2 = ',E4E2
C
C====================================
C              Calculate the J terms.
C====================================
C
               T1 = SECOND()
               CALL CCHO_JTERM(WORK(KIINT),WORK(KJINT),WORK(KOINT),
     &                         WORK(KT2VO),FOCKD,
     &                         NUMCHO,CHOELE,WORK(KENDB2),LWRKB2,E4J,
     &                         IB1,ISYMB,NUMIB,
     &                         FBATCH,IPRINT.GE.INFO,NONA)
               T1 = SECOND() - T1
               TIMEJ = TIMEJ + T1
C
C====================================
C              Calculate the G1 term.
C====================================
C
               T1 = SECOND()
               CALL CCHO_G1TERM(WORK(KIINT),WORK(KJINT),WORK(KT2VO),
     &                          T1AM,WORK(KIAJB),FOCKD,
     &                          NUMCHO,CHOELE,WORK(KENDB2),LWRKB2,E5G1,
     &                          IB1,ISYMB,NUMIB,
     &                          FBATCH,IPRINT.GE.INFO,NONA)
               T1 = SECOND() - T1
               TIMG1 = TIMG1 + T1
C
C====================================
C              Calculate the G terms.
C====================================
C
               T1 = SECOND()
               CALL CCHO_GTERM(WORK(KIINT),WORK(KJINT),WORK(KT2VO),
     &                         T1AM,WORK(KIAJB),FOCKD,
     &                         NUMCHO,CHOELE,WORK(KENDB2),LWRKB2,E5G,
     &                         IB1,ISYMB,NUMIB,
     &                         FBATCH,IPRINT.GE.INFO,NONA)
c    &                         EYS,EXR,EWQ,EVP)
               T1 = SECOND() - T1
               TIMEG = TIMEG + T1
C
C------------------------------------
C              Calculate the F2 term.
C------------------------------------
C
               T1 = SECOND()
               CALL CCHO_FTERM2(WORK(KIINT),WORK(KJINT),
     &                          WORK(KOINT),WORK(KT2VO),FOCKD,NUMCHO,
     &                          CHOELE,WORK(KENDB2),LWRKB2,
     &                          E4F2,IB1,ISYMB,NUMIB,
     &                          IPRINT.GE.INFO)
               T1 = SECOND() - T1
               TIMF2 = TIMF2 + T1
C
               IF (IPRINT .GE. INFO) THEN
                  WRITE(LUPRI,'(3X,A,F10.2,A)')
     &            'Accumulated A-time     : ',TIMA1+TIMA2,' seconds'
                  WRITE(LUPRI,'(3X,A,F10.2,A)')
     &            'Accumulated B-time     : ',TIMEB,' seconds'
                  WRITE(LUPRI,'(3X,A,F10.2,A)')
     &            'Accumulated D-time     : ',TIMED,' seconds'
                  WRITE(LUPRI,'(3X,A,F10.2,A)')
     &            'Accumulated E1-time    : ',TIME1,' seconds'
                  WRITE(LUPRI,'(3X,A,F10.2,A)')
     &            'Accumulated E2-time    : ',TIME2,' seconds'
                  WRITE(LUPRI,'(3X,A,F10.2,A)')
     &            'Accumulated F2-time    : ',TIMF2,' seconds'
                  WRITE(LUPRI,'(3X,A,F10.2,A)')
     &            'Accumulated G1-time    : ',TIMG1,' seconds'
                  WRITE(LUPRI,'(3X,A,F10.2,A)')
     &            'Accumulated G-time     : ',TIMEG,' seconds'
                  WRITE(LUPRI,'(3X,A,F10.2,A)')
     &            'Accumulated I-time     : ',TIMEI,' seconds'
                  WRITE(LUPRI,'(3X,A,F10.2,A)')
     &            'Accumulated J-time     : ',TIMEJ,' seconds'
                  TTOT = SECOND() - TITOT
                  WRITE(LUPRI,'(3X,A,F10.2,A,/)')
     &            'Total time used so far : ',TTOT,' seconds'
                  CALL FLSHFO(LUPRI)
               ENDIF
C
            ENDIF
           TLAST = SECOND()
C
C========================================================
C          Write restart information at the end of batch.
C========================================================
C
         C4ORV = E4AT + E4BT + E4D + E4E1 + E4E2 + E4F2 + E4I + E4J
         C4ORV = 2.0D0 * C4ORV + OLD4V
         C5ORV = E5G1 + E5G
         C5ORV = 2.0D0 * C5ORV + OLD5V
C
         WRITE(LUCHOU,'(A,I2,2I4)') 'B-symmetry, first and last B :',
     &                               ISYMB, IB1, IB2
         WRITE(LUCHOU,'(A,D20.10)') 'A  contribution :', E4AT
         WRITE(LUCHOU,'(A,D20.10)') 'B  contribution :', E4BT
         WRITE(LUCHOU,'(A,D20.10)') 'D  contribution :', E4D 
         WRITE(LUCHOU,'(A,D20.10)') 'E1 contribution :', E4E1
         WRITE(LUCHOU,'(A,D20.10)') 'E2 contribution :', E4E2
         WRITE(LUCHOU,'(A,D20.10)') 'F2 contribution :', E4F2
         WRITE(LUCHOU,'(A,D20.10)') 'G1 contribution :', E5G1
         WRITE(LUCHOU,'(A,D20.10)') 'G  contribution :', E5G 
         WRITE(LUCHOU,'(A,D20.10)') 'I  contribution :', E4I 
         WRITE(LUCHOU,'(A,D20.10)') 'J  contribution :', E4J 
         WRITE(LUCHOU,'(A,20A)') '-----------------',('-',IK = 1,20)
C
         WRITE(LUCHOU,'(A,D20.10)') '4-V correction  :', C4ORV
         WRITE(LUCHOU,'(A,D20.10)') '5-V correction  :', C5ORV
         WRITE(LUCHOU,*)
         WRITE(LUCHOU,*)
         CALL FLSHFO(LUCHOU)
C
         END DO
C
         WRITE(LUCHOU,'(9X,A,I2,//)') 
     &       'Completed contribution from virtuals of symmetry :',ISYMB
C
         IF (RSTVIR) THEN
            RSTVIR = .FALSE.                  ! No restart in next symmetry
            IFVIOR = 1
         END IF
C
  100 CONTINUE
C
      WRITE(LUCHOU,*)
      WRITE(LUCHOU,*)
      WRITE(LUCHOU,'(9X,A)') '======================================'
      WRITE(LUCHOU,'(9X,A)') 'Virtual contribution to CHOPT finished'
      WRITE(LUCHOU,'(9X,A)') '======================================'
      WRITE(LUCHOU,*)
      WRITE(LUCHOU,*)
      CALL FLSHFO(LUCHOU)
C
      IF (IPRINT .GE. INFO) THEN
         TIMEV = SECOND() - TIMEV
         WRITE(LUPRI,'(A,F10.2,A,/,A,/)')
     &   'Contributions from virtual integrals completed in ',TIMEV,
     &   ' seconds',
     &   'Starting contributions from occupied integrals'
         CALL FLSHFO(LUPRI)
      ENDIF
C
C===========================
C     Calculate the H terms.
C===========================
C
  222 CONTINUE           !H term
C
      WRITE(LUCHOU,'(A)') 'About to calculate H-term'
      WRITE(LUCHOU,'(A)') '.........'
      CALL FLSHFO(LUCHOU)
C 
      TIMEH = SECOND()
      CALL CCHO_HTERM(FOCKD,T1AM,WORK(KOINT),WORK(KT2VO),WORK(KIAJB),
     &                WORK(KEND1),LWRK1,CHOELE,NUMCHO,E5H,
     &                FBATCH,IPRINT.GE.INFO,NONB,NONA)
      TIMEH = SECOND() - TIMEH
C
      C5ORO = C5ORO + TWO*E5H
      WRITE(LUCHOU,'(A,D20.10)') 'H  contribution :', E5H 
      WRITE(LUCHOU,'(A,D20.10)') '5-O correction  :', C5ORO
      WRITE(LUCHOU,*)
      CALL FLSHFO(LUCHOU)
C
  223 CONTINUE           !H1 term
C
      WRITE(LUCHOU,'(A)') 'About to calculate H1-term'
      WRITE(LUCHOU,'(A)') '.........'
      CALL FLSHFO(LUCHOU)
C 
      TIMH1 = SECOND()
      CALL CCHO_H1TERM(FOCKD,T1AM,WORK(KOINT),WORK(KT2VO),WORK(KIAJB),
     &                 WORK(KEND1),LWRK1,CHOELE,NUMCHO,E5H1,
     &                 FBATCH,IPRINT.GE.INFO,NONB,NONA)
      TIMH1 = SECOND() - TIMH1
C
      C5ORO = C5ORO + TWO*E5H1
      WRITE(LUCHOU,'(A,D20.10)') 'H1 contribution :', E5H1
      WRITE(LUCHOU,'(A,D20.10)') '5-O correction  :', C5ORO
      WRITE(LUCHOU,*)
      CALL FLSHFO(LUCHOU)
C
      IF (IPRINT .GE. INFO) THEN
         WRITE(LUPRI,'(A,F10.2,A)')
     &   'Time for H-terms: ',TIMEH+TIMH1,' seconds'
         CALL FLSHFO(LUPRI)
      ENDIF
C
C
C=================================
C           Calculate the F1 term.
C=================================
C
  224 CONTINUE           !F1 term
C
      KENDB2 = KEND1
      LWRKB2 = LWRK1
C
      WRITE(LUCHOU,'(A)') 'About to calculate F1-term'
      WRITE(LUCHOU,'(A)') '.........'
      CALL FLSHFO(LUCHOU)
C 
      TIMF1 = SECOND()
      CALL CCHO_FTERM1(WORK(KOINT),WORK(KT2VO),FOCKD,NUMCHO,
     &                 CHOELE,WORK(KENDB2),LWRKB2,E4F1,
     &                 FBATCH,IPRINT.GE.INFO,NONB,NONA)
      TIMF1 = SECOND() - TIMF1
C
      C4ORO = C4ORO + TWO*E4F1
      WRITE(LUCHOU,'(A,D20.10)') 'F1 contribution :', E4F1
      WRITE(LUCHOU,'(A,D20.10)') '4-O correction  :', C4ORO
      WRITE(LUCHOU,*)
      CALL FLSHFO(LUCHOU)
C
      IF (IPRINT .GE. INFO) THEN
         WRITE(LUPRI,'(A,F10.2,A)')
     &   'Time for F1-term: ',TIMF1,' seconds'
         CALL FLSHFO(LUPRI)
      ENDIF
C
C=================================
C           Calculate the C1 term.
C=================================
C
  225 CONTINUE           !C1 term
C
      KENDB2 = KEND1
      LWRKB2 = LWRK1
C
      WRITE(LUCHOU,'(A)') 'About to calculate C1-term'
      WRITE(LUCHOU,'(A)') '.........'
      CALL FLSHFO(LUCHOU)
C 
      TIMC1 = SECOND()
      CALL CC_CHOC1TERM(WORK(KOINT),WORK(KT2VO),FOCKD,NUMCHO,
     &                 CHOELE,WORK(KENDB2),LWRKB2,E4C1,
     &                 FBATCH,IPRINT.GE.INFO,NONB,NONA)
      TIMC1 = SECOND() - TIMC1
C
      C4ORO = C4ORO + TWO*E4C1
      WRITE(LUCHOU,'(A,D20.10)') 'C1 contribution :', E4C1
      WRITE(LUCHOU,'(A,D20.10)') '4-O correction  :', C4ORO
      WRITE(LUCHOU,*)
      CALL FLSHFO(LUCHOU)
C
*     CALL C1_DEBUG(WORK(KOINT),WORK(KT2VO),FOCKD,NUMCHO,
*    &                 CHOELE,WORK(KENDB2),LWRKB2,THRCHO,
*    &                 NRHFT,NVIRT)
C
C=================================
C           Calculate the C2 term.
C=================================
C
  226 CONTINUE           !C2 term
C
      KENDB2 = KEND1
      LWRKB2 = LWRK1
C
      WRITE(LUCHOU,'(A)') 'About to calculate C2-term'
      WRITE(LUCHOU,'(A)') '.........'
      CALL FLSHFO(LUCHOU)
C 
      TIMC2 = SECOND()
      CALL CC_CHOC2TERM(WORK(KOINT),WORK(KT2VO),FOCKD,NUMCHO,
     &                 CHOELE,WORK(KENDB2),LWRKB2,E4C2,
     &                 FBATCH,IPRINT.GE.INFO,NONB,NONA)
      TIMC2 = SECOND() - TIMC2
C
      C4ORO = C4ORO + TWO*E4C2
      WRITE(LUCHOU,'(A,D20.10)') 'C2 contribution :', E4C2
      WRITE(LUCHOU,'(A,D20.10)') '4-O correction  :', C4ORO
      WRITE(LUCHOU,*)
      CALL FLSHFO(LUCHOU)
C
*     CALL C2_DEBUG(WORK(KOINT),WORK(KT2VO),FOCKD,NUMCHO,
*    &                 CHOELE,WORK(KENDB2),LWRKB2,THRCHO,
*    &                 NRHFT,NVIRT)
C
      WRITE(LUCHOU,'(A,D20.10)') 'C1 contribution :', E4C1
      WRITE(LUCHOU,'(A,D20.10)') 'C2 contribution :', E4C2
      WRITE(LUCHOU,'(A,D20.10)') 'F1 contribution :', E4F1
      WRITE(LUCHOU,'(A,D20.10)') 'H1 contribution :', E5H1
      WRITE(LUCHOU,'(A,D20.10)') 'H  contribution :', E5H
      WRITE(LUCHOU,'(A,20A)') '-----------------',('-',IK = 1,20)
C
      WRITE(LUCHOU,'(A,D20.10)') '4-O correction  :', C4ORO
      WRITE(LUCHOU,'(A,D20.10)') '5-O correction  :', C5ORO
      WRITE(LUCHOU,*)
      WRITE(LUCHOU,*)
      CALL FLSHFO(LUCHOU)

C
C
      WRITE(LUCHOU,*)
      WRITE(LUCHOU,*)
      WRITE(LUCHOU,'(9X,A)') '======================================='
      WRITE(LUCHOU,'(9X,A)') 'Occupied contribution to CHOPT finished'
      WRITE(LUCHOU,'(9X,A)') '======================================='
      WRITE(LUCHOU,*)
      WRITE(LUCHOU,*)
C
      WRITE(LUCHOU,'(A,D20.10)') 'E(4) correction :', C4ORO + C4ORV
      WRITE(LUCHOU,'(A,D20.10)') 'E(5) correction :', C5ORO + C5ORV
      WRITE(LUCHOU,'(A,20A)') '-----------------',('-',IK = 1,20)
      WRITE(LUCHOU,'(A,D20.10)') 'E(T) correction :', 
     &                           C4ORO + C4ORV + C5ORO + C5ORV 
      WRITE(LUCHOU,*)
      WRITE(LUCHOU,*)
      CALL FLSHFO(LUCHOU)
C
      IF (IPRINT .GE. INFO) THEN
         WRITE(LUPRI,'(A,F10.2,A)')
     &   'Time for C-terms: ',TIMC1+TIMC2,' seconds'
         CALL FLSHFO(LUPRI)
      ENDIF
C
C===============================================
C     Calculate total CCSD(T) energy correction.
C===============================================
C
      CALL HEADER('CC_CHOPT: CCSD(T) Energy Correction',-1)
C
c     write(LUPRI,*)
c     write(LUPRI,*) '   From GTERM:'
c     write(LUPRI,*) '   EYS = ',EYS
c     write(LUPRI,*) '   EXR = ',EXR
c     write(LUPRI,*) '   EWQ = ',EWQ
c     write(LUPRI,*) '   EVP = ',EVP
c     write(LUPRI,*)
C
      E4TH = E4AT + E4BT + E4C1 + E4C2 + E4D + E4E1 + E4E2 + E4F1 + E4F2
     &     + E4I + E4J
      E5TH = E5G1 + E5G + E5H + E5H1
      ETOT = E4TH + E5TH
C
      C4TOT = C4ORV + C4ORO
      C5TOT = C5ORV + C5ORO
C
      XENA = TWO*E4AT
      XENB = TWO*E4BT
      XENC = TWO*(E4C1 + E4C2)
      XEND = TWO*E4D
      XENE = TWO*(E4E1 + E4E2)
      XENF = TWO*(E4F1 + E4F2)
      XENG = TWO*(E5G1 + E5G)
      XENH = TWO*(E5H1 + E5H)
      XENI = TWO*E4I
      XENJ = TWO*E4J
      XEN4 = TWO*E4TH
      XEN5 = TWO*E5TH
C
      WRITE(LUPRI,*) 'ENERGIES/2 FOR EACH CHOLESKY VECTOR:'
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) 'A TERM: '
      WRITE(LUPRI,'(4D20.10)') (ENERGA(I),I=1,MXCHVE)
      WRITE(LUPRI,*) 'B TERM: ' 
      WRITE(LUPRI,'(4D20.10)') (ENERGB(I),I=1,MXCHVE)
      WRITE(LUPRI,*) 'C TERM: ' 
      WRITE(LUPRI,'(4D20.10)') (ENERGC(I),I=1,MXCHVE)
      WRITE(LUPRI,*) 'D TERM: ' 
      WRITE(LUPRI,'(4D20.10)') (ENERGD(I),I=1,MXCHVE)
      WRITE(LUPRI,*) 'E TERM: ' 
      WRITE(LUPRI,'(4D20.10)') (ENERGE(I),I=1,MXCHVE)
      WRITE(LUPRI,*) 'F TERM: ' 
      WRITE(LUPRI,'(4D20.10)') (ENERGF(I),I=1,MXCHVE)
      WRITE(LUPRI,*) 'G TERM: ' 
      WRITE(LUPRI,'(4D20.10)') (ENERGG(I),I=1,MXCHVE)
      WRITE(LUPRI,*) 'H TERM: ' 
      WRITE(LUPRI,'(4D20.10)') (ENERGH(I),I=1,MXCHVE)
      WRITE(LUPRI,*) 'I TERM: ' 
      WRITE(LUPRI,'(4D20.10)') (ENERGI(I),I=1,MXCHVE)
      WRITE(LUPRI,*) 'J TERM: ' 
      WRITE(LUPRI,'(4D20.10)') (ENERGJ(I),I=1,MXCHVE)
      DO I = 1,MXCHVE
         ENERGX(I) = ENERGX(I)
     &             +(ENERGA(I)
     &             + ENERGB(I)
     &             + ENERGC(I)
     &             + ENERGD(I)
     &             + ENERGE(I)
     &             + ENERGF(I)
     &             + ENERGG(I)
     &             + ENERGH(I)
     &             + ENERGI(I)
     &             + ENERGJ(I))*TWO
      ENDDO
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) 'Total contribiution of each vector ' 
      WRITE(LUPRI,*)

C-tbp WRITE(LUPRI,'(4D20.10)') (ENERGX(I),I=1,MXCHVE)
      DO I = 1,MXCHVE
         IF (I .LT. 10) THEN
            WRITE(LUPRI,'(A,I1,A,D20.10)')
     &      'ECHOPT_',I,':',ENERGX(I)
         ELSE IF (I .LT. 100) THEN
            WRITE(LUPRI,'(A,I2,A,D20.10)')
     &      'ECHOPT_',I,':',ENERGX(I)
         ELSE IF (I .LT. 1000) THEN
            WRITE(LUPRI,'(A,I3,A,D20.10)')
     &      'ECHOPT_',I,':',ENERGX(I)
         ELSE IF (I .LT. 10000) THEN
            WRITE(LUPRI,'(A,I4,A,D20.10)')
     &      'ECHOPT_',I,':',ENERGX(I)
         ELSE
            WRITE(LUPRI,*)
     &      'ECHOPT_',I,':',ENERGX(I)
         END IF
      END DO

     
      SUMENERG = ZERO
      DO I = 1,MXCHVE
         SUMENERG = SUMENERG
     &            + ENERGA(I)
     &            + ENERGB(I)
     &            + ENERGC(I)
     &            + ENERGD(I)
     &            + ENERGE(I)
     &            + ENERGF(I)
     &            + ENERGG(I)
     &            + ENERGH(I)
     &            + ENERGI(I)
     &            + ENERGJ(I)
      ENDDO
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(A,D20.10,/)') '--- TOTAL ENERGY: ',TWO*SUMENERG
      WRITE(LUPRI,*)
C
      IF (RSTH.OR.RSTH1.OR.RSTF1.OR.RSTC1.OR.RSTC2.OR.RSTVIR) THEN
         WRITE(LUPRI,'(/,10X,A,/,A,/)') 
     &               'Final results from CHOPT (restart included)',
     &               '==========================================='
         WRITE(LUPRI,'(10X,A,D20.10)') '4th order correction', C4TOT
         WRITE(LUPRI,'(10X,A,D20.10)') '5th order correction', C5TOT
         WRITE(LUPRI,'(10X,A,D20.10)') 'Total (T) correction', 
     &                                                 C4TOT + C5TOT
      END IF
C      
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: A    ',E4AT*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: B    ',E4BT*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: C1   ',E4C1*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: C2   ',E4C2*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: D    ',E4D*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: E1   ',E4E1*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: E2   ',E4E2*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: F1   ',E4F1*2.0d0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: F2   ',E4F2*2.0d0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: J    ',E4J*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: I    ',E4I*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: G    ',E5G*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: G1   ',E5G1*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: H    ',E5H*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: H1   ',E5H1*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: E4th ',E4TH*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: E5th ',E5TH*2.0D0
         WRITE(LUPRI,'(10X,A,D20.10)') 'Contribution: ETOT ',ETOT*2.0D0
      END IF
C
      XENA = E4AT
C
      IF (IPRINT .GT. 6) THEN
         WRITE(LUPRI,'(/,10X,A,F10.2)') 'TIMEA: ',TIMA1 + TIMA2
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMEB: ',TIMEB
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMC1: ',TIMC1
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMC2: ',TIMC2
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMED: ',TIMED
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIME1: ',TIME1
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIME2: ',TIME2
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMF1: ',TIMF1
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMF2: ',TIMF2
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMEJ: ',TIMEJ
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMEI: ',TIMEI
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMEG: ',TIMEG
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMG1: ',TIMG1
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMEH: ',TIMEH
         WRITE(LUPRI,'(10X,A,F10.2)')   'TIMH1: ',TIMH1
         WRITE(LUPRI,'(10X,A,F10.2,/)') 'TITOT: ',SECOND() - TITOT
      END IF
c    
C De-commented by Domenico
      IF (IPRINT .GT. 10) THEN
         WRITE(LUPRI,*) 'TIMINGS FOR EACH CHOLESKY VECTOR:'
         WRITE(LUPRI,*) 'SET UP: ' 
         WRITE(LUPRI,'(4F15.1)') (SCNDSX(I),I=1,MXCHVE)
         WRITE(LUPRI,*) 'A TERM: '
         WRITE(LUPRI,'(4F15.1)') (SCNDSA(I),I=1,MXCHVE)
         WRITE(LUPRI,*) 'B TERM: ' 
         WRITE(LUPRI,'(4F15.1)') (SCNDSB(I),I=1,MXCHVE)
         WRITE(LUPRI,*) 'C TERM: ' 
         WRITE(LUPRI,'(4F15.1)') (SCNDSC(I),I=1,MXCHVE)
         WRITE(LUPRI,*) 'D TERM: ' 
         WRITE(LUPRI,'(4F15.1)') (SCNDSD(I),I=1,MXCHVE)
         WRITE(LUPRI,*) 'E TERM: ' 
         WRITE(LUPRI,'(4F15.1)') (SCNDSE(I),I=1,MXCHVE)
         WRITE(LUPRI,*) 'F TERM: ' 
         WRITE(LUPRI,'(4F15.1)') (SCNDSF(I),I=1,MXCHVE)
         WRITE(LUPRI,*) 'G TERM: ' 
         WRITE(LUPRI,'(4F15.1)') (SCNDSG(I),I=1,MXCHVE)
         WRITE(LUPRI,*) 'H TERM: ' 
         WRITE(LUPRI,'(4F15.1)') (SCNDSH(I),I=1,MXCHVE)
         WRITE(LUPRI,*) 'I TERM: ' 
         WRITE(LUPRI,'(4F15.1)') (SCNDSI(I),I=1,MXCHVE)
         WRITE(LUPRI,*) 'J TERM: ' 
         WRITE(LUPRI,'(4F15.1)') (SCNDSJ(I),I=1,MXCHVE)
      END IF
C De-commented by Domenico
C     SUMSCNDS = ZERO
C    DO I = 1,MXCHVE
C       SUMSCNDS = SUMSCNDS
C   &            + SCNDSA(I)
C   &            + SCNDSB(I)
C   &            + SCNDSC(I)
C   &            + SCNDSD(I)
C   &            + SCNDSE(I)
C   &            + SCNDSF(I)
C   &            + SCNDSG(I)
C   &            + SCNDSH(I)
C   &            + SCNDSI(I)
C  &            + SCNDSJ(I)
C   &            + SCNDSX(I)
C    ENDDO
C    WRITE(LUPRI,'(A,F15.1,/)') '--- TOTAL TIME: ',SUMSCNDS
C
C---------------------
C     Close new files.
C---------------------
C
      CALL WCLOSE2(LUVI1,FILVI1,'KEEP')
      CALL WCLOSE2(LUVI2,FILVI2,'KEEP')
C
      CALL GETTIM(TEND,WEND)
      TCPU = TEND -TSTART
      TWAL = WEND -WSTART
      WRITE (LUPRI,'(//A,2F10.2,A/)')
     *      ' CPU and wall time in CC_CHOPT :',TCPU,TWAL,' seconds'
C
      CALL QEXIT('CC_CHOPT')
      RETURN
      END
C  /* Deck ccho_traint */
      SUBROUTINE CCHO_TRAINT(WORK,LWORK,OCINT,LUVI1,FILVI1,LUVI2,FILVI2)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C     - occupied part included, TBP, June 2002.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION WORK(LWORK), OCINT(*)
      CHARACTER*7 FILVI1, FILVI2
C
      CHARACTER*6 FN3VI
      CHARACTER*8 FN3VI2,FNTOC
      PARAMETER (FN3VI  = 'CC3_VI', FN3VI2 = 'CC3_VI12')
C
*     dimension temp(10)
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <ccinftap.h>
#include <cc_cho.h>
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KCMO = 1
      KEND = KCMO  + NLAMDS
      LWRK = LWORK - KEND
C
      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient memory in CCHO_TRAINT')
      ENDIF
C
C----------------------------------------------
C     Read MO-coefficients from interface file.
C----------------------------------------------
C
      LUSIRG = -1
      CALL GPOPEN(LUSIRG,'SIRIFC','OLD',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND LUSIRG
C
      CALL MOLLAB('TRCCINT ',LUSIRG,LUPRI)
      READ (LUSIRG)
C
      READ (LUSIRG)
      READ (LUSIRG) (WORK(KCMO+I-1), I=1,NLAMDS)
C
      CALL GPCLOSE(LUSIRG,'KEEP')
C
C---------------------------------------
C     Reorder the MO-coefficient matrix.
C---------------------------------------
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND),LWRK)
C
C-----------------------------
C     Open old integral files.
C-----------------------------
C
      LU3VI  = 0
      LU3VI2 = 0
      CALL WOPEN2(LU3VI,FN3VI,64,0)
      CALL WOPEN2(LU3VI2,FN3VI2,64,0)
C
C-----------------------------------------------------
C     Transform integrals in CC3_VI (from ccinftap.h).
C-----------------------------------------------------
C
      IF (SKIVI1) THEN
         WRITE(LUPRI,'(//,10X,A,/,10X,A,/)') 
     &         'Using existing CHO_VI1 file', 
     &         '---------------------------'
      ELSE
         WRITE(LUCHOU,'(A)') 'About to build CHO_VI1 file'
         WRITE(LUCHOU,'(A)') '............'
         CALL CCHO_TRINT1(WORK(KCMO),FN3VI,LU3VI,FILVI1,LUVI1,
     &                    WORK(KEND),LWRK)
         WRITE(LUCHOU,'(A,/)') 'CHO_VI1 file built up'
C
         CALL WCLOSE2(LU3VI,FN3VI,'DELETE')
      END IF
C
C-------------------------------------------------------
COLD: Transform integrals in CC3_VI2.
CNEW: Transform integrals in CC3_VI12 (from ccinftap.h).
C-------------------------------------------------------
C
      IF (SKIVI2) THEN
         WRITE(LUPRI,'(/,10X,A,/,10X,A,//)') 
     &         'Using existing CHO_VI2 file', 
     &         '---------------------------'
      ELSE
         WRITE(LUCHOU,'(A)') 'About to build CHO_VI2 file'
         WRITE(LUCHOU,'(A)') '............'
         CALL CCHO_TRINT1(WORK(KCMO),FN3VI2,LU3VI2,FILVI2,LUVI2,
     &                    WORK(KEND),LWRK)
         WRITE(LUCHOU,'(A,//)') 'CHO_VI2 file built up'
C
         CALL WCLOSE2(LU3VI2,FN3VI2,'DELETE')
      END IF
C
C----------------
C     Allocation.
C----------------
C
      KINTOC = KEND
      KEND1  = KINTOC + NTOTOC(1)
      LWRK1  = LWORK  - KEND1 + 1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Insufficient memory in CCHO_TRAINT'
         WRITE(LUPRI,*) 'Need (more than): ',KEND1-1
         WRITE(LUPRI,*) 'Available       : ',LWORK
         CALL QUIT('Insufficient memory in CCHO_TRAINT')
      ENDIF
C
C-------------------------------------------------------
C     Read (ia|j delta) integrals (file should be open).
C-------------------------------------------------------
C
      IOFF   = 1
      IF (NTOTOC(1) .GT. 0) THEN
         LUTOC = 0
         FNTOC = 'CCSDT_OC'
         CALL WOPEN2(LUTOC,FNTOC,64,0)
C        CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(1))
C
Casm     Never read more than 2Gb (268435456 words) in a single shot
C
         LRDMAX = 268435456 / 4 - 1
         NRDBAT = (NTOTOC(1) - 1) / LRDMAX + 1
         LRD1 = 0
         LRD2 = 0
         KOFF = KINTOC
c        write(lupri,*) 'NTOTOC(1),nrdbat',NTOTOC(1),nrdbat
         DO IRDBAT = 1,NRDBAT
            LRD3 = LRDMAX
            LRD1 = LRD1 + LRDMAX
            IF (LRD1 .GT. NTOTOC(1)) LRD3 = NTOTOC(1) - LRD2
c           write(lupri,*) 'ioff,irdbat,lrd1,lrd2,lrd3',
c    &                      ioff,irdbat,lrd1,lrd2,lrd3
c           call  flshfo(lupri)
            CALL GETWA2(LUTOC,FNTOC,WORK(KOFF),IOFF,LRD3)
            IOFF = IOFF + LRD3
            LRD2 = LRD2 + LRD3
            KOFF = KOFF + LRD3
c           write(lupri,*) 'irdbat,koff,ldr2,ldr3,ntotoc :',
c    &                      irdbat,koff,ldr2,ldr3,ntotoc(1)
c           call  flshfo(lupri)
         END DO
         IF (LRD2 .NE. NTOTOC(1)) THEN
            WRITE(LUPRI,*) 'Error reading CCSDT_OC'
            CALL QUIT('Error reading CCSDT_OC')
         END IF
         CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      ENDIF
C
C-------------------------------------------------
C     Transform to (ia|j k) (sorted  as (ji,k,a)).
C-------------------------------------------------
C
      CALL CCSDT_TROCC(WORK(KINTOC),OCINT,WORK(KCMO),WORK(KEND1),LWRK1)
C
C------------------------
C     Resort to (ij,k,a).
C------------------------
C
      KSCR  = 1
      KEND2 = KSCR  + NTRAOC(1)
      LWRK2 = LWORK - KEND2 + 1
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Error in CCHO_TRAINT:'
         WRITE(LUPRI,*)
     &   'Insufficient memory for resorting occupied integrals'
         WRITE(LUPRI,*) 'Need     : ',KEND2-1
         WRITE(LUPRI,*) 'Available: ',LWORK
         CALL QUIT('Insufficient memory in CCHO_TRAINT')
      ENDIF
C
      CALL DCOPY(NTRAOC(1),OCINT,1,WORK(KSCR),1)
      CALL CCHO_RSRTOI(WORK(KSCR),OCINT,1)
C
      RETURN
      END
C  /* Deck ccho_rsrtoi */
      SUBROUTINE CCHO_RSRTOI(XJIKA,XIJKA,ISYM)
C
C     JLC, BFR, TBP, HK, and AS, October 2002.
C
C     Purpose: resort occupied integrals (ia|jk) from
C              (ji,k,a) to (ij,k,a) ordering.
C
C              ISYM is the symmetry of the integrals.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XJIKA(*), XIJKA(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      DO ISYMA = 1,NSYM
C
         ISYIJK = MULD2H(ISYMA,ISYM)
         ISYJIK = ISYIJK
C
         DO ISYMK = 1,NSYM
C
            IF (NRHF(ISYMK) .LE. 0) GOTO 200
C
            ISYMIJ = MULD2H(ISYMK,ISYIJK)
            ISYMJI = ISYMIJ
C
            DO A = 1,NVIR(ISYMA)
               DO K = 1,NRHF(ISYMK)
C
                  DO ISYMJ = 1,NSYM
C
                     ISYMI = MULD2H(ISYMJ,ISYMIJ)
C
                     IF (NRHF(ISYMI) .LE. 0) GOTO 100
C
                     DO J = 1,NRHF(ISYMJ)
C
                        KJIKA = ISJIKA(ISYJIK,ISYMA)
     &                        + NMAJIK(ISYJIK)*(A - 1)
     &                        + ISJIK(ISYMJI,ISYMK)
     &                        + NMATIJ(ISYMJI)*(K - 1)
     &                        + IMATIJ(ISYMJ,ISYMI)
     &                        + J
                        KIJKA = ISJIKA(ISYIJK,ISYMA)
     &                        + NMAJIK(ISYIJK)*(A - 1)
     &                        + ISJIK(ISYMIJ,ISYMK)
     &                        + NMATIJ(ISYMIJ)*(K - 1)
     &                        + IMATIJ(ISYMI,ISYMJ)
     &                        + NRHF(ISYMI)*(J - 1)
     &                        + 1
C
                        CALL DCOPY(NRHF(ISYMI),XJIKA(KJIKA),NRHF(ISYMJ),
     &                                         XIJKA(KIJKA),1)
C
                     ENDDO
C
  100                CONTINUE
C
                  ENDDO
C
               ENDDO
            ENDDO
C
  200       CONTINUE
C
         ENDDO
C
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_trint1 */
      SUBROUTINE CCHO_TRINT1(CMO,FILNM1,LUNIT1,FILNM2,LUNIT2,
     &                       WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      CHARACTER*(*) FILNM1,FILNM2
C
      DIMENSION CMO(*)
      DIMENSION WORK(LWORK)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
C
      ISYINT = ISYMOP
C
C-------------------------------
C     Loop over symmetries of b.
C-------------------------------
C
      DO 100 ISYMB = 1,NSYM
C
         NVIRB = NVIR(ISYMB)
         IF (NVIRB .EQ. 0) GOTO 100
C
         ISYEMD = MULD2H(ISYINT,ISYMB)
C
         LENMIN = NCKA(ISYEMD) + NCKATR(ISYEMD)
         NUMB   = MIN(NVIRB,LWORK/LENMIN)
         IF (NUMB .EQ. 0) CALL QUIT('Not enough space in CCHO_TRINT1')
C
Casm     Apparently, it is not possible to read more than 2 Gb (268435456 dw)
C
         MXDAL1 = 268435456 / NCKA(ISYEMD)
         MXDAL2 = 268435456 / NCKATR(ISYEMD)
         MXDALF = MIN(MXDAL1,MXDAL2)
         NUMB   = MIN(NUMB,MXDALF)
C
         NBATB  = (NVIRB-1)/NUMB + 1
C
C------------------------
C        Batch structure.
C------------------------
C
         IB2 = 0
         DO IBATB = 1,NBATB
C
            IB1 = IB2 + 1
            IB2 = IB2 + NUMB
            IF (IB2 .GT. NVIRB) THEN
               IB2  = NVIRB
               NUMB = IB2 - IB1 + 1
            END IF
C
C-----------------------------
C           Memory allocation.
C-----------------------------
C
            KUNTRA = 1 
            KTRANS = KUNTRA + NCKA(ISYEMD)*NUMB
            KEND   = KTRANS + NCKATR(ISYEMD)*NUMB
C
            LWRK = LWORK - KEND
            IF (LWRK .LT. 0) THEN
               CALL QUIT('Not enough space in CCHO_TRINT1')
            ENDIF
C
C----------------------------------------
C           Read untransformed integrals.
C----------------------------------------
C
            NCKA8  = NCKA(ISYEMD)
            IB8    = IB1
            IOFF   = ICKAD(ISYEMD,ISYMB) + NCKA8*(IB8-1) + 1
            LENGTH = NCKA(ISYEMD)*NUMB
C
            IF (LENGTH .GT. 0) THEN
               CALL GETWA2(LUNIT1,FILNM1,WORK(KUNTRA),IOFF,LENGTH)
*              xnrm = dnorm2(length,work(kuntra),1)
*              write(lupri,*) LUNIT1,FILNM1
*              write(lupri,*) ioff,': ',xnrm
*              write(lupri,'(4F15.10)') (work(ii),ii=kuntra,kuntra+39)
            END IF
C
C---------------------------------------
C           Loop over orbitals in batch.
C---------------------------------------
C
            DO B = IB1,IB2
C
               IB = B - IB1 + 1
C
C--------------------------------------------
C              Transform each symmetry block.
C--------------------------------------------
C
               DO ISYMD = 1,NSYM
C
                  ISYMEM = MULD2H(ISYEMD,ISYMD)
C
                  NTOTEM = MAX(NT1AM(ISYMEM),1)
                  NBASD  = MAX(NBAS(ISYMD),1)
C
                  KOFF1 = KUNTRA + NCKA(ISYEMD)*(IB-1)
     &                  + ICKA(ISYMEM,ISYMD)
                  KOFF2 = ILMVIR(ISYMD) + 1
                  KOFF3 = KTRANS + NCKATR(ISYEMD)*(IB-1)
     &                  + ICKATR(ISYMEM,ISYMD)
C
                  CALL DGEMM('N','N',NT1AM(ISYMEM),NVIR(ISYMD),
     &                       NBAS(ISYMD),ONE,WORK(KOFF1),NTOTEM,
     &                       CMO(KOFF2),NBASD,ZERO,WORK(KOFF3),NTOTEM)
C
               END DO
C
            END DO
C
C---------------------------------------
C           Write transformed integrals.
C---------------------------------------
C
            NCKA8  = NCKATR(ISYEMD)
            IB8    = IB1
            IOFF   = ICKBD(ISYEMD,ISYMB) + NCKA8*(IB8-1) + 1
            LENGTH = NCKATR(ISYEMD)*NUMB
C
            IF (LENGTH .GT. 0) THEN
               CALL PUTWA2(LUNIT2,FILNM2,WORK(KTRANS),IOFF,LENGTH)
*              xnrm = dnorm2(length,work(ktrans),1)
*              write(lupri,*) LUNIT2,FILNM2
*              write(lupri,*) ioff,': ',xnrm
*              write(lupri,'(4F15.10)') (work(ii),ii=ktrans,ktrans+39)
            ENDIF
C
         END DO
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccho_t2a */
      SUBROUTINE CCHO_T2A(T2VO,T2AM,ISYMT2)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION T2VO(*),T2AM (*)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      DO ISYMIJ = 1,NSYM
         ISYMAB = MULD2H(ISYMIJ,ISYMT2)
         DO ISYMI = 1,NSYM
            ISYMJ = MULD2H(ISYMI,ISYMIJ)
            DO I = 1,NRHF(ISYMI)
               DO J = 1,NRHF(ISYMJ)
                  DO ISYMA = 1,NSYM
                     ISYMB  = MULD2H(ISYMA,ISYMAB)
                     ISYMAI = MULD2H(ISYMA,ISYMI)
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
                     IF (ISYMAI .NE. ISYMBJ) CALL QUIT('Problem')
                     DO A = 1,NVIR(ISYMA)
                        DO B = 1,NVIR(ISYMB)
                           NAI = IT1AM(ISYMA,ISYMI) 
     &                         + NVIR(ISYMA)*(I-1) + A
                           NBJ = IT1AM(ISYMB,ISYMJ) 
     &                         + NVIR(ISYMB)*(J-1) + B
                           IF (ISYMAI .EQ. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ) 
     &                              + INDEX(NAI,NBJ)
                           ELSE
                              IF (ISYMAI .LT. ISYMBJ) THEN
                                 NAIBJ = IT2AM(ISYMAI,ISYMBJ) 
     &                                 + NT1AM(ISYMAI)*(NBJ-1) + NAI
                              ELSE
                                 NAIBJ = IT2AM(ISYMBJ,ISYMAI) 
     &                                 + NT1AM(ISYMBJ)*(NAI-1) + NBJ
                              ENDIF
                           ENDIF
                           NIJ = IMATIJ(ISYMI,ISYMJ)
     &                         + NRHF(ISYMI)*(J-1) + I
                           NAB = IMATAB(ISYMA,ISYMB)
     &                         + NVIR(ISYMA)*(B-1) + A
                           NABIJ = IT2VO(ISYMAB,ISYMIJ) 
     &                           + NMATAB(ISYMAB)*(NIJ-1) + NAB
                           T2VO(NABIJ) = T2AM(NAIBJ)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_ttild1 */
      SUBROUTINE CCHO_TTILD1(T2VO,T2TIL,NA,IA1,ISYMA,NB,IB1,
     &                       ISYMB,NIJ,ISYMT2)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION T2VO(*),T2TIL(NIJ,NA,NB)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYMAB = MULD2H(ISYMA,ISYMB)
      ISYMIJ = MULD2H(ISYMAB,ISYMT2)
C
      DO B = 1,NB
         DO A = 1,NA
            DO IJ = 1,NIJ
               IA    = IA1 + A - 1
               IB    = IB1 + B - 1
               NAB   = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(IB-1) + IA
               NBA   = IMATAB(ISYMB,ISYMA) + NVIR(ISYMB)*(IA-1) + IB
               NABIJ = IT2VO(ISYMAB,ISYMIJ) 
     &               + NMATAB(ISYMAB)*(IJ-1) + NAB
               NBAIJ = IT2VO(ISYMAB,ISYMIJ) 
     &               + NMATAB(ISYMAB)*(IJ-1) + NBA
               T2TIL(IJ,A,B)=TWO*T2VO(NABIJ)-T2VO(NBAIJ)
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_ttild2 */
      SUBROUTINE CCHO_TTILD2(T2TIL,CHOVEC,NA,NB,NIJ)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      INTEGER A,B
C
      DIMENSION T2TIL(NIJ,NA,NB),CHOVEC(NIJ,NB)
C
C
      DO B = 1,NB
         DO A = 1,NA
            DO IJ = 1,NIJ
               T2TIL(IJ,A,B)=T2TIL(IJ,A,B)*CHOVEC(IJ,B)
            ENDDO
         ENDDO
      ENDDO

C
      RETURN
      END
C  /* Deck ccho_y1 */
      SUBROUTINE CCHO_Y1(XMAT,Y1MAT,Y2MAT,NA,NCD,NC,ISYMC,ND,ISYMD)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION XMAT(NCD,NA),Y1MAT(NC,NA,ND),Y2MAT(NC,NA,ND)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
C
      DO A = 1,NA
         DO D = 1,ND
            DO C = 1,NC
C
               ICD = IMATAB(ISYMC,ISYMD) + NC*(D-1) + C
               IDC = IMATAB(ISYMD,ISYMC) + ND*(C-1) + D
C
               Y1MAT(C,A,D) = XMAT(ICD,A)
               Y2MAT(C,A,D) = XMAT(IDC,A)
C
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C
C  /* Deck ccho_rdint */
      SUBROUTINE CCHO_RDINT(LUNIT,FILNAM,XINT,IB1,NUMIB,ISYMB)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      CHARACTER*7 FILNAM
      DIMENSION XINT(*)
C
C
      ISYEMD = MULD2H(ISYMB,ISYMOP)
C
      IOFF   = ICKBD(ISYEMD,ISYMB) + NCKATR(ISYEMD)*(IB1 - 1) + 1
      LENGTH = NCKATR(ISYEMD)*NUMIB
      IF (LENGTH .GT. 0) THEN
         CALL GETWA2(LUNIT,FILNAM,XINT,IOFF,LENGTH)
      ENDIF
C
      RETURN
      END
C
C  /* Deck ccho_rdkin */
      SUBROUTINE CCHO_RDKIN(LUNIT,FILNAM,XKINT,XLINT,WORK,LWORK,
     &                      IA1,NUMIA,ISYMA,IOFSET)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C     XINT is packed under symmetries of C. Inside each C symmetry
C     block, the packing is (em,c,a)
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      CHARACTER*7 FILNAM
      DIMENSION XKINT(*),XLINT(*),WORK(LWORK)
      DIMENSION IOFSET(8)
C
C
      ISYEMC = MULD2H(ISYMA,ISYMOP)
C
      DO A = 1,NUMIA
C
         IA = IA1 + A - 1
C
         IOFF   = ICKBD(ISYEMC,ISYMA) + NCKATR(ISYEMC)*(IA - 1) + 1
         LENGTH = NCKATR(ISYEMC)
         IF (LWORK .LT. LENGTH)
     &   CALL QUIT('Not enough space in CCHO_RDKIN')
C
         IF (LENGTH .GT. 0) THEN
            CALL GETWA2(LUNIT,FILNAM,WORK,IOFF,LENGTH)
         END IF
C
         DO ISYMC = 1,NSYM
C
            ISYMEM = MULD2H(ISYEMC,ISYMC)
            LENEMC = NT1AM(ISYMEM)*NVIR(ISYMC)
C
            KOFF1 = ICKATR(ISYMEM,ISYMC) + 1
            KOFF2 = IOFSET(ISYMC) + LENEMC*(A-1) + 1
C
            CALL DCOPY(LENEMC,WORK(KOFF1),1,XKINT(KOFF2),1)
C
            DO C = 1,NVIR(ISYMC)
               DO ISYMM = 1,NSYM
                  ISYME  = MULD2H(ISYMM,ISYMEM) 
                  ISYMCM = MULD2H(ISYMC,ISYMM)
                  DO M = 1,NRHF(ISYMM)
                     DO E = 1,NVIR(ISYME)
                        KOFF3 = IOFSET(ISYMC) 
     &                        + NT1AM(ISYMEM)*NVIR(ISYMC)*(A-1) 
     &                        + NT1AM(ISYMEM)*(C-1) + IT1AM(ISYME,ISYMM)
     &                        + NVIR(ISYME)*(M-1) + E
                        KOFF4 = ICKATR(ISYMCM,ISYME)+NT1AM(ISYMCM)*(E-1)
     &                        + IT1AM(ISYMC,ISYMM)+NVIR(ISYMC)*(M-1)+C
                        XLINT(KOFF3) = WORK(KOFF4)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
C
         END DO
      END DO
C
      RETURN
      END
C
C  /* Deck ccho_init */
      SUBROUTINE CCHO_INIT(FOCKD,CHOELE,MAXCHO,NUMCHO,WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C     CCSD(T) energy correction with Cholesky decomposition.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION FOCKD(*),CHOELE(MAXCHO)
C
      DIMENSION WORK(LWORK)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      NDIMP = (NRHFT+1)*NRHFT*NVIRT/2 + (NVIRT+1)*NVIRT*NRHFT/2
C
      KFOCC = 1
      KFVIR = KFOCC + NRHFT
      KOMEG = KFVIR + NVIRT
      KEND1 = KOMEG + NDIMP
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCHO_INIT')
      ENDIF
C
      IO = 0
      IV = 0
      DO ISYM = 1,NSYM
         DO I = 1,NRHF(ISYM)
            KOFF1 = IRHF(ISYM) + I
            WORK(KFOCC+IO) = FOCKD(KOFF1)
            IO = IO + 1
         ENDDO
         DO A = 1,NVIR(ISYM)
            KOFF2 = IVIR(ISYM) + A
            WORK(KFVIR+IV) = FOCKD(KOFF2)
            IV = IV + 1
         ENDDO
      ENDDO
C
      CALL CCHO_INIT1(WORK(KFOCC),WORK(KFVIR),WORK(KOMEG),CHOELE,
     *                MAXCHO,NUMCHO)
C
      RETURN
      END
C
C  /* Deck ccho_init1 */
      SUBROUTINE CCHO_INIT1(FOCC,FVIR,OMEG,CHOELE,MAXCHO,NUMCHO)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C     CCSD(T) energy correction with Cholesky decomposition.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-15)
C
      DIMENSION FOCC(*),FVIR(*),OMEG(*),CHOELE(MAXCHO)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      NDIMP = (NRHFT+1)*NRHFT*NVIRT/2 + (NVIRT+1)*NVIRT*NRHFT/2
C
      IND = 0
      DO I = 1,NRHFT
         DO J = 1,I
            DO B = 1,NVIRT
               IND = IND + 1
               OMEG(IND) = FVIR(B)-FOCC(I)-FOCC(J)
            ENDDO
         ENDDO
      ENDDO 
      DO I = 1,NRHFT
         DO A = 1,NVIRT
            DO B = 1,A
               IND = IND + 1
               OMEG(IND) = FVIR(A) + FVIR(B) - FOCC(I)
            ENDDO
         ENDDO
      ENDDO
C
      DO N = 1,MAXCHO
         Y = ZERO
         DO P = 1,NDIMP
            X = ONE/(TWO*OMEG(P))
            DO M = 1,N-1
               X = X*((OMEG(P)-CHOELE(M))/(OMEG(P)+CHOELE(M)))**TWO
            ENDDO
            IF (X .GT. Y) THEN
               Y = X
               CHOELE(N) = OMEG(P)
            ENDIF
         ENDDO
         NUMCHO = N
         IF (Y .LT. THRES) THEN
C
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Cholesky decomposition'
            WRITE(LUPRI,*) '----------------------'
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) '# vectors: ',NUMCHO
            WRITE(LUPRI,*) 'Maximum diagonal element',Y
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)
C
            RETURN 
         END IF
      ENDDO
C
      RETURN
      END
C
C  /* Deck ccho_vec */
      SUBROUTINE CCHO_VEC(FOCKD,CHOELE,CHOVEC,ICHO,NUMCHO,NUMIB,IB1,
     &                    ISYMB,ISYMA,ISYEMC,ISYMIJ,ISYIJC)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C     Construct vector for Cholesky decomposition.
C     The order of elements is:
C       1.- ijb for fixed #b from the batch over b (numib*nij)
C       2.- emc for fixed ISYEMC  (nckatr(isyemc))
C       3.- ijc for fixed ISYIJC  (nmaija(isyijc))
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),CHOVEC(*)
      DIMENSION IOF2(8)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
      NIJ = NMATIJ(ISYMIJ)
C
      IF (ICHO .EQ. 1) THEN
C
C--------------------------------
C     Calculate the first vector.
C--------------------------------
C
C     First part
C     ----------
C
      IND = 0
      DO B = 1,NUMIB
         KOFFB = IVIR(ISYMB) + IB1 + B - 1
         DO ISYMJ = 1,NSYM
            ISYMI = MULD2H(ISYMIJ,ISYMJ)
            DO J = 1,NRHF(ISYMJ)
               KOFFJ = IRHF(ISYMJ) + J
               DO I = 1,NRHF(ISYMI)
                  KOFFI = IRHF(ISYMI) + I
                  OME = FOCKD(KOFFB) - FOCKD(KOFFI) - FOCKD(KOFFJ)
C
                  IND = IND + 1
                  CHOVEC(IND) = SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
               END DO
            END DO
         END DO
      END DO
C
C     Second part.
C     ------------
C
      DO ISYMC = 1,NSYM
         ISYMEM = MULD2H(ISYEMC,ISYMC)
         DO C = 1,NVIR(ISYMC)
            KOFFC = IVIR(ISYMC) + C
            DO ISYMM = 1,NSYM
               ISYME = MULD2H(ISYMEM,ISYMM)
               DO M = 1,NRHF(ISYMM)
                  KOFFM = IRHF(ISYMM) + M
                  DO E = 1,NVIR(ISYME)
                     KOFFE = IVIR(ISYME) + E
                     OME = FOCKD(KOFFC) + FOCKD(KOFFE) - FOCKD(KOFFM)
C
                     IND = IND + 1
                     CHOVEC(IND) = SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C     Third part.
C     -----------
C
      DO ISYMC = 1,NSYM
         JSYMIJ = MULD2H(ISYIJC,ISYMC)
         DO C = 1,NVIR(ISYMC)
            KOFFC = IVIR(ISYMC) + C
            DO ISYMJ = 1,NSYM
               ISYMI = MULD2H(JSYMIJ,ISYMJ)
               DO J = 1,NRHF(ISYMJ)
                  KOFFJ = IRHF(ISYMJ) + J
                  DO I = 1,NRHF(ISYMI)
                     KOFFI = IRHF(ISYMI) + I
                     OME = FOCKD(KOFFC) - FOCKD(KOFFI) - FOCKD(KOFFJ)
C
                     IND = IND + 1
                     CHOVEC(IND) = SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      ELSE
C
C------------------------------------
C     Construct the updating vector.
C------------------------------------
C
C     First part.
C     -----------
C
      IND = 0
      DO B = 1,NUMIB
         KOFFB = IVIR(ISYMB) + IB1 + B - 1
         DO ISYMJ = 1,NSYM
            ISYMI = MULD2H(ISYMIJ,ISYMJ)
            DO J = 1,NRHF(ISYMJ)
               KOFFJ = IRHF(ISYMJ) + J
               DO I = 1,NRHF(ISYMI)
                  KOFFI = IRHF(ISYMI) + I
                  OME = FOCKD(KOFFB) - FOCKD(KOFFI) - FOCKD(KOFFJ)
C
                  IND = IND + 1
                  CHOVEC(IND) = (OME-CHOELE(ICHO-1))/(OME+CHOELE(ICHO))
C
               END DO
            END DO
         END DO
      END DO
C
C     Second part.
C     ------------
C
      DO ISYMC = 1,NSYM
         ISYMEM = MULD2H(ISYEMC,ISYMC)
         DO C = 1,NVIR(ISYMC)
            KOFFC = IVIR(ISYMC) + C
            DO ISYMM = 1,NSYM
               ISYME = MULD2H(ISYMEM,ISYMM)
               DO M = 1,NRHF(ISYMM)
                  KOFFM = IRHF(ISYMM) + M
                  DO E = 1,NVIR(ISYME)
                     KOFFE = IVIR(ISYME) + E
                     OME = FOCKD(KOFFC) + FOCKD(KOFFE) - FOCKD(KOFFM)
C
                     IND = IND + 1
                     CHOVEC(IND) = (OME-CHOELE(ICHO-1))
     &                           / (OME+CHOELE(ICHO))
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C     Third part.
C     -----------
C
      DO ISYMC = 1,NSYM
         JSYMIJ = MULD2H(ISYIJC,ISYMC)
         DO C = 1,NVIR(ISYMC)
            KOFFC = IVIR(ISYMC) + C
            DO ISYMJ = 1,NSYM
               ISYMI = MULD2H(JSYMIJ,ISYMJ)
               DO J = 1,NRHF(ISYMJ)
                  KOFFJ = IRHF(ISYMJ) + J
                  DO I = 1,NRHF(ISYMI)
                     KOFFI = IRHF(ISYMI) + I
                     OME = FOCKD(KOFFC) - FOCKD(KOFFI) - FOCKD(KOFFJ)
C
                     IND = IND + 1
                     CHOVEC(IND) = (OME-CHOELE(ICHO-1))
     &                           / (OME+CHOELE(ICHO))
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      NDIM = IND
      FACT = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
      CALL DSCAL(NDIM,FACT,CHOVEC,1)
C
      END IF
C
      RETURN
      END
C
C  /* Deck ccho_sclint */
      SUBROUTINE CCHO_SCLINT(XKINT,XLINT,CHOVEC,NA,NB)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XKINT(NB,NA),XLINT(NB,NA),CHOVEC(NB) 
C
      DO I = 1,NA
         DO J = 1,NB
            XKINT(J,I) = XKINT(J,I)*CHOVEC(J)
            XLINT(J,I) = XLINT(J,I)*CHOVEC(J)
         ENDDO
      ENDDO
C
      RETURN
      END
C
C  /* Deck ccho_tbijc */
      SUBROUTINE CCHO_TBIJC(T2VO,T2BIJC,NB,IB1,ISYMB,ISYIJC,NIJC)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION T2VO(*),T2BIJC(NB,NIJC)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DO ISYMC = 1,NSYM
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYMIJ = MULD2H(ISYMC,ISYIJC)
         NIJ = NMATIJ(ISYMIJ)
         DO C = 1,NVIR(ISYMC)
            DO IJ = 1,NIJ
               IJC = IMAIJA(ISYMIJ,ISYMC) 
     &             + NMATIJ(ISYMIJ)*(C-1) + IJ
               DO B = 1,NB
                  IB = IB1 + B - 1
C
                  NBC   = IMATAB(ISYMB,ISYMC) + NVIR(ISYMB)*(C-1) + IB
                  NBCIJ = IT2VO(ISYMBC,ISYMIJ) 
     &                  + NMATAB(ISYMBC)*(IJ-1) + NBC
C
                  T2BIJC(B,IJC) = T2VO(NBCIJ)
C
               END DO
            END DO
         END DO
      END DO

C
      RETURN
      END
C
C  /* Deck ccho_tijca */
      SUBROUTINE CCHO_TIJCA(T2VO,T2IJCA,NA,IA1,ISYMA,ISYIJC,NIJC)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION T2VO(*),T2IJCA(NIJC,NA)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DO A = 1,NA
         IA  = IA1 + A - 1
         DO ISYMC = 1,NSYM
            ISYMAC = MULD2H(ISYMA,ISYMC)
            ISYMIJ = MULD2H(ISYMC,ISYIJC)
            NIJ = NMATIJ(ISYMIJ)
            DO C = 1,NVIR(ISYMC)
               NAC = IMATAB(ISYMA,ISYMC) + NVIR(ISYMA)*(C-1) + IA
               NCA = IMATAB(ISYMC,ISYMA) + NVIR(ISYMC)*(IA-1) + C
               DO IJ = 1,NIJ
                  IJC   = IMAIJA(ISYMIJ,ISYMC) 
     &                  + NMATIJ(ISYMIJ)*(C-1) + IJ
C
                  NACIJ = IT2VO(ISYMAC,ISYMIJ) 
     &                  + NMATAB(ISYMAC)*(IJ-1) + NAC
                  NCAIJ = IT2VO(ISYMAC,ISYMIJ) 
     &                  + NMATAB(ISYMAC)*(IJ-1) + NCA
C
                  T2IJCA(IJC,A) = TWO*T2VO(NACIJ) - T2VO(NCAIJ)
C
               END DO
            END DO
         END DO
      END DO

C
      RETURN
      END
C
C  /* Deck ccho_at2sc */
      SUBROUTINE CCHO_AT2SC(NIJC,NA,CHOVEC,T2IJCA)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      INTEGER A
C
      DIMENSION CHOVEC(NIJC),T2IJCA(NIJC,NA)
C
C
      DO A = 1,NA
         DO IJC = 1,NIJC
            T2IJCA(IJC,A) = T2IJCA(IJC,A)*CHOVEC(IJC)
         END DO
      END DO
C
      RETURN
      END
C
C  /* Deck ccho_aterm */
      SUBROUTINE CCHO_ATERM(AMAT,XJINT,XKINT,XLINT,SCRINT,XAMAT,TBIJC,
     &                      TIJCA,NA,NB,NIJC,ISYMA,ISYMB,ISYMIJC,ACOR)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION AMAT(NB,NA),XAMAT(NB,NA),TBIJC(NB,NIJC),TIJCA(NIJC,NA)
      DIMENSION XJINT(*),XKINT(*),XLINT(*),SCRINT(*)
      DIMENSION IOFSET(8)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
C
C-----------------
C     Construct X.
C-----------------
C
      CALL DGEMM('N','N',NB,NA,NIJC,ONE,TBIJC,NB,TIJCA,NIJC,
     &           ZERO,XAMAT,NB)
C
      ISYEMF = MULD2H(ISYMOP,ISYMB)
      NUMEMF = NCKATR(ISYEMF)
C
      ICOUNT = 0
      DO ISYMF = 1,NSYM
         ISYMEM = MULD2H(ISYEMF,ISYMF)
         IOFSET(ISYMF) = ICOUNT
         ICOUNT = ICOUNT + NVIR(ISYMF)*NT1AM(ISYMEM)*NA
      END DO
C
      DO A = 1,NA
C
C---------------------------------------------------
C        Resort KINT and LINT in SCRINT to get TCME.
C---------------------------------------------------
C
C
         DO ISYMF = 1,NSYM
C
            ISYMEM = MULD2H(ISYEMF,ISYMF)
            LENEMF = NVIR(ISYMF)*NT1AM(ISYMEM)
C
            KOFF1 = IOFSET(ISYMF) + LENEMF*(A-1) + 1
            KOFF2 = ICKATR(ISYMEM,ISYMF) + 1
C
            CALL DCOPY(LENEMF,XKINT(KOFF1),1,SCRINT(KOFF2),1)
C
            CALL DSCAL(LENEMF,TWO,SCRINT(KOFF2),1)
            CALL DAXPY(LENEMF,-ONE,XLINT(KOFF1),1,SCRINT(KOFF2),1)
C
         END DO
C
C-----------------------------
C        Construct AMAT(NB,A).
C-----------------------------
C
         CALL DGEMV('T',NUMEMF,NB,ONE,XJINT,NUMEMF,SCRINT,1,
     &               ZERO,AMAT(1,A),1)
C
      END DO
C
C--------------------------------
C     Add contribution to energy.
C--------------------------------
C
      NAB = NA*NB
      ACOR = - DDOT(NAB,AMAT,1,XAMAT,1)
C
      RETURN
      END
C
C  /* Deck ccho_eterm */
      SUBROUTINE CCHO_ETERM(XIINT,XJINT,T2VO,FOCKD,NUMCHO,CHOELE,
     *                      WORK,LWORK,E4DE,IB1,ISYMB,NUMIB,
     *                      FBATCH,PRINT,NONI)
C
C     Javier Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     Calculate (part of E-terms):
C
C        E4DE = Sum(b) Sum(cia) [ V(ci,a;b) * D(ai,c;b)
C                               + W(ci,a;b) * E(ai,c;b) ]
C
C     where
C
C        V(ci,a;b) = - Sum(em) d(em,b) * (em|ba) * s(em,ci)
C                    - Sum(em) d(em,b) * (bm|ea) * s(ei,cm)
C
C        W(ci,a;b) = 2 Sum(em) d(em,b) * (em|ba) * s(em,ci)
C                    - Sum(em) d(em,b) * (bm|ea) * s(em,ci)
C
C        D(ai,c;b) =   Sum(dl) d(dl,i) * (dl|bc) * s(dl,ai)
C                    - Sum(dl) d(dl,i) * (bl|dc) * t(dl,ai)
C
C        E(ai,c;b) = - Sum(dl) d(dl,i) * (bl|dc) * t(di,al)
C
C    and s(ai,bj) = 2 t(ai,bj) - t(aj,bi), and d(ai,j) and d(ai,b) denotes
C    the occupied and virtual parts of the Cholesky decomposition of the
C    orbital energy denominator.
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XIINT(*),XJINT(*),T2VO(*),CHOELE(*),FOCKD(*)
      DIMENSION WORK(LWORK)

c      DIMENSION IOFF2(8),IOFF3(8),IOFF4(8)
      DIMENSION IOFF2(8),IOFF3(8)
      LOGICAL FBATCH,PRINT
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C 
      TLAST = SECOND()
C 
C
      IF (PRINT) THEN
         TIMT = SECOND()
         WRITE(LUPRI,'(6X,A,/,6X,A,/)')
     &   'Calculation of the E1 term:',
     &   '==========================='
      ENDIF
C
      ISYMT2 = 1
C
C     I(em,d,b) = (bm|ed) and J(em,d,b) = (em|bd) 
C     -------------------------------------------
C
      ISINTI = ISYMB
      ISINTJ = ISYMB
C
      NTOTT1 = 0
      DO ISYM = 1,NSYM
         NTOTT1 = NTOTT1 + NT1AM(ISYM)
      ENDDO
C
      MAXAC = -1
      DO ISYMA = 1,NSYM
         DO ISYMC = 1,NSYM
            IAC = NVIR(ISYMA)*NVIR(ISYMC)
            IF (MAXAC .LT. IAC) MAXAC = IAC
         ENDDO
      ENDDO
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KCHOV = 1
      KINTI = KCHOV + NTOTT1*NUMIB 
      KINTJ = KINTI + NCKATR(ISINTI)
      KEND1 = KINTJ + NCKATR(ISINTJ)
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in CCHO_ETERM')
      ENDIF
C
      DO 100 ISYMI = 1,NSYM
C
         IF (NRHF(ISYMI) .EQ. 0) GOTO 100
C
C        Symmetry of amplitudes subblocks
C        --------------------------------
C
         ISYT21 = ISYMI
         ISYT22 = ISYMI
         ISYT23 = ISYMI
         ISYT24 = ISYMI
         ISYT25 = ISYMI
C
         ISYMBI = MULD2H(ISYMI,ISYMB)
         ISYMAC = ISYMBI
C
C-------------------------------
C        Batch over the I index.
C-------------------------------
C
         LENWRK = 5*NCKATR(ISYT21) + 4*MAXAC + NTOTT1
         IF (FBATCH) THEN
            NEFI  = MIN(NONI,NRHF(ISYMI))
            LEFF  = NEFI*LENWRK + 1
            LWRKB = MIN(LWRK1,LEFF)
         ELSE
            LWRKB = LWRK1
         ENDIF
         NUMI = MIN(NRHF(ISYMI),LWRKB/LENWRK)
         IF (NUMI .EQ. 0) THEN
            WRITE(LUPRI,*) 'NUMI .EQ. 0 in CC_ETERM'
            CALL QUIT('Not enough space in CC_ETERM')
         END IF
C
         NBATI = (NRHF(ISYMI)-1)/NUMI + 1
C
         IF (PRINT) THEN
            WRITE(LUPRI,'(9X,A,I1,A,/,9X,A)')
     &      'Batch over I, symmetry ',ISYMI,':',
     &      '-------------------------'
            WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10)')
     &      'Minimum work space required   : ',LENWRK,
     &      'Work space available for batch: ',LWRKB
            WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10,/)')
     &      'Number of occupied orbitals   : ',NRHF(ISYMI),
     &      'Required number of I-batches  : ',NBATI
         ENDIF
C
         II2 = 0
         DO IBATI = 1,NBATI
C
            II1 = II2 + 1
            II2 = II2 + NUMI
            IF (II2 .GT. NRHF(ISYMI)) II2 = NRHF(ISYMI)
            NUMII = II2 - II1 + 1
C
            IF (PRINT) THEN
               WRITE(LUPRI,'(12X,A,I10,A,/,12X,A)')
     &         'I-batch number ',IBATI,':',
     &         '--------------------------'
               WRITE(LUPRI,'(12X,A,I10,1X,I10,/)')
     &         'First and last I: ',II1,II2
            ENDIF
C
C-----------------------------
C           Some index arrays.
C-----------------------------
C
            ICOUN2 = 0
            ICOUN3 = 0
c            ICOUN4 = 0
            DO ISYMC = 1,NSYM
C
               ISYMCI = MULD2H(ISYMC,ISYMI)
               ISYMEM = MULD2H(ISYMCI,ISYMT2)
c               ISYMA  = MULD2H(ISYMC,ISYMAC)
C
               IOFF2(ISYMC) = ICOUN2
               IOFF3(ISYMC) = ICOUN3
c               IOFF4(ISYMC) = ICOUN4
C
               ICOUN2 = ICOUN2 + NT1AM(ISYMC)
               ICOUN3 = ICOUN3 + NT1AM(ISYMEM)*NVIR(ISYMC)*NUMII
c               ICOUN4 = ICOUN4 + NVIR(ISYMA)*NVIR(ISYMC)*NUMII
C
            END DO
C
C------------------------------
C           Dynamic allocation.
C------------------------------
C
            KCHOO = KEND1
            KT21  = KCHOO + NTOTT1*NUMII
            KT22  = KT21  + NCKATR(ISYT21)*NUMII
            KT23  = KT22  + NCKATR(ISYT22)*NUMII
            KT24  = KT23  + NCKATR(ISYT23)*NUMII
            KT25  = KT24  + NCKATR(ISYT24)*NUMII
            KVMAT = KT25  + NCKATR(ISYT25)*NUMII
            KWMAT = KVMAT + MAXAC*NUMII
            KDMAT = KWMAT + MAXAC*NUMII
            KEMAT = KDMAT + MAXAC*NUMII
            KEND2 = KEMAT + MAXAC*NUMII
C
            LWRK2 = LWORK - KEND2 + 1
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient memory in CCHO_ETERM')
            ENDIF
C
C-----------------------------------
C           Extract relevant arrays.
C-----------------------------------
C
            CALL CCHO_DET2(T2VO,WORK(KT21),WORK(KT22),WORK(KT23),
     &                     WORK(KT24),
     &                     II1,NUMII,ISYMI,IOFF3,ISYMT2)
            CALL DCOPY(NCKATR(ISYT21)*NUMII,WORK(KT21),1,WORK(KT25),1)
C
            DO ICHO = 1,NUMCHO
C
               DECOR = ZERO
C
C---------------------------------------------
C              Calculate Cholesky information.
C---------------------------------------------
C
               CALL CCHO_DECHO(FOCKD,CHOELE,NUMCHO,ICHO,
     &                         WORK(KCHOO),WORK(KCHOV),IB1,NUMIB,ISYMB,
     &                         II1,NUMII,ISYMI)
C
C-----------------------------------------------------------
C              Scale the T2 amplitudes with cholesky update.
C-----------------------------------------------------------
C
               CALL CCHO_SCT2O(WORK(KT21),WORK(KT23),WORK(KT24),
     &                         WORK(KCHOO),NUMII,ISYT21,
     &                         NTOTT1,IOFF2,IOFF3)
C
               DO B = 1,NUMIB
C
C-----------------------------------
C                 Extract integrals.
C-----------------------------------
C
                  CALL CCHO_DEIN(XIINT,XJINT,WORK(KINTI),WORK(KINTJ),
     &                           ISINTI,B)
C
C---------------------------------
C                 Scale integrals.
C---------------------------------
C
                  CALL CCHO_DESCL(WORK(KINTI),WORK(KINTJ),WORK(KCHOV),
     &                            ISINTI,B,IOFF2,NTOTT1)
C
                  DO ISYMA = 1,NSYM
C
                     ISYMC  = MULD2H(ISYMA,ISYMAC)
                     ISYMEM = MULD2H(ISYMA,ISINTI)
                     ISYMDL = MULD2H(ISYMC,ISINTI)
C
C--------------------------------------------------------
C                    Calculate the V and W intermediates.
C--------------------------------------------------------
C
                     KOFT5 = KT25  + IOFF3(ISYMC)
                     KOFFJ = KINTJ + ICKATR(ISYMEM,ISYMA)
C
                     NEM   = NT1AM(ISYMEM)
                     NTOEM = MAX(NEM,1)
                     NCI   = NVIR(ISYMC)*NUMII
                     NTOCI = MAX(NCI,1)
C
                     CALL DGEMM('T','N',NCI,NVIR(ISYMA),NEM,
     &                          ONE,WORK(KOFT5),NTOEM,WORK(KOFFJ),NEM,
     &                          ZERO,WORK(KVMAT),NTOCI)
C
                     NCIA  = NCI*NVIR(ISYMA)
                     CALL DCOPY(NCIA,WORK(KVMAT),1,WORK(KWMAT),1)
C
                     KOFT2 = KT22  + IOFF3(ISYMC)
                     KOFFI = KINTI + ICKATR(ISYMEM,ISYMA)
C
                     CALL DGEMM('T','N',NCI,NVIR(ISYMA),NEM,
     &                          XMONE,WORK(KOFT2),NTOEM,WORK(KOFFI),NEM,
     &                          XMONE,WORK(KVMAT),NTOCI)
C
                     CALL DGEMM('T','N',NCI,NVIR(ISYMA),NEM,
     &                          XMONE,WORK(KOFT5),NTOEM,WORK(KOFFI),NEM,
     &                          TWO,WORK(KWMAT),NTOCI)
C
C--------------------------------------------------------
C                    Calculate the D and E intermediates.
C--------------------------------------------------------
C
                     KOFT1 = KT21 + IOFF3(ISYMA)
                     KOFFJ = NCKATR(ISINTJ)*(B - 1)
     &                     + ICKATR(ISYMDL,ISYMC) + 1
C
                     NDL   = NT1AM(ISYMDL)
                     NTODL = MAX(NDL,1)
                     NAI   = NVIR(ISYMA)*NUMII
                     NTOAI = MAX(NAI,1)
C
                     CALL DGEMM('T','N',NAI,NVIR(ISYMC),NDL,
     &                          ONE,WORK(KOFT1),NTODL,
     &                          XJINT(KOFFJ),NTODL,
     &                          ZERO,WORK(KDMAT),NTOAI)
C
                     KOFT3 = KT23 + IOFF3(ISYMA)
                     KOFFI = NCKATR(ISINTI)*(B - 1)
     &                     + ICKATR(ISYMDL,ISYMC) + 1
C
                     CALL DGEMM('T','N',NAI,NVIR(ISYMC),NDL,
     &                          XMONE,WORK(KOFT3),NTODL,
     &                          XIINT(KOFFI),NTODL,
     &                          ONE,WORK(KDMAT),NTOAI)
C
                     KOFT4 = KT24 + IOFF3(ISYMA)
C
                     CALL DGEMM('T','N',NAI,NVIR(ISYMC),NDL,
     &                          XMONE,WORK(KOFT4),NTODL,
     &                          XIINT(KOFFI),NTODL,
     &                          ZERO,WORK(KEMAT),NTOAI)
C
C--------------------------------------------------
C                    Calculate energy contribution.
C--------------------------------------------------
C
                     DO A = 1,NVIR(ISYMA)
                        DO I = 1,NUMII
                           DO C = 1,NVIR(ISYMC)
C
                              NCIA = NVIR(ISYMC)*NUMII*(A - 1)
     &                             + NVIR(ISYMC)*(I - 1) + C
                              NAIC = NVIR(ISYMA)*NUMII*(C - 1)
     &                             + NVIR(ISYMA)*(I - 1) + A
C
                              KOFFV = KVMAT + NCIA - 1
                              KOFFD = KDMAT + NAIC - 1
                              KOFFW = KWMAT + NCIA - 1
                              KOFFE = KEMAT + NAIC - 1
C
                              DECOR = DECOR + WORK(KOFFV)*WORK(KOFFD) 
     &                                      + WORK(KOFFW)*WORK(KOFFE)
C
                           ENDDO
                        ENDDO
                     ENDDO
C
                  ENDDO
C
               ENDDO
C
               E4DE = E4DE + DECOR 
C 
               ENERGE(ICHO) = ENERGE(ICHO) + DECOR
C 
C Decommented by Domenico
              TNOW = SECOND()
              DELTAT = TNOW - TLAST
              TLAST = TNOW
              SCNDSE(ICHO) = SCNDSE(ICHO) 
     &             + DELTAT
C Decommented by Domenico 
C
c              WRITE(LUPRI,'(A,4I4,2D15.6)')'ISC,ISI,IBI,IV,DCO,DT',
c    &                  ISYMC,ISYMI,IBATI,ICHO,DECOR,E4DT
C
c              IF (PRINT) THEN
c                 WRITE(LUPRI,'(15X,A,I3,A,/,15X,A)')
c    &            'Status after Cholesky vector',ICHO,':',
c    &            '--------------------------------'
c                 IF (ABS(DECOR) .LT. THRCHO) THEN
c                    WRITE(LUPRI,'(15X,A)') 'E1 term converged'
c                 ELSE
c                    WRITE(LUPRI,'(15X,A)') 'E1 term not converged'
c                 ENDIF
c                 TIM = SECOND() - TIMT
c                 WRITE(LUPRI,'(15X,A,F10.2,A,/)')
c    &            'Accumulated E1-time: ',TIM,' seconds'
c              ENDIF
C
               IF (ABS(DECOR) .LT. THRCHO) GOTO 999 
C
            ENDDO 
C
  999       CONTINUE
C
         ENDDO
  100 CONTINUE
C
      RETURN
      END
C
C  /* Deck ccho_det2 */
      SUBROUTINE CCHO_DET2(T2VO,T21,T22,T23,T24,II1,NUMII,ISYMI,IOFF3,
     &                     ISYMT2)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (TWO = 2.0D0)
C
      DIMENSION T2VO(*),T21(*),T22(*), T23(*), T24(*)
      DIMENSION IOFF3(8)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DO I = 1,NUMII
C
         II = II1 + I - 1
C
         DO ISYMC = 1,NSYM
C
            ISYMCI = MULD2H(ISYMC,ISYMI)
            ISYMEM = MULD2H(ISYMCI,ISYMT2)
C
            DO C = 1,NVIR(ISYMC)
               DO ISYMM = 1,NSYM
C
                  ISYME  = MULD2H(ISYMEM,ISYMM)
                  ISYMEC = MULD2H(ISYME,ISYMC)
                  ISYMMI = MULD2H(ISYMM,ISYMI)
C
                  DO M = 1,NRHF(ISYMM)
                     DO E = 1,NVIR(ISYME)
C
                        NEC = IMATAB(ISYME,ISYMC) 
     &                      + NVIR(ISYME)*(C-1) + E
                        NCE = IMATAB(ISYMC,ISYME)
     &                      + NVIR(ISYMC)*(E-1) + C
                        NMI = IMATIJ(ISYMM,ISYMI)
     &                      + NRHF(ISYMM)*(II-1) + M
C
                        NCI = NVIR(ISYMC)*(I-1) + C
                        NEM = IT1AM(ISYME,ISYMM) + NVIR(ISYME)*(M-1) + E
C
                        NECMI = IT2VO(ISYMEC,ISYMMI)
     &                        + NMATAB(ISYMEC)*(NMI-1) + NEC
                        NCEMI = IT2VO(ISYMEC,ISYMMI)
     &                        + NMATAB(ISYMEC)*(NMI-1) + NCE
C
                        NEMCI = IOFF3(ISYMC) + NT1AM(ISYMEM)*(NCI-1)
     &                        + NEM
C
                        T21(NEMCI) = TWO*T2VO(NECMI) - T2VO(NCEMI)
                        T22(NEMCI) = TWO*T2VO(NCEMI) - T2VO(NECMI)
                        T23(NEMCI) = T2VO(NECMI)
                        T24(NEMCI) = T2VO(NCEMI)
C
                     END DO
                  END DO
C
               END DO
            END DO
C
         END DO
C
      END DO
C
C
      RETURN
      END
C
C  /* Deck ccho_dein */
      SUBROUTINE CCHO_DEIN(XIINT,XJINT,XINTD,XINTW,ISINTD,C)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION XIINT(*),XJINT(*),XINTD(*),XINTW(*)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
C---------------------------
C     Extract the integrals.
C---------------------------
C
      LEN  = NCKATR(ISINTD)
      KOFF = NCKATR(ISINTD)*(C-1) + 1
C
      CALL DCOPY(LEN,XIINT(KOFF),1,XINTD,1)
      CALL DCOPY(LEN,XJINT(KOFF),1,XINTW,1)
C
C      CALL DSCAL(LEN,TWO,XINTW,1)
C      CALL DAXPY(LEN,-ONE,XINTD,1,XINTW,1)
C
      RETURN
      END
C
C  /* Deck ccho_decho */
      SUBROUTINE CCHO_DECHO(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                      IC1,NUMIC,ISYMC,II1,NUMII,ISYMI)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C     Construct vector for Cholesky decomposition:
C
C     For the occupied part, OCCHO(emi), only the update vector
C     For the virtual  part, VICHO(emc), the actual vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C-------------------
C     Occupied part.
C-------------------
C
      IF (ICHO .EQ. 1) THEN
C
C        First vector.
C        --------------
C
         IND = 0
         DO I = 1,NUMII
            KOFFI = IRHF(ISYMI) + II1 + I - 1
            DO ISYMEM = 1,NSYM
               DO ISYMM = 1,NSYM
                  ISYME = MULD2H(ISYMM,ISYMEM)
                  DO M = 1,NRHF(ISYMM)
                     KOFFM = IRHF(ISYMM) + M
                     DO E = 1,NVIR(ISYME)
                        KOFFE = IVIR(ISYME) + E
                        OME = FOCKD(KOFFE)-FOCKD(KOFFI)-FOCKD(KOFFM)
C
                        IND = IND + 1
                        OCCHO(IND) = 
     &                       SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     ENDDO
                  END DO
               END DO
            END DO
         END DO
C
      ELSE
C
C        Updating vector.
C        ----------------
C
         IND = 0
         DO I = 1,NUMII
            KOFFI = IRHF(ISYMI) + II1 + I - 1
            DO ISYMEM = 1,NSYM
               DO ISYMM = 1,NSYM
                  ISYME = MULD2H(ISYMM,ISYMEM)
                  DO M = 1,NRHF(ISYMM)
                     KOFFM = IRHF(ISYMM) + M
                     DO E = 1,NVIR(ISYME)
                        KOFFE = IVIR(ISYME) + E
                        OME = FOCKD(KOFFE)-FOCKD(KOFFI)-FOCKD(KOFFM)
C
                        IND = IND + 1
                        OCCHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     ENDDO
                  END DO
               END DO
            END DO
         END DO
C
         NDIMOC = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMOC,FACTOR,OCCHO,1)
C
      END IF
C
C------------------
C     Virtual part.
C------------------
C
C     Get the omega's.
C     ----------------
C
      IND = 0
      DO C = 1,NUMIC
         KOFFC = IVIR(ISYMC) + IC1 + C - 1
         DO ISYMEM = 1,NSYM
            DO ISYMM = 1,NSYM
               ISYME = MULD2H(ISYMEM,ISYMM)
               DO M = 1,NRHF(ISYMM)
                  KOFFM = IRHF(ISYMM) + M
                  DO E = 1,NVIR(ISYME)
                     KOFFE = IVIR(ISYME) + E
C
                     IND = IND + 1
                     VICHO(IND) = FOCKD(KOFFC) + FOCKD(KOFFE) 
     &                          - FOCKD(KOFFM)
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      NDIMVI = IND
C
C
C     Construct the vector.
C     ---------------------
C
      DO P = 1,NDIMVI
C
         OMEGA = VICHO(P)
C
         VICHO(P) = SQRT(TWO*CHOELE(ICHO))/(OMEGA+CHOELE(ICHO))
C
         DO JCHO = 1,ICHO-1
            VICHO(P) = VICHO(P) 
     &               *(OMEGA-CHOELE(JCHO))/(OMEGA+CHOELE(JCHO))
         END DO
C
      END DO
C
      RETURN
      END
C
C  /* Deck ccho_sct2o */
      SUBROUTINE CCHO_SCT2O(T21,T23,T24,OCCHO,NUMII,ISYT2O,NTOTT1,
     &                      IOFF2,IOFF3)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C     Scale T2 subblocks with occupied part of Cholesky vector.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION T21(*),T23(*),T24(*),OCCHO(*)
      DIMENSION IOFF2(8),IOFF3(8)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
C
      DO ISYMA = 1,NSYM
C
         ISYMEM = MULD2H(ISYT2O,ISYMA)
C
         DO I = 1,NUMII
            DO A = 1,NVIR(ISYMA)
C
               NAI = NVIR(ISYMA)*(I-1) + A
C
               DO NEM = 1,NT1AM(ISYMEM)
C
                  NEMAI = IOFF3(ISYMA) + NT1AM(ISYMEM)*(NAI-1) + NEM
                  NEMI  = NTOTT1*(I-1) + IOFF2(ISYMEM) + NEM
C
                  T21(NEMAI) = T21(NEMAI)*OCCHO(NEMI)
                  T23(NEMAI) = T23(NEMAI)*OCCHO(NEMI)
                  T24(NEMAI) = T24(NEMAI)*OCCHO(NEMI)
C
               END DO
            END DO
         END DO
      END DO
C
      RETURN
      END
C
C  /* Deck ccho_descl */
      SUBROUTINE CCHO_DESCL(XINTI,XINTJ,VICHO,ISINT,B,IOFF2,NTOTT1)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C     Scale integrals with Cholesky vectors.
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XINTI(*),XINTJ(*),VICHO(*) 
      DIMENSION IOFF2(8)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
C
      DO ISYMA = 1,NSYM
C
         ISYMEM = MULD2H(ISINT,ISYMA)
C
         DO A = 1,NVIR(ISYMA)
            DO NEM = 1,NT1AM(ISYMEM)
C
               NEMA = ICKATR(ISYMEM,ISYMA) + NT1AM(ISYMEM)*(A-1) + NEM
               NEMB = NTOTT1*(B-1) + IOFF2(ISYMEM) + NEM
C
               XINTI(NEMA) = XINTI(NEMA)*VICHO(NEMB)
               XINTJ(NEMA) = XINTJ(NEMA)*VICHO(NEMB)
C
            END DO
         END DO
C
      END DO
C
      RETURN
      END
C  /* Deck ccho_eterm2 */
      SUBROUTINE CCHO_ETERM2(XIINT,XJINT,XOINT,T2VO,FOCKD,NUMCHO,CHOELE,
     *                       WORK,LWORK,E4DE,IB1,ISYMB,NUMIB,
     *                       FBATCH,PRINT,NONI)
C
C     Javier Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch and Alfredo Sanchez, June 2002
C
C
C     Calculate (part of E-terms):
C
C        E4DE = Sum(b) Sum(ijk) [ Q(kj,i;b) * K(jk,i;b)
C                               + R(kj,i;b) * K(kj,i;b) ]
C
C     where
C
C        Q(kj,i;b) = 2 Sum(ac) t(ac,kj) * S(ac,i;b)
C
C        R(kj,i;b) = 2 Sum(ac) t(ac,kj) * U(ac,i;b)
C
C        S(ac,i;b) =   Sum(dl) d(dl,a) * s(dl,ai) * (dl|bc)
C                    - Sum(dl) d(dl,a) * t(dl,ai) * (bl|dc)
C
C        U(ac,i;b) = - Sum(dl) d(dl,a) * t(di,al) * (bl|dc)
C
C     and s(ai,bj) = 2 t(ai,bj) - t(aj,bi), and d(ai,j) and d(ai,b) denotes
C
C     and s(ai,bj) = 2 t(ai,bj) - t(aj,bi), and d(ai,j) and d(ai,b) denotes
C     the occupied and virtual parts of the Cholesky decomposition of the
C     orbital energy denominator.
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XIINT(*),XJINT(*),XOINT(*),T2VO(*),CHOELE(*),FOCKD(*)
      DIMENSION WORK(LWORK)
      DIMENSION IOFF3(8)
      LOGICAL FBATCH,PRINT
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C 
      dimension scndsce(maxcho)
      TLAST = SECOND()
C 
C
      IF (PRINT) THEN
         TIMT = SECOND()
         WRITE(LUPRI,'(6X,A,/,6X,A,/)')
     &   'Calculation of the E2 term:',
     &   '==========================='
      ENDIF
C
      ISYMT2 = 1
C
C     I(em,d,b) = (bm|ed) and J(em,d,b) = (em|bd) 
C     K(jk,i,b) = 2*(jb,ki) - (kb,ji).
C     -------------------------------------------
C
      ISINTI = ISYMB
      ISINTJ = ISYMB
      ISINTK = ISYMB
C
      KEND1  = 1
      LWRK1  = LWORK
C
      DO 100 ISYMI = 1,NSYM
C
         IF (NRHF(ISYMI) .EQ. 0) GOTO 100
C
C        Symmetry of amplitudes subblocks
C        --------------------------------
C
         ISYT21 = ISYMI
         ISYT22 = ISYMI
         ISYT23 = ISYMI
C
         ISYMBI = MULD2H(ISYMI,ISYMB)
         ISYMAC = ISYMBI
         ISYMKJ = MULD2H(ISYMAC,ISYMT2)
C
C-------------------
C        Allocation.
C-------------------
C
         KCHOV = KEND1
         KCHOO = KCHOV + NCKATR(ISYT21)
         KEND2 = KCHOO + NMATIJ(ISYMKJ)*NUMIB
         LWRK2 = LWORK - KEND2 + 1
C
         IF (LWRK2 .LE. 0) THEN
            WRITE(LUPRI,*) 'Insufficient memory in CCHO_ETERM2'
            WRITE(LUPRI,*) 'Need (more than): ',KEND2-1
            WRITE(LUPRI,*) 'Available       : ',LWORK
            CALL QUIT('Insufficient memory in CCHO_ETERM2')
         ENDIF
C
C        Find largest subblock of AC matrix.
C        -----------------------------------
C
         MAXAC = -1
         DO ISYMC = 1,NSYM
            ISYMA = MULD2H(ISYMC,ISYMAC)
            MAXAC = MAX(MAXAC,NVIR(ISYMA)*NVIR(ISYMC))
         ENDDO
C
C-------------------------------
C        Batch over the I index.
C-------------------------------
C
         LENWRK = 3*NCKATR(ISYT21) + 2*NMATAB(ISYMAC) + MAXAC
     &          + 3*NMATIJ(ISYMKJ)
         IF (FBATCH) THEN
            NEFI  = MIN(NONI,NRHF(ISYMI))
            LEFF  = NEFI*LENWRK + 1
            LWRK2 = MIN(LWRK2,LEFF)
         ENDIF
         NUMI = MIN(NRHF(ISYMI),LWRK2/LENWRK)
         IF (NUMI .EQ. 0) THEN
            WRITE(LUPRI,*) 'NUMI .EQ. 0 in CC_ETERM2'
            CALL QUIT('Not enough space in CC_ETERM2')
         END IF
C
         NBATI = (NRHF(ISYMI)-1)/NUMI + 1
C
         IF (PRINT) THEN
            WRITE(LUPRI,'(9X,A,I1,A,/,9X,A)')
     &      'Batch over I, symmetry ',ISYMI,':',
     &      '-------------------------'
            WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10)')
     &      'Minimum work space required   : ',LENWRK,
     &      'Work space available for batch: ',LWRK2
            WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10,/)')
     &      'Number of occupied orbitals   : ',NRHF(ISYMI),
     &      'Required number of I-batches  : ',NBATI
         ENDIF
C
         II2 = 0
         DO IBATI = 1,NBATI
C
            II1 = II2 + 1
            II2 = II2 + NUMI
            IF (II2 .GT. NRHF(ISYMI)) II2 = NRHF(ISYMI)
            NUMII = II2 - II1 + 1
C
            IF (PRINT) THEN
               WRITE(LUPRI,'(12X,A,I10,A,/,12X,A)')
     &         'I-batch number ',IBATI,':',
     &         '--------------------------'
               WRITE(LUPRI,'(12X,A,I10,1X,I10,/)')
     &         'First and last I: ',II1,II2
            ENDIF
C
C-----------------------------
C           Some index arrays.
C-----------------------------
C
            ICOUN3 = 0
            DO ISYMA = 1,NSYM
C
               ISYMAI = MULD2H(ISYMA,ISYMI)
               ISYMDL = MULD2H(ISYMAI,ISYMT2)
C
               IOFF3(ISYMA) = ICOUN3
C
               ICOUN3 = ICOUN3 + NT1AM(ISYMDL)*NVIR(ISYMA)*NUMII
C
            END DO
C
C------------------------------
C           Dynamic allocation.
C------------------------------
C
            KT21  = KEND2
            KT22  = KT21  + NCKATR(ISYT21)*NUMII
            KT23  = KT22  + NCKATR(ISYT22)*NUMII
            KSMAT = KT23  + NCKATR(ISYT23)*NUMII
            KUMAT = KSMAT + NMATAB(ISYMAC)*NUMII
            KQMAT = KUMAT + NMATAB(ISYMAC)*NUMII
            KRMAT = KQMAT + NMATIJ(ISYMKJ)*NUMII
            KKINT = KRMAT + NMATIJ(ISYMKJ)*NUMII
            KSSCR = KKINT + NMATIJ(ISYMKJ)*NUMII
            KEND3 = KSSCR + MAXAC*NUMII
            LWRK3 = LWORK - KEND3 + 1
C
            IF (LWRK3 .LT. 0) THEN
               WRITE(LUPRI,*) 'Batching problem in CCHO_ETERM2 !!!!'
               CALL QUIT('Batching bug in CCHO_ETERM2 !!!!')
            ENDIF
C
C---------------------------------
C           Extract T2 amplitudes.
C---------------------------------
C
            CALL CCHO_DET22(T2VO,WORK(KT21),WORK(KT22),WORK(KT23),
     &                      II1,NUMII,ISYMI,IOFF3,ISYMT2)
C
            DO ICHO = 1,NUMCHO
C
               DECOR = ZERO
C
C---------------------------------------------
C              Calculate Cholesky information.
C---------------------------------------------
C
               CALL CCHO_DECHO2(FOCKD,CHOELE,NUMCHO,ICHO,
     &                          WORK(KCHOO),WORK(KCHOV),IB1,NUMIB,ISYMB,
     &                          ISYMKJ,ISYMI)
C
C-----------------------------------------------------------
C              Scale the T2 amplitudes with Cholesky update.
C-----------------------------------------------------------
C
               CALL CCHO_SCT2O2(WORK(KT21),WORK(KT22),WORK(KT23),
     &                          WORK(KCHOV),NUMII,ISYT21,IOFF3)
C
               DO B = 1,NUMIB
C
                  IB = IB1 + B - 1
C
C-------------------------------------------------
C                 Extract 2CME occupied integrals.
C-------------------------------------------------
C
                  CALL CCHO_DEINO(XOINT,WORK(KKINT),II1,NUMII,ISYMI,IB,
     &                            ISYMB)
C
C------------------------------------------
C                 Scale occupied integrals.
C------------------------------------------
C
                  KOFFC = KCHOO + NMATIJ(ISYMKJ)*(B - 1)
                  CALL CCHO_DESCL2(WORK(KKINT),WORK(KOFFC),NUMII,ISYMI,
     &                             ISINTK)
C
                  DO ISYMC = 1,NSYM
C
                     ISYMA  = MULD2H(ISYMC,ISYMAC)
                     ISYMDL = MULD2H(ISYMA,ISYT21)
C
C--------------------------------
C                    Calculate S.
C--------------------------------
C
                     KOFT1 = KT21 + IOFF3(ISYMA)
                     KOFFJ = NCKATR(ISYMB)*(B - 1)
     &                     + ICKATR(ISYMDL,ISYMC) + 1
C
                     NC    = NVIR(ISYMC)
                     NDL   = NT1AM(ISYMDL)
                     NTODL = MAX(NDL,1)
                     NAI   = NVIR(ISYMA)*NUMII
                     NTOAI = MAX(NAI,1)
C
                     CALL DGEMM('T','N',NAI,NC,NDL,
     &                          ONE,WORK(KOFT1),NTODL,
     &                          XJINT(KOFFJ),NTODL,
     &                          ZERO,WORK(KSSCR),NTOAI)
C
                     KOFT2 = KT22 + IOFF3(ISYMA)
                     KOFFI = KOFFJ
C
                     CALL DGEMM('T','N',NAI,NC,NDL,
     &                          XMONE,WORK(KOFT2),NTODL,
     &                          XIINT(KOFFI),NTODL,
     &                          ONE,WORK(KSSCR),NTOAI)
C
                     DO C = 1,NVIR(ISYMC)
                        DO I = 1,NUMII
                           KOFF1 = KSSCR + NVIR(ISYMA)*NUMII*(C - 1)
     &                           + NVIR(ISYMA)*(I - 1)
                           KOFF2 = KSMAT
     &                           + NMATAB(ISYMAC)*(I - 1)
     &                           + IMATAB(ISYMA,ISYMC)
     &                           + NVIR(ISYMA)*(C - 1)
                           CALL DCOPY(NVIR(ISYMA),WORK(KOFF1),1,
     &                                            WORK(KOFF2),1)
                        ENDDO
                     ENDDO
C
C--------------------------------
C                    Calculate U.
C--------------------------------
C
                     KOFT3 = KT23 + IOFF3(ISYMA)
C
                     CALL DGEMM('T','N',NAI,NC,NDL,
     &                          XMONE,WORK(KOFT3),NTODL,
     &                          XIINT(KOFFI),NTODL,
     &                          ZERO,WORK(KSSCR),NTOAI)
C
                     DO C = 1,NVIR(ISYMC)
                        DO I = 1,NUMII
                           KOFF1 = KSSCR + NVIR(ISYMA)*NUMII*(C - 1)
     &                           + NVIR(ISYMA)*(I - 1)
                           KOFF2 = KUMAT
     &                           + NMATAB(ISYMAC)*(I - 1)
     &                           + IMATAB(ISYMA,ISYMC)
     &                           + NVIR(ISYMA)*(C - 1)
                           CALL DCOPY(NVIR(ISYMA),WORK(KOFF1),1,
     &                                            WORK(KOFF2),1)
                        ENDDO
                     ENDDO
C
                  ENDDO
C
                  KOFFT = IT2VO(ISYMAC,ISYMKJ) + 1 
C
                  NKJ   = NMATIJ(ISYMKJ)
                  NTOKJ = MAX(1,NKJ)
                  NAC   = NMATAB(ISYMAC)
                  NTOAC = MAX(1,NAC)
C
                  CALL DGEMM('T','N',NKJ,NUMII,NAC,
     &                       ONE,T2VO(KOFFT),NTOAC,
     &                       WORK(KSMAT),NTOAC,
     &                       ZERO,WORK(KQMAT),NTOKJ)
C
                  CALL DGEMM('T','N',NKJ,NUMII,NAC,
     &                       ONE,T2VO(KOFFT),NTOAC,
     &                       WORK(KUMAT),NTOAC,
     &                       ZERO,WORK(KRMAT),NTOKJ)
C
                  DO I = 1,NUMII
                     DO ISYMJ = 1,NSYM
C
                        ISYMK = MULD2H(ISYMKJ,ISYMJ)
C
                        DO J = 1,NRHF(ISYMJ)
                           DO K = 1,NRHF(ISYMK)
C
                              NKJ  = IMATIJ(ISYMK,ISYMJ)
     &                             + NRHF(ISYMK)*(J-1) + K
                              NJK  = IMATIJ(ISYMJ,ISYMK)
     &                             + NRHF(ISYMJ)*(K-1) + J
C
                              NKJI = NMATIJ(ISYMKJ)*(I-1) + NKJ
                              NJKI = NMATIJ(ISYMKJ)*(I-1) + NJK
C
                              KOFFQ = KQMAT + NKJI -1
                              KOFFR = KRMAT + NKJI -1
                              KOFF1 = KKINT + NJKI -1
                              KOFF2 = KKINT + NKJI -1
C
                              DECOR = DECOR 
     &                              + WORK(KOFFQ)*WORK(KOFF1)
     &                              + WORK(KOFFR)*WORK(KOFF2)
C
                           ENDDO  
                        ENDDO
C
                     ENDDO
                  ENDDO
C
               ENDDO
C
               DECOR = TWO*DECOR
               E4DE  = E4DE + DECOR 
C 
               ENERGE(ICHO) = ENERGE(ICHO) + DECOR
C 
C Decommented by Domenico
              TNOW = SECOND()
              DELTAT = TNOW - TLAST
              TLAST = TNOW
              SCNDSE(ICHO)= SCNDSE(ICHO) + DELTAT
C Decommented by Domenico
C 
C
c              WRITE(LUPRI,'(A,4I4,2D15.6)')'ISC,ISI,IBI,IV,DCO,DT',
c    &                  ISYMC,ISYMI,IBATI,ICHO,DECOR,E4DT
C
c              IF (PRINT) THEN
c                 WRITE(LUPRI,'(15X,A,I3,A,/,15X,A)')
c    &            'Status after Cholesky vector',ICHO,':',
c    &            '--------------------------------'
c                 IF (ABS(DECOR) .LT. THRCHO) THEN
c                    WRITE(LUPRI,'(15X,A)') 'E2 term converged'
c                 ELSE
c                    WRITE(LUPRI,'(15X,A)') 'E2 term not converged'
c                 ENDIF
c                 TIM = SECOND() - TIMT
c                 WRITE(LUPRI,'(15X,A,F10.2,A,/)')
c    &            'Accumulated E2-time: ',TIM,' seconds'
c              ENDIF
C
               IF (ABS(DECOR) .LT. THRCHO) GOTO 999 
C
            ENDDO 
C
  999       CONTINUE
C
         ENDDO
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccho_det22 */
      SUBROUTINE CCHO_DET22(T2VO,T21,T22,T23,II1,NUMII,ISYMI,IOFF3,
     &                      ISYMT2)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (TWO = 2.0D0)
C
      DIMENSION T2VO(*),T21(*),T22(*), T23(*)
      DIMENSION IOFF3(8)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DO I = 1,NUMII
C
         II = II1 + I - 1
C
         DO ISYMC = 1,NSYM
C
            ISYMCI = MULD2H(ISYMC,ISYMI)
            ISYMEM = MULD2H(ISYMCI,ISYMT2)
C
            DO C = 1,NVIR(ISYMC)
               DO ISYMM = 1,NSYM
C
                  ISYME  = MULD2H(ISYMEM,ISYMM)
                  ISYMEC = MULD2H(ISYME,ISYMC)
                  ISYMMI = MULD2H(ISYMM,ISYMI)
C
                  DO M = 1,NRHF(ISYMM)
                     DO E = 1,NVIR(ISYME)
C
                        NEC = IMATAB(ISYME,ISYMC) 
     &                      + NVIR(ISYME)*(C-1) + E
                        NCE = IMATAB(ISYMC,ISYME)
     &                      + NVIR(ISYMC)*(E-1) + C
                        NMI = IMATIJ(ISYMM,ISYMI)
     &                      + NRHF(ISYMM)*(II-1) + M
C
                        NCI = NVIR(ISYMC)*(I-1) + C
                        NEM = IT1AM(ISYME,ISYMM) + NVIR(ISYME)*(M-1) + E
C
                        NECMI = IT2VO(ISYMEC,ISYMMI)
     &                        + NMATAB(ISYMEC)*(NMI-1) + NEC
                        NCEMI = IT2VO(ISYMEC,ISYMMI)
     &                        + NMATAB(ISYMEC)*(NMI-1) + NCE
C
                        NEMCI = IOFF3(ISYMC) + NT1AM(ISYMEM)*(NCI-1)
     &                        + NEM
C
                        T21(NEMCI) = TWO*T2VO(NECMI) - T2VO(NCEMI)
                        T22(NEMCI) = T2VO(NECMI)
                        T23(NEMCI) = T2VO(NCEMI)
C
                     END DO
                  END DO
C
               END DO
            END DO
C
         END DO
C
      END DO
C
C
      RETURN
      END
C  /* Deck ccho_sct2o2 */
      SUBROUTINE CCHO_SCT2O2(T21,T22,T23,VCCHO,NUMII,ISYT2O,IOFF3)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C     Scale T2 subblocks with virtual part of Cholesky vector.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION T21(*),T22(*),T23(*),VCCHO(*)
      DIMENSION IOFF3(8)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DO ISYMA = 1,NSYM
C
         ISYMDL = MULD2H(ISYT2O,ISYMA)
C
         DO I = 1,NUMII
            DO A = 1,NVIR(ISYMA)
C
               NAI = NVIR(ISYMA)*(I-1) + A
C
               DO NDL = 1,NT1AM(ISYMDL)
C
                  NDLAI = IOFF3(ISYMA) + NT1AM(ISYMDL)*(NAI-1) + NDL
                  NDLA  = ICKATR(ISYMDL,ISYMA) + NT1AM(ISYMDL)*(A - 1)
     &                  + NDL
C
                  T21(NDLAI) = T21(NDLAI)*VCCHO(NDLA)
                  T22(NDLAI) = T22(NDLAI)*VCCHO(NDLA)
                  T23(NDLAI) = T23(NDLAI)*VCCHO(NDLA)
C
               END DO
C
            END DO
         END DO
C
      END DO
C
      RETURN
      END
C  /* Deck ccho_deino */
      SUBROUTINE CCHO_DEINO(XOINT,XKINT,I1,NUMI,ISYMI,B,ISYMB)
C
C     JLC, BFR, TBP, HK AND AS, June 2002.
C
C     Extract 2CME of (bj|ik) integrals (ordered as (jk,i,b)).
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XOINT(*), XKINT(*)
#include <ccorb.h>
#include <ccsdsym.h>

      PARAMETER (TWO = 2.00D0)

      ISYMBI = MULD2H(ISYMB,ISYMI)
      ISYMJK = ISYMBI
      ISYJKI = ISYMB

      DO II = 1,NUMI

         I = I1 + II - 1

         DO ISYMK = 1,NSYM

            ISYMJ = MULD2H(ISYMK,ISYMJK)

            DO K = 1,NRHF(ISYMK)
               DO J = 1,NRHF(ISYMJ)

                  JK = IMATIJ(ISYMJ,ISYMK) + NRHF(ISYMJ)*(K - 1) + J
                  KJ = IMATIJ(ISYMK,ISYMJ) + NRHF(ISYMK)*(J - 1) + K

                  NJKIB = ISJIKA(ISYJKI,ISYMB)
     &                  + NMAJIK(ISYJKI)*(B - 1)
     &                  + ISJIK(ISYMJK,ISYMI)
     &                  + NMATIJ(ISYMJK)*(I - 1) + JK
                  NKJIB = ISJIKA(ISYJKI,ISYMB)
     &                  + NMAJIK(ISYJKI)*(B - 1)
     &                  + ISJIK(ISYMJK,ISYMI)
     &                  + NMATIJ(ISYMJK)*(I - 1) + KJ

                  NJKI  = NMATIJ(ISYMJK)*(II - 1) + JK

                  XKINT(NJKI) = TWO*XOINT(NJKIB) - XOINT(NKJIB)

               ENDDO
            ENDDO

         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck ccho_descl2 */
      SUBROUTINE CCHO_DESCL2(XKINT,CHOO,NUMI,ISYMI,ISINTK)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Scale occupied integrals.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XKINT(*), CHOO(*)
#include <ccorb.h>
#include <ccsdsym.h>

      ISYMJK = MULD2H(ISYMI,ISINTK)

      DO I = 1,NUMI
         DO JK = 1,NMATIJ(ISYMJK)

            JKI = NMATIJ(ISYMJK)*(I - 1) + JK

            XKINT(JKI) = XKINT(JKI)*CHOO(JK)

         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck ccho_decho2 */
      SUBROUTINE CCHO_DECHO2(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                      IB1,NUMBB,ISYMB,ISYMJK,ISYDLA)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Construct vector for Cholesky decomposition:
C
C     For the occupied part, OCCHO(jk,#b)  actual vector
C     For the virtual  part, VICHO(emc)    update vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C-------------------
C     Occupied part.
C-------------------
C
      IND = 0
      DO B = 1,NUMBB
         DO ISYMK = 1,NSYM
            ISYMJ = MULD2H(ISYMK,ISYMJK)
            DO K = 1,NRHF(ISYMK)
               DO J = 1,NRHF(ISYMJ)
C
                  KOFFK = IRHF(ISYMK) + K
                  KOFFJ = IRHF(ISYMJ) + J
                  KOFFB = IVIR(ISYMB) + IB1 + B -1
C
                  IND = IND + 1
                  OCCHO(IND) = FOCKD(KOFFB)-FOCKD(KOFFJ)-FOCKD(KOFFK)
C
               END DO
            END DO
         END DO
      END DO
C
      NDIMOC = IND
C
C
C     Construct the vector.
C     ---------------------
C
      DO P = 1,NDIMOC
C
         OMEGA = OCCHO(P)
C
C
         OCCHO(P) = SQRT(TWO*CHOELE(ICHO))/(OMEGA+CHOELE(ICHO))
C
         DO JCHO = 1,ICHO-1
            OCCHO(P) = OCCHO(P)
     &               *(OMEGA-CHOELE(JCHO))/(OMEGA+CHOELE(JCHO))
         END DO
C
      END DO
C
C------------------
C     Virtual part.
C------------------
C
      IF (ICHO .EQ. 1) THEN
C
C        Get the omega's.
C        ----------------
C
         IND = 0
         DO ISYMC = 1,NSYM
            ISYMEM = MULD2H(ISYMC,ISYDLA)
            DO C = 1,NVIR(ISYMC)
               KOFFC = IVIR(ISYMC) + C
               DO ISYMM = 1,NSYM
                  ISYME = MULD2H(ISYMEM,ISYMM)
                  DO M = 1,NRHF(ISYMM)
                     KOFFM = IRHF(ISYMM) + M
                     DO E = 1,NVIR(ISYME)
                        KOFFE = IVIR(ISYME) + E
C
                        OME = FOCKD(KOFFE)-FOCKD(KOFFM)+FOCKD(KOFFC)
C
                        IND = IND + 1
                        VICHO(IND) =
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
      ELSE
C
         IND = 0
         DO ISYMC = 1,NSYM
            ISYMEM = MULD2H(ISYMC,ISYDLA)
            DO C = 1,NVIR(ISYMC)
               KOFFC = IVIR(ISYMC) + C
               DO ISYMM = 1,NSYM
                  ISYME = MULD2H(ISYMEM,ISYMM)
                  DO M = 1,NRHF(ISYMM)
                     KOFFM = IRHF(ISYMM) + M
                     DO E = 1,NVIR(ISYME)
                        KOFFE = IVIR(ISYME) + E
C
                        IND = IND + 1
                        OME = FOCKD(KOFFC)+FOCKD(KOFFE)-FOCKD(KOFFM)
C
                        VICHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
         NDIMVI = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMVI,FACTOR,VICHO,1)
C
      ENDIF
C
      RETURN
      END
C  /* Deck ccho_dterm */
      SUBROUTINE CCHO_DTERM(XJINT,XKINT,T2VO,FOCKD,NUMCHO,CHOELE,
     *                      WORK,LWORK,E4D,IA1,ISYMA,NUMIA,
     *                      FBATCH,PRINT,NONI)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Calculate D-term:
C
C        E4D = E4D + 2 * Sum(a) Sum(i) [ V(a,i) * W(a,i) ]
C
C     where
C
C        V(a,i) = Sum(cjb) d(cjb) * (jc|ab) * s(cj,bi)
C
C        W(a,i) = Sum(dlk) d(dlk) * s(dl,ak) * (dl|ik)
C
C     and s(ai,bj) = 2 t(ai,bj) - t(aj,bi), and d(aij) and d(aib) denote
C     the occupied and virtual parts of the Cholesky decomposition of the
C     orbital energy denominator.
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XJINT(*),XKINT(*),T2VO(*),CHOELE(*),FOCKD(*)
      DIMENSION WORK(LWORK)
      DIMENSION IOFF3(8)
      LOGICAL FBATCH,PRINT
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C 
      TLAST = SECOND()
C 
C
      IF (PRINT) TIMT = SECOND()
C
      ISYDLK = ISYMA
      ISYCJB = ISYMA
      ISYMI  = ISYMA
C
C-----------------------------
C     Return if nothing to do.
C-----------------------------
C
      IF ((NRHF(ISYMI).LE.0) .OR. (NCKATR(ISYCJB).LE.0) .OR.
     &    (NCKI(ISYDLK).LE.0) .OR. (NUMIA.LE.0)) RETURN
C
      IF (PRINT) THEN
         WRITE(LUPRI,'(6X,A,/,6X,A,/)')
     &   'Calculation of the D term:',
     &   '=========================='
      ENDIF
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KT21  = 1
      KCHOO = KT21  + NCKI(ISYDLK)*NUMIA 
      KCHOV = KCHOO + NCKI(ISYDLK)
      KEND1 = KCHOV + NCKATR(ISYCJB)
      LWRK1 = LWORK - KEND1 + 1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in CCHO_DTERM')
      ENDIF
C
C----------------------
C     Extract s(dlk,a).
C----------------------
C
      CALL CCHO_TDLKA(T2VO,WORK(KT21),IA1,NUMIA,ISYMA)
C
C------------------------
C     Batch over i index.
C------------------------
C
      LENWRK = NCKATR(ISYCJB) + NCKI(ISYDLK)
     &       + 2*NUMIA
      IF (FBATCH) THEN
         NEFI  = MIN(NRHF(ISYMI),NONI)
         LEFF  = NEFI*LENWRK + 1
         LWRK1 = MIN(LWRK1,LEFF)
      ENDIF
      NUMI = MIN(NRHF(ISYMI),LWRK1/LENWRK)
      IF (NUMI .EQ. 0) THEN
         WRITE(LUPRI,*) 'NUMI .EQ. 0 in CCHO_DTERM'
         CALL QUIT('Not enough space in CCHO_DTERM')
      END IF
C
      NBATI = (NRHF(ISYMI)-1)/NUMI + 1
C
      IF (PRINT) THEN
         WRITE(LUPRI,'(6X,A,I1,A,/,6X,A)')
     &   'Batch over I, symmetry ',ISYMI,':',
     &   '-------------------------'
         WRITE(LUPRI,'(6X,A,I10,/,6X,A,I10)')
     &   'Minimum work space required   : ',LENWRK,
     &   'Work space available for batch: ',LWRK1
         WRITE(LUPRI,'(6X,A,I10,/,6X,A,I10,/)')
     &   'Number of occupied orbitals   : ',NRHF(ISYMI),
     &   'Required number of I-batches  : ',NBATI
      ENDIF
C
      II2 = 0
      DO IBATI = 1,NBATI
C
         II1 = II2 + 1
         II2 = II2 + NUMI
         IF (II2 .GT. NRHF(ISYMI)) II2 = NRHF(ISYMI)
         NUMII = II2 - II1 + 1
C
         IF (PRINT) THEN
            WRITE(LUPRI,'(9X,A,I10,A,/,9X,A)')
     &      'I-batch number ',IBATI,':',
     &      '--------------------------'
            WRITE(LUPRI,'(9X,A,I10,1X,I10,/)')
     &      'First and last I: ',II1,II2
         ENDIF
C
C---------------------------
C        Dynamic allocation.
C---------------------------
C
         KT22  = KEND1
         KKINT = KT22  + NCKATR(ISYCJB)*NUMII
         KVMAT = KKINT + NCKI(ISYDLK)*NUMII
         KWMAT = KVMAT + NUMIA*NUMII
         KEND2 = KWMAT + NUMIA*NUMII
         LWRK2 = LWORK - KEND2 + 1
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Batching bug in CCHO_DTERM !!!!')
         ENDIF
C
C-------------------------
C        Extract s(cjb,i).
C-------------------------
C
         CALL CCHO_TCJBI(T2VO,WORK(KT22),II1,NUMII,ISYMI)
C
C-------------------------
C        Extract K(dlk,i).
C-------------------------
C
         CALL CCHO_KDLKI(XKINT,WORK(KKINT),II1,NUMII,ISYMI)        
C
C----------------------
C        Cholesky loop.
C----------------------
C
         DO ICHO = 1,NUMCHO
C
C---------------------------------------
C           Get Cholesky update vectors.
C---------------------------------------
C
            CALL CCHO_DECHO3(FOCKD,CHOELE,NUMCHO,ICHO,WORK(KCHOO),
     &                       WORK(KCHOV),ISYDLK,ISYCJB,1,1)

c     ocnorm = dsqrt(ddot(NCKI(ISYDLK),WORK(KCHOO),1,WORK(KCHOO),1))
c     write(LUPRI,*) '   DTERM:'
c     write(LUPRI,*) '          ICHO,ISYDLK           : ',ICHO,ISYDLK
c     write(LUPRI,*) '          Norm of d(dl,k) update: ',ocnorm
C
C---------------------------------------------------
C           Scale amplitudes and occupied integrals.
C---------------------------------------------------
C
            CALL CCHO_SCVEC(WORK(KT22),WORK(KCHOV),NCKATR(ISYCJB),NUMII)
            CALL CCHO_SCVEC(WORK(KKINT),WORK(KCHOO),NCKI(ISYDLK),NUMII)
C
C-----------------------------
C           Calculate V and W.
C-----------------------------
C
            CALL DGEMM('T','N',NUMIA,NUMII,NCKATR(ISYCJB),
     &                 ONE,XJINT,NCKATR(ISYCJB),
     &                 WORK(KT22),NCKATR(ISYCJB),
     &                 ZERO,WORK(KVMAT),NUMIA)
C
            CALL DGEMM('T','N',NUMIA,NUMII,NCKI(ISYDLK),
     &                 ONE,WORK(KT21),NCKI(ISYDLK),
     &                 WORK(KKINT),NCKI(ISYDLK),
     &                 ZERO,WORK(KWMAT),NUMIA)
C
C-----------------------------------------
C           Calculate energy contribution.
C-----------------------------------------
C
            EDCOR = TWO*DDOT(NUMIA*NUMII,WORK(KVMAT),1,WORK(KWMAT),1)
            E4D   = E4D + EDCOR
C 
              ENERGD(ICHO) = ENERGD(ICHO) + EDCOR
C 
C Decommented by Domenico
           TNOW = SECOND()
           DELTAT = TNOW - TLAST
           TLAST = TNOW
            SCNDSD(ICHO) = SCNDSD(ICHO) 
     &                          + DELTAT
C 
C Decommented by Domenico
C
c     write(LUPRI,*) '   DTERM: E4D = ',E4D,' (ICHO = ',ICHO,')'
C
c           IF (PRINT) THEN
c              WRITE(LUPRI,'(12X,A,I3,A,/,12X,A)')
c    &         'Status after Cholesky vector',ICHO,':',
c    &         '--------------------------------'
c              IF (ABS(EDCOR) .LT. THRCHO) THEN
c                 WRITE(LUPRI,'(12X,A)') 'D term converged'
c              ELSE
c                 WRITE(LUPRI,'(12X,A)') 'D term not converged'
c              ENDIF
c              TIM = SECOND() - TIMT
c              WRITE(LUPRI,'(12X,A,F10.2,A,/)')
c    &         'Accumulated D-time: ',TIM,' seconds'
c           ENDIF
C
            IF (ABS(EDCOR) .LT. THRCHO) GOTO 999
C
         ENDDO
C
  999    CONTINUE
C
      ENDDO
C
      RETURN
      END
C
C  /* Deck ccho_a1term */
c     SUBROUTINE CCHO_A1TERM(XJINT,T2VO,FOCKD,NUMCHO,CHOELE,
c    *                      WORK,LWORK,E4A1,THRCHO)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Calculate A1-term:
C
C        E4A1 = E4A1 + 2 Sum(i) Sum(jkl) [ W(kl,j) * V(kl,j) ]
C
C     where REMEMBER to CHANGE !!!!!!!!!!!!!!!!!!!!!!!!!!
C
C        V(a,i) = Sum(cjb) d(cjb) * (jc|ab) * s(cj,bi)
C
C        W(a,i) = Sum(dlk) d(dlk) * s(dl,ak) * (dl|ik)
C
C     and s(ai,bj) = 2 t(ai,bj) - t(aj,bi), and d(aij) and d(aib) denote
C     the occupied and virtual parts of the Cholesky decomposition of the
C     orbital energy denominator.
C
c#include <implicit.h>
c#include <priunit.h>
C
c     DIMENSION XJINT(*),T2VO(*),CHOELE(*),FOCKD(*)
c     DIMENSION WORK(LWORK)
c     DIMENSION IOFF3(8)
C
c#include <ccorb.h>
c#include <ccsdsym.h>
C
c     PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
C-------------------------------------
C     Start loop over symmetries of I.
C-------------------------------------
C
c     DO ISYMI = 1,NSYM
C
c        ISYABJ = ISYMI
C
C------------------------
C        Dynamic allocation.
C------------------------
C
c        KT21  = 1
c        KCHOO = KT21  + NCKI(ISYDLK)
c        KCHOV = KCHOO + NCKATR(ISYCJB)
c        KEND1 = KCHOV + NCKATR(ISYDLK)
c        LWRK1 = LWORK - KEND1 + 1
c        IF (LWRK1 .LT. 0) THEN
c           CALL QUIT('Insufficient memory in CCHO_DTERM')
c        ENDIF
C
C----------------------
C        Split work space.
C----------------------
C
c        LWRKIJ = LWORK/2
C
C---------------------------
C        Batch over i index.
C---------------------------
C
c        LEWRKI = NCKATR(ISYCJB) + NCKI(ISYDLK)
c    &          + 2*NUMIA
C
c        NUMI = MIN(NRHF(ISYMI),LWRKIJ/LEWRKI)
c        IF (NUMI .EQ. 0) THEN
c           WRITE(LUPRI,*) 'NUMI .EQ. 0 in CC_A1TERM'
c           CALL QUIT('Not enough space in CC_A1TERM')
c        END IF
C
c        NBATI = (NRHF(ISYMI)-1)/NUMI + 1
C
c        II2 = 0
c        DO IBATI = 1,NBATI
C
c           II1 = II2 + 1
c           II2 = II2 + NUMI
c           IF (II2 .GT. NRHF(ISYMI)) II2 = NRHF(ISYMI)
c           NUMII = II2 - II1 + 1
C
C------------------------------
C           Dynamic allocation.
C------------------------------
C
c           KT22  = KEND1
c           KKINT = KT22  + NCKATR(ISYCJB)*NUMII
c           KVMAT = KKINT + NCKI(ISYDLK)*NUMII
c           KWMAT = KVMAT + NUMIA*NUMII
c           KEND2 = KWMAT + NUMIA*NUMII
c           LWRK2 = LWORK - KEND2 + 1
c           IF (LWRK2 .LT. 0) THEN
c              CALL QUIT('Batching problem in CCHO_A1TERM')
c           ENDIF
C
C-------------------------------------------
C           Start loop over symmetries of J.
C-------------------------------------------
C
c           DO ISYMJ = 1,NSYM
C
C---------------------------------
C              Batch over index j.
C---------------------------------
C
c              LEWRKJ = NCKATR(ISYCJB) + NCKI(ISYDLK)
c    &                + 2*NUMIA
C
c              NUMJ = MIN(NRHF(ISYMJ),LWRK2/LEWRKJ)
c              IF (NUMJ .EQ. 0) THEN
c                 WRITE(LUPRI,*) 'NUMJ .EQ. 0 in CC_A1TERM'
c                 CALL QUIT('Not enough space in CC_A1TERM')
c              END IF
C
c              NBATJ = (NRHF(ISYMJ)-1)/NUMJ + 1
C
c              JJ2 = 0
c              DO IBATJ = 1,NBATJ
C
c                 JJ1 = JJ2 + 1
c                 JJ2 = JJ2 + NUMJ
c                 IF (JJ2 .GT. NRHF(ISYMJ)) JJ2 = NRHF(ISYMJ)
c                 NUMJJ = JJ2 - JJ1 + 1
C
C------------------------------------
C                 Dynamic allocation.
C------------------------------------
C
C
C-------------------------------
C                 Cholesky loop.
C-------------------------------
C
c                 DO ICHO = 1,NUMCHO
C
C---------------------------------------
C                    Get Cholesky update vectors.
C---------------------------------------
C
c                    CALL CCHO_DECHO3(FOCKD,CHOELE,NUMCHO,ICHO,
c    &                   WORK(KCHOO),WORK(KCHOV),ISYCJB,ISYDLK,1,1)
C
C------------------------------------------------------------
C                    Scale amplitudes and occupied integrals.
C------------------------------------------------------------
C
c           CALL CCHO_SCVEC(WORK(KT22),WORK(KCHOV),NCKATR(ISYCJB),NUMII)
c           CALL CCHO_SCVEC(WORK(KKINT),WORK(KCHOO),NCKI(ISYDLK),NUMII)
C
C-----------------------------
C           Calculate V and W.
C-----------------------------
C
c           CALL DGEMM('T','N',NUMIA,NUMII,NCKATR(ISYCJB),
c    &                 ONE,XJINT,NCKATR(ISYCJB),
c    &                 WORK(KT22),NCKATR(ISYCJB),
c    &                 ZERO,WORK(KVMAT),NUMIA)
C
c           CALL DGEMM('T','N',NUMIA,NUMII,NCKI(ISYDLK),
c    &                 ONE,WORK(KT21),NCKI(ISYDLK),
c    &                 WORK(KKINT),NCKI(ISYDLK),
c    &                 ZERO,WORK(KWMAT),NUMIA)
C
C-----------------------------------------
C           Calculate energy contribution.
C-----------------------------------------
C
c           EA1COR = TWO*DDOT(NUMIA*NUMII,WORK(KVMAT),1,WORK(KWMAT),1)
c           E4A1   = E4D + EA1COR
C
c           IF (ABS(EA1COR) .LT. THRCHO) GOTO 999
C
c        ENDDO
C
c 999    CONTINUE
C
c     ENDDO
C
c     RETURN
c     END
C
C  /* Deck ccho_scvec */
      SUBROUTINE CCHO_SCVEC(XMAT,SCAL,NI,NJ)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Scale XMAT.
C
#include <implicit.h>
C
      DIMENSION XMAT(NI,NJ), SCAL(NI)
C
      DO J = 1,NJ
         DO I = 1,NI
            XMAT(I,J) = XMAT(I,J)*SCAL(I)
         ENDDO
      ENDDO
C
      RETURN
      END
C
C  /* Deck ccho_kdlki */
      SUBROUTINE CCHO_KDLKI(XKINT,XKDLKI,II1,NUMII,ISYMI)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract K(lk,i,d) as K(dlk,i)
C
#include <implicit.h>
C
      PARAMETER (TWO = 2.0D0)
C
      DIMENSION XKINT(*),XKDLKI(*)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYDLK = ISYMI
      DO ISYMD = 1,NSYM
C
         ISYMLK = MULD2H(ISYMD,ISYDLK)
         ISYLKI = ISYMD
C
         DO D = 1,NVIR(ISYMD)
            DO I = 1,NUMII
C
               II = II1 + I - 1
C
               DO ISYMK = 1,NSYM
C
                  ISYML  = MULD2H(ISYMK,ISYMLK)
                  ISYMDL = MULD2H(ISYML,ISYMD)
C
                  DO K = 1,NRHF(ISYMK)
                     DO L = 1,NRHF(ISYML)
C
                        NLK   = IMATIJ(ISYML,ISYMK)
     &                        + NRHF(ISYML)*(K-1) + L
                        NDL   = IT1AM(ISYMD,ISYML)
     &                        + NVIR(ISYMD)*(L-1) + D
C
                        NLKI  = ISJIK(ISYMLK,ISYMI)
     &                        + NMATIJ(ISYMLK)*(II-1) + NLK
                        NDLK  = ICKI(ISYMDL,ISYMK)
     &                        + NT1AM(ISYMDL)*(K-1) + NDL
C
                        NLKID = ISJIKA(ISYLKI,ISYMD)
     &                        + NMAJIK(ISYLKI)*(D-1) + NLKI
                        NDLKI = NCKI(ISYDLK)*(I-1) + NDLK
C
                        XKDLKI(NDLKI) = XKINT(NLKID)
C
                     ENDDO
                  ENDDO
C
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO
 
      RETURN
      END
C  /* Deck ccho_tdlka */
      SUBROUTINE CCHO_TDLKA(T2VO,TDLKA,IA1,NA,ISYMA)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract T2 amplitudes as s(dlk,a)
C     s(dlk,a) = 2 t(da,lk) - t(da,kl)
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (XMONE = -1.0D0, TWO = 2.0D0)
C
      DIMENSION T2VO(*),TDLKA(*)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYDLK = ISYMA
C
      DO ISYMLK = 1,NSYM
         ISYMDA = ISYMLK
         ISYMD  = MULD2H(ISYMDA,ISYMA)
         DO ISYMK = 1,NSYM
            ISYML = MULD2H(ISYMLK,ISYMK)
            ISYMDL = MULD2H(ISYMD,ISYML)
            DO K = 1,NRHF(ISYMK)
               DO L = 1,NRHF(ISYML)
                  NLK = IMATIJ(ISYML,ISYMK)
     &                + NRHF(ISYML)*(K-1) + L
                  NKL = IMATIJ(ISYMK,ISYML)
     &                + NRHF(ISYMK)*(L-1) + K 
                  DO A = 1,NA
                     IAA   = IA1 + A - 1
                     KDALK = IT2VO(ISYMDA,ISYMLK)
     &                     + NMATAB(ISYMDA)*(NLK-1)
     &                     + IMATAB(ISYMD,ISYMA)
     &                     + NVIR(ISYMD)*(IAA-1) + 1
                     KDAKL = IT2VO(ISYMDA,ISYMLK)
     &                     + NMATAB(ISYMDA)*(NKL-1)
     &                     + IMATAB(ISYMD,ISYMA)
     &                     + NVIR(ISYMD)*(IAA-1) + 1
                     KDLKA = NCKI(ISYDLK)*(A-1)
     &                     + ICKI(ISYMDL,ISYMK)
     &                     + NT1AM(ISYMDL)*(K-1)
     &                     + IT1AM(ISYMD,ISYML)
     &                     + NVIR(ISYMD)*(L-1) + 1
                     CALL DCOPY(NVIR(ISYMD),T2VO(KDALK),1,
     &                                      TDLKA(KDLKA),1)
                     CALL DSCAL(NVIR(ISYMD),TWO,TDLKA(KDLKA),1)
                     CALL DAXPY(NVIR(ISYMD),XMONE,T2VO(KDAKL),1,
     &                                            TDLKA(KDLKA),1)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C
C  /* Deck ccho_tcjbi */
      SUBROUTINE CCHO_TCJBI(T2VO,TCJBI,II1,NUMII,ISYMI)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract T2 amplitudes as s(cjb,i)
C     s(cjb,i) = 2 t(cb,ji) - t(bc,ji)
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (XMONE = -1.0D0, TWO = 2.0D0)
C
      DIMENSION T2VO(*),TCJBI(*)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYCJB = ISYMI
C
      DO ISYMJI = 1,NSYM
C
         ISYMJ  = MULD2H(ISYMJI,ISYMI)
         ISYMCB = ISYMJI
         ISYMIJ = ISYMJI
C
         DO I = 1,NUMII
C
            II = II1 + I - 1
C
            DO J = 1,NRHF(ISYMJ)
C
               JI = IMATIJ(ISYMJ,ISYMI)
     &            + NRHF(ISYMJ)*(II - 1) + J
               IJ = IMATIJ(ISYMI,ISYMJ)
     &            + NRHF(ISYMI)*(J - 1)  + II
C
               DO ISYMB = 1,NSYM
C
                  ISYMC  = MULD2H(ISYMB,ISYMCB)
                  ISYMCJ = MULD2H(ISYMB,ISYCJB)
C
                  DO B = 1,NVIR(ISYMB)
C
                     KCBJI = IT2VO(ISYMCB,ISYMJI)
     &                     + NMATAB(ISYMCB)*(JI - 1)
     &                     + IMATAB(ISYMC,ISYMB)
     &                     + NVIR(ISYMC)*(B - 1) + 1
                     KCBIJ = IT2VO(ISYMCB,ISYMIJ)
     &                     + NMATAB(ISYMCB)*(IJ - 1)
     &                     + IMATAB(ISYMC,ISYMB) 
     &                     + NVIR(ISYMC)*(B - 1) + 1
                     KCJBI = NCKATR(ISYCJB)*(I - 1)
     &                     + ICKATR(ISYMCJ,ISYMB)
     &                     + NT1AM(ISYMCJ)*(B - 1)
     &                     + IT1AM(ISYMC,ISYMJ)
     &                     + NVIR(ISYMC)*(J - 1) + 1
C
                     CALL DCOPY(NVIR(ISYMC),T2VO(KCBJI),1,
     &                                      TCJBI(KCJBI),1)
                     CALL DSCAL(NVIR(ISYMC),TWO,TCJBI(KCJBI),1)
                     CALL DAXPY(NVIR(ISYMC),XMONE,T2VO(KCBIJ),1,
     &                                            TCJBI(KCJBI),1)
C
                  ENDDO
C
               ENDDO
C
            ENDDO
C
         ENDDO
C
      ENDDO
C old:
c     DO ISYMJI = 1,NSYM
c        ISYMJ = MULD2H(ISYMJI,ISYMI)
c        ISYMCB = ISYMJI
c        DO I = 1,NUMII
c           II = II1 + I - 1
c           DO J = 1,NRHF(ISYMJ)
c              NJI = IMATIJ(ISYMJ,ISYMI)
c    &             + NRHF(ISYMJ)*(II-1)
c    &             + J
c              DO ISYMB = 1,NSYM
c                 ISYMC = MULD2H(ISYMCB,ISYMB)
c                 ISYMCJ = MULD2H(ISYMC,ISYMJ)
c                 DO B = 1,NVIR(ISYMB)
c                    DO C = 1,NVIR(ISYMC)
c                       NCB = IMATAB(ISYMC,ISYMB)
c    &                      + NVIR(ISYMC)*(B-1) + C
c                       NBC = IMATAB(ISYMB,ISYMC)
c    &                      + NVIR(ISYMB)*(C-1) + B
c                       NCJ = IT1AM(ISYMC,ISYMJ) 
c    &                      + NVIR(ISYMC)*(J-1) + C
c                       NCJB = ICKATR(ISYMCJ,ISYMB)
c    &                       + NT1AM(ISYMCJ)*(B-1)
c    &                       + NCJ
c                       NCJBI = NCKATR(ISYCJB)*(I-1)
c    &                        + NCJB
c                       NCBJI = IT2VO(ISYMCB,ISYMJI)
c    &                        + NMATAB(ISYMCB)*(NJI-1)
c    &                        + NCB
c                       NBCJI = IT2VO(ISYMCB,ISYMJI)
c    &                        + NMATAB(ISYMCB)*(NJI-1)
c    &                        + NBC
c                       TCJBI(NCJBI) = TWO*T2VO(NCBJI)
c    &                               - T2VO(NBCJI)
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO 
c     ENDDO
C
      RETURN
      END
C
C  /* Deck ccho_decho4 */
      SUBROUTINE CCHO_DECHO4(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                      ISYCJB,ISYDLK)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C
C     Construct vector for Cholesky decomposition:
C
C     For the occupied part, OCCHO(dlk)  actual vector
C     For the virtual  part, VICHO(cjb)  actual vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C-------------------
C     Occupied part.
C-------------------
C
      IND = 0
      DO ISYMK = 1,NSYM
         DO K = 1,NRHF(ISYMK)
            ISYMDL = MULD2H(ISYMK,ISYDLK)
            DO ISYML = 1,NSYM
               ISYMD = MULD2H(ISYML,ISYMDL)
               DO L = 1,NRHF(ISYML)
                  DO D = 1,NVIR(ISYMD)
C
                     KOFFK = IRHF(ISYMK) + K
                     KOFFL = IRHF(ISYML) + L
                     KOFFD = IVIR(ISYMD) + D
C
                     OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFK)
C
                     IND = IND + 1
                     OCCHO(IND) = 
     &                   SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                  END DO
               END DO
            END DO
         END DO
      ENDDO
C
      NDIMOC = IND
C
C
C     Construct the vector.
C     ---------------------
C
      DO P = 1,NDIMOC
C
         OMEGA = OCCHO(P)
C
C
         OCCHO(P) = SQRT(TWO*CHOELE(ICHO))/(OMEGA+CHOELE(ICHO))
C
         DO JCHO = 1,ICHO-1
            OCCHO(P) = OCCHO(P) 
     &               *(OMEGA-CHOELE(JCHO))/(OMEGA+CHOELE(JCHO))
         END DO
C
      END DO
C
C------------------
C     Virtual part.
C------------------
C
C
      IND = 0
      DO ISYMB = 1,NSYM
         DO B = 1,NVIR(ISYMB)
            ISYMCJ = MULD2H(ISYMB,ISYCJB)
            DO ISYMJ = 1,NSYM
               ISYMC = MULD2H(ISYMC,ISYMCJ)
               DO J = 1,NRHF(ISYMJ)
                  DO C = 1,NVIR(ISYMC)
C
                     KOFFB = IVIR(ISYMK) + B
                     KOFFJ = IRHF(ISYMJ) + J
                     KOFFC = IVIR(ISYMC) + C
C
                     OME = FOCKD(KOFFC)-FOCKD(KOFFJ)+FOCKD(KOFFB)
C
                     IND = IND + 1
                     VICHO(IND) = 
     &                   SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                  END DO
               END DO
            END DO
         END DO
      ENDDO
C
      NDIMVI = IND
C
C
C     Construct the vector.
C     ---------------------
C
      DO P = 1,NDIMVI
C
         OMEGA = VICHO(P)
C
C
         VICHO(P) = SQRT(TWO*CHOELE(ICHO))/(OMEGA+CHOELE(ICHO))
C
         DO JCHO = 1,ICHO-1
            VICHO(P) = VICHO(P) 
     &               *(OMEGA-CHOELE(JCHO))/(OMEGA+CHOELE(JCHO))
         END DO
C
      END DO
C
      RETURN
      END
C
C  /* Deck ccho_decho3 */
      SUBROUTINE CCHO_DECHO3(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                       ISYDLK,ISYCJB,IOPTO,IOPTV)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Construct vector for Cholesky decomposition:
C
C     IF (IOPTO .EQ. 1):
C        For the occupied part, OCCHO(dlk)  update vector
C     IF (IOPTV .EQ. 1):
C        For the virtual  part, VICHO(cjb)  update vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C-------------------
C     Occupied part.
C-------------------
C
      IF (IOPTO .EQ. 1) THEN
C
         IF (ICHO .EQ. 1) THEN
            IND = 0
            DO ISYMK = 1,NSYM
               ISYMDL = MULD2H(ISYMK,ISYDLK)
               DO K = 1,NRHF(ISYMK)
                  DO ISYML = 1,NSYM
                     ISYMD = MULD2H(ISYML,ISYMDL)
                     DO L = 1,NRHF(ISYML)
                        DO D = 1,NVIR(ISYMD)
C
                           KOFFK = IRHF(ISYMK) + K
                           KOFFL = IRHF(ISYML) + L
                           KOFFD = IVIR(ISYMD) + D
C
                           OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFK)
C
                           IND = IND + 1
                           OCCHO(IND) = 
     &                         SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                        END DO
                     END DO
                  END DO
               END DO
            ENDDO
            NDIMOC = IND
         ELSE
            IND = 0
            DO ISYMK = 1,NSYM
               ISYMDL = MULD2H(ISYMK,ISYDLK)
               DO K = 1,NRHF(ISYMK)
                  DO ISYML = 1,NSYM
                     ISYMD = MULD2H(ISYML,ISYMDL)
                     DO L = 1,NRHF(ISYML)
                        DO D = 1,NVIR(ISYMD)
C
                           KOFFK = IRHF(ISYMK) + K
                           KOFFL = IRHF(ISYML) + L
                           KOFFD = IVIR(ISYMD) + D
C
                           OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFK)
C
                           IND = IND + 1
                           OCCHO(IND) =
     &                          (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                        END DO
                     END DO
                  END DO
               END DO
            ENDDO
C
            NDIMOC = IND
            FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
            CALL DSCAL(NDIMOC,FACTOR,OCCHO,1)
C
         ENDIF
C
      ENDIF
C
C------------------
C     Virtual part.
C------------------
C
      IF (IOPTV .EQ. 1) THEN
C
         IF (ICHO .EQ. 1) THEN
            IND = 0
            DO ISYMB = 1,NSYM
               ISYMCJ = MULD2H(ISYMB,ISYCJB)
               DO B = 1,NVIR(ISYMB)
                  DO ISYMJ = 1,NSYM
                     ISYMC = MULD2H(ISYMJ,ISYMCJ)
                     DO J = 1,NRHF(ISYMJ)
                        DO C = 1,NVIR(ISYMC)
C
                           KOFFB = IVIR(ISYMB) + B
                           KOFFJ = IRHF(ISYMJ) + J
                           KOFFC = IVIR(ISYMC) + C
C
                           OME = FOCKD(KOFFC)-FOCKD(KOFFJ)+FOCKD(KOFFB)
C
                           IND = IND + 1
                           VICHO(IND) = 
     &                         SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                        END DO
                     END DO
                  END DO
               END DO
            ENDDO
            NDIMVI = IND
C
         ELSE
C
            IND = 0
            DO ISYMB = 1,NSYM
               ISYMCJ = MULD2H(ISYMB,ISYCJB)
               DO B = 1,NVIR(ISYMB)
                  DO ISYMJ = 1,NSYM
                     ISYMC = MULD2H(ISYMJ,ISYMCJ)
                     DO J = 1,NRHF(ISYMJ)
                        DO C = 1,NVIR(ISYMC)
C
                           KOFFB = IVIR(ISYMB) + B
                           KOFFJ = IRHF(ISYMJ) + J
                           KOFFC = IVIR(ISYMC) + C
C
                           OME = FOCKD(KOFFC)-FOCKD(KOFFJ)+FOCKD(KOFFB)
C
                           IND = IND + 1
                           VICHO(IND) =
     &                          (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                        END DO
                     END DO
                  END DO
               END DO
            ENDDO
C
            NDIMVI = IND
            FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
            CALL DSCAL(NDIMVI,FACTOR,VICHO,1)
C
         ENDIF
C
      ENDIF

c     if (iopto .eq. 1) then
c        onorm = dsqrt(ddot(NDIMOC,occho,1,occho,1))
c        write(LUPRI,*) '  DECHO3: Dim. of occupied Chol. vec.: ',NDIMOC
c        write(LUPRI,*) '  DECHO3: Norm of occupied Chol. vec.: ',onorm
c     endif
c     if (ioptv .eq. 1) then
c        vnorm = dsqrt(ddot(NDIMVI,vicho,1,vicho,1))
c        write(LUPRI,*) '  DECHO3: Dim. of virtual  Chol. vec.: ',NDIMVI
c        write(LUPRI,*) '  DECHO3: Norm of virtual  Chol. vec.: ',vnorm
c     endif
C
      RETURN
      END
C  /* Deck ccho_jterm */
      SUBROUTINE CCHO_JTERM(XIINT,XJINT,XKINT,T2VO,FOCKD,NUMCHO,CHOELE,
     *                      WORK,LWORK,E4J,IA1,ISYMA,NUMIA,
     *                      FBATCH,PRINT,NONJ)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Calculate J-term (with identical contributions from the I-term
C     giving rise to the factors of 2 on the V and W intermediates):
C
C        E4J = E4J + Sum(a) Sum(ibj) [ 2 * V(b,ij;a) * Q(bj,i;a)
C                                    +     S(b,ji;a) * Q(bj,i;a)
C                                    + 2 * W(b,ij;a) * R(bj,i;a)
C                                    +     U(b,ji;a) * R(bj,i;a) ]
C
C     where
C
C        V(b,ij;a) =   Sum(dl) d(dlj) * (dj|il) * (la|bd)
C
C        S(b,ij;a) = - Sum(mn) d(bmn) * (bn|im) * (ma|nj)
C
C        W(b,ij;a) =   Sum(dl) d(dlj) * L(dlij) * (ld|ba)
C                    - Sum(dl) d(dlj) * (dl|ij) * (la|bd)
C
C        U(b,ji;a) =   Sum(mn) d(bmn) * (bn|im) * (na|mj)
C
C        Q(bj,i;a) =   Sum(ck) d(cka) * s(ck,ai) * s(ck,bj)
C
C        R(bj,i;a) =   Sum(ck) d(cka) * t(ck,ai) * s(ck,bj)
C                    - Sum(ck) d(cka) * t(ci,ak) * s(cj,bk)
C
C     and s(ai,bj) = 2 t(ai,bj) - t(aj,bi), and d(aij) and d(aib) denote
C     the occupied and virtual parts of the Cholesky decomposition of the
C     orbital energy denominator.
C
C     NOTE: The scaling of the integrals in the S and U intermediates
C           has potential operation count O^4V^2 (worst case: when
C           NUMIJ=1 for all symmetries).
C           Alternatively, one might calculate Q again and include
C           this term in the A term.....
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XIINT(*),XJINT(*),XKINT(*),T2VO(*),FOCKD(*)
      DIMENSION WORK(LWORK)
      LOGICAL   FBATCH,PRINT
C
      INTEGER IOFF2(8),IOFF3(8), IOFF4(8)
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C 
      TLAST = SECOND()
C 
C
      IF (PRINT) THEN
         TIMT = SECOND()
         WRITE(LUPRI,'(6X,A,/,6X,A,/)')
     &   'Calculation of the J term:',
     &   '=========================='
      ENDIF

c     lckai  = NCKI(ISYMA)*NUMIA
c     ldlba  = NCKATR(ISYMA)*NUMIA
c     ljild  = NTRAOC(1)
c     xinorm = dsqrt(ddot(ldlba,XIINT,1,XIINT,1))
c     xjnorm = dsqrt(ddot(ldlba,XJINT,1,XJINT,1))
c     xknorm = dsqrt(ddot(ljild,XKINT,1,XKINT,1))
c     xtnorm = dsqrt(ddot(NT2SQ(1),T2VO,1,T2VO,1))
c     write(LUPRI,*)
c     write(LUPRI,*) '   Entering JTERM:'
c     write(LUPRI,*) '   ==============='
c     write(LUPRI,*) '   IA1,NUMIA,ISYMA    : ',IA1,NUMIA,ISYMA
c     write(LUPRI,*) '   Norm of I(dl,b;#a) : ',xinorm
c     write(LUPRI,*) '   Norm of J(dl,b;#a) : ',xjnorm
c     write(LUPRI,*) '   Norm of K(ji,l;d)  : ',xknorm
c     write(LUPRI,*) '   Norm of T2VO       : ',xtnorm
c     write(LUPRI,*)

C
C--------------------------------------------------------------
C     Assign symmetries (T2VO and integrals assumed tot. sym.).
C--------------------------------------------------------------
C
      ISYCKI = ISYMA
      ISYDLB = ISYMA
C
C------------------------------
C     Set up index array IOFF2.
C------------------------------
C
      ICOUN2 = 0
      DO ISYM = 1,NSYM
         IOFF2(ISYM) = ICOUN2
         ICOUN2 = ICOUN2 + NT1AM(ISYM)
      ENDDO
      NTOT1 = ICOUN2

c     write(LUPRI,*) '   JTERM: NTOT1: ',NTOT1
C
C----------------
C     Allocation.
C----------------
C
      KCHOV = 1
      KT21  = KCHOV + NTOT1*NUMIA
      KT22  = KT21  + NCKI(ISYCKI)*NUMIA
      KT23  = KT22  + NCKI(ISYCKI)*NUMIA
      KT24  = KT23  + NCKI(ISYCKI)
      KEND0 = KT24  + NCKI(ISYCKI)
      LWRK0 = LWORK - KEND0 + 1
C
      IF (LWRK0 .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A)')
     &   'Insufficient memory in CCHO_JTERM'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND0-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in CCHO_JTERM')
      ENDIF
C
C----------------------------------
C     Find max. number of virtuals.
C----------------------------------
C
      MAXB = -1
      DO ISYM = 1,NSYM
         MAXB = MAX(MAXB,NVIR(ISYM))
      ENDDO

c     write(LUPRI,*) '   JTERM: MAXB: ',MAXB
C
C---------------------------------------------------------------------
C     Extract T21(ck,i;a) = T2VO(ca,ki) and T22(ck,i;a) = T2VO(ca,ik).
C---------------------------------------------------------------------
C
      CALL CCHO_JXT2(T2VO,WORK(KT21),WORK(KT22),IA1,NUMIA,ISYMA)

c     len    = NCKI(ISYMA)*NUMIA
c     t1norm = dnorm2(LEN,WORK(KT21),1)
c     t2norm = dnorm2(LEN,WORK(KT22),1)
c     write(LUPRI,*) '   JTERM: Norm of T21(ck,i;#a): ',t1norm
c     write(LUPRI,*) '   JTERM: Norm of T22(ck,i;#a): ',t2norm
C
C----------------------------
C     Loop over j-symmetries.
C----------------------------
C
      DO ISYMJ = 1,NSYM
C
         IF (NRHF(ISYMJ) .EQ. 0) GO TO 1000
C
         ISYDLI = ISYMJ
         ISYCKB = ISYMJ
         ISYMAJ = MULD2H(ISYMJ,ISYMA)
         ISYMBI = ISYMAJ
         ISYMMN = ISYMAJ
C
C-------------------
C        Allocation.
C-------------------
C
         KCHOL = KEND0
         KKIN4 = KCHOL + NMATIJ(ISYMMN)*NVIRT
         KKIN5 = KKIN4 + NMATIJ(ISYMMN)*MAXB
         KEND1 = KKIN5 + NMATIJ(ISYMMN)*MAXB
         LWRK1 = LWORK - KEND1 + 1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A)')
     &      'Insufficient memory in CCHO_JTERM'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in CCHO_JTERM')
         ENDIF
C
C---------------------------------
C        Find largest bi-subblock.
C---------------------------------
C
         MAXBI = -1
         DO ISYMI = 1,NSYM
            ISYMB = MULD2H(ISYMI,ISYMBI)
            MAXBI = MAX(MAXBI,NVIR(ISYMB)*NRHF(ISYMI))
         ENDDO

c     write(LUPRI,*) '   JTERM: MAXBI,ISYMBI: ',MAXBI,ISYMBI
C
C---------------------
C        Set up batch.
C---------------------
C
         MINMEM = 2*NCKATR(ISYCKB) + 3*NCKI(ISYDLI) + 4*MAXBI
     &          + 2*MAXB + NTOT1
         IF (FBATCH) THEN
            NEFJ  = MIN(NRHF(ISYMJ),NONJ)
            LEFF  = NEFJ*MINMEM + 1
            LWRK1 = MIN(LWRK1,LEFF)
         ENDIF
         NUMJ   = MIN(LWRK1/MINMEM,NRHF(ISYMJ))
C
         IF (NUMJ .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A)')
     &      'Insufficient memory in CCHO_JTERM'
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Calculated number of occupieds than can be treated is ',
     &      NUMJ
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Memory available for batch: ',LWRK1,
     &      'Minimum memory required   : ',MINMEM
            CALL QUIT('Insufficient memory in CCHO_JTERM')
         ENDIF
C
         NBATJ = (NRHF(ISYMJ) - 1)/NUMJ + 1
C
         IF (PRINT) THEN
            WRITE(LUPRI,'(9X,A,I1,A,/,9X,A)')
     &      'Batch over J, symmetry ',ISYMJ,':',
     &      '-------------------------'
            WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10)')
     &      'Minimum work space required   : ',MINMEM,
     &      'Work space available for batch: ',LWRK1
            WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10,/)')
     &      'Number of occupied orbitals   : ',NRHF(ISYMJ),
     &      'Required number of J-batches  : ',NBATJ
         ENDIF

c     write(LUPRI,*) '   JTERM: ISYMJ,NRHF(ISYMJ),NBATJ: ',
c    &           ISYMJ,NRHF(ISYMJ),NBATJ
C
C-------------------
C        Batch loop.
C-------------------
C
         DO IBATJ = 1,NBATJ
C
            NUMIJ = NUMJ
            IF (IBATJ .EQ. NBATJ) THEN
               NUMIJ = NRHF(ISYMJ) - NUMJ*(NBATJ - 1)
            ENDIF
C
            J1 = NUMJ*(IBATJ - 1) + 1
C
            IF (PRINT) THEN
               WRITE(LUPRI,'(12X,A,I10,A,/,12X,A)')
     &         'J-batch number ',IBATJ,':',
     &         '--------------------------'
               WRITE(LUPRI,'(12X,A,I10,1X,I10,/)')
     &         'First and last J: ',J1,J1+NUMIJ-1
            ENDIF
C
            KKIN1 = KEND1
            KKIN2 = KKIN1 + NCKI(ISYDLI)*NUMIJ
            KKIN3 = KKIN2 + NCKI(ISYDLI)*NUMIJ
            KT25  = KKIN3 + NCKI(ISYDLI)*NUMIJ
            KT26  = KT25  + NCKATR(ISYCKB)*NUMIJ
            KVMAT = KT26  + NCKATR(ISYCKB)*NUMIJ
            KWMAT = KVMAT + MAXBI*NUMIJ
            KQMAT = KWMAT + MAXBI*NUMIJ
            KRMAT = KQMAT + MAXBI*NUMIJ
            KSMAT = KRMAT + MAXBI*NUMIJ
            KUMAT = KSMAT + MAXB*NUMIJ
            KCHOO = KUMAT + MAXB*NUMIJ
            KEND2 = KCHOO + NTOT1*NUMIJ
            LWRK2 = LWORK - KEND2 + 1
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,'(//,5X,A)')
     &         'Bug detected in j-batching in CCHO_JTERM'
               CALL QUIT('Error in CCHO_JTERM')
            ENDIF
C
C--------------------------------------------------------------
C           Set up index array IOFF3: pointer to ck,b#j blocks.
C           Set up index array IOFF4: pointer to dl,i#j blocks.
C--------------------------------------------------------------
C
            ICOUN3 = 0
            ICOUN4 = 0
C
            DO ISYMB = 1,NSYM
C
               ISYMI  = ISYMB
               ISYMCK = MULD2H(ISYMB,ISYCKB)
               ISYMDL = MULD2H(ISYMI,ISYDLI)
C
               IOFF3(ISYMB) = ICOUN3
               IOFF4(ISYMI) = ICOUN4
C
               ICOUN3 = ICOUN3 + NT1AM(ISYMCK)*NVIR(ISYMB)*NUMIJ
               ICOUN4 = ICOUN4 + NT1AM(ISYMDL)*NRHF(ISYMI)*NUMIJ
C
            ENDDO

c     write(LUPRI,*) '   JTERM: IOFF2: ',(IOFF2(JSYM),JSYM=1,NSYM)
c     write(LUPRI,*) '   JTERM: IOFF3: ',(IOFF3(JSYM),JSYM=1,NSYM)
c     write(LUPRI,*) '   JTERM: IOFF4: ',(IOFF4(JSYM),JSYM=1,NSYM)
C
C---------------------------------------------------------
C           Extract occupied integrals according to IOFF4.
C---------------------------------------------------------
C
            CALL CCHO_JXK(XKINT,WORK(KKIN1),WORK(KKIN2),J1,NUMIJ,ISYMJ,
     &                    IOFF4)

c     len  = NCKI(ISYDLI)*NUMIJ
c     xk1nrm = dnorm2(len,WORK(KKIN1),1)
c     xk2nrm = dnorm2(len,WORK(KKIN2),1)
c     write(LUPRI,*) '   JTERM: Norm of K1(dl,i;#j): ',xk1nrm
c     write(LUPRI,*) '   JTERM: Norm of K2(dl,i;#j): ',xk2nrm
C
C----------------------------------------------------------
C           Set up 2CME of occ. integrals (IOFF4 ordering).
C----------------------------------------------------------
C
            LENTOT = NCKI(ISYDLI)*NUMIJ
            CALL DCOPY(LENTOT,WORK(KKIN2),1,WORK(KKIN3),1)
            CALL DSCAL(LENTOT,TWO,WORK(KKIN3),1)
            CALL DAXPY(LENTOT,XMONE,WORK(KKIN1),1,WORK(KKIN3),1)

c           len  = NCKI(ISYDLI)*NUMIJ
c           xk3nrm = dnorm2(len,WORK(KKIN1),1)
c           write(LUPRI,*) '   JTERM: Norm of K3(dl,i;#j): ',xk3nrm
C
C---------------------------------------------------------
C           Extract 2CME amplitudes, and their transposed.
C           Order according to IOFF3.
C---------------------------------------------------------
C
            CALL CCHO_JXT2CME(T2VO,WORK(KT25),WORK(KT26),J1,NUMIJ,ISYMJ,
     &                        IOFF3)

c           len    = NCKATR(ISYCKB)*NUMIJ
c           t5norm = dnorm2(len,WORK(KT25),1)
c           t6norm = dnorm2(len,WORK(KT26),1)
c           write(LUPRI,*) '   JTERM: Norm of T5(ck,b;#j): ',t5norm
c           write(LUPRI,*) '   JTERM: Norm of T6(ck,b;#j): ',t6norm
C
C-------------------------
C           Cholesky loop.
C-------------------------
C
            DO ICHO = 1,NUMCHO
C
               EJCOR = ZERO
C
C--------------------------------------------
C              Get Cholesky vectors.
C              (mn,b : actual, fixed ISYMMN)
C              (Occ. : update, vir. : actual)
C--------------------------------------------
C
               CALL CCHO_DECHO5(FOCKD,CHOELE,NUMCHO,ICHO,WORK(KCHOL),
     &                          WORK(KCHOO),WORK(KCHOV),ISYMMN,
     &                          IA1,NUMIA,ISYMA,J1,NUMIJ,ISYMJ)

c     call chk_vec1(FOCKD,CHOELE,NUMCHO,ICHO,
c    &              WORK(KCHOL),WORK(KCHOO),WORK(KCHOV),ISYMMN,
c    &              IA1,NUMIA,ISYMA,J1,NUMIJ,ISYMJ,WORK(KEND2),LWRK2)
C
C-------------------------------------------------------
C              Scale integrals with occ. update vectors.
C-------------------------------------------------------
C
               DO ISYMI = 1,NSYM
                  ISYMDL = MULD2H(ISYMI,ISYDLI)
                  DO J = 1,NUMIJ
                     KOFFO  = KCHOO + NTOT1*(J - 1) + IOFF2(ISYMDL)
                     KOFFK1 = KKIN1 + IOFF4(ISYMI)
     &                      + NT1AM(ISYMDL)*NRHF(ISYMI)*(J - 1)
                     KOFFK2 = KKIN2 + IOFF4(ISYMI)
     &                      + NT1AM(ISYMDL)*NRHF(ISYMI)*(J - 1)
                     KOFFK3 = KKIN3 + IOFF4(ISYMI)
     &                      + NT1AM(ISYMDL)*NRHF(ISYMI)*(J - 1)
                     CALL CCHO_SCVEC(WORK(KOFFK1),WORK(KOFFO),
     &                               NT1AM(ISYMDL),NRHF(ISYMI))
                     CALL CCHO_SCVEC(WORK(KOFFK2),WORK(KOFFO),
     &                               NT1AM(ISYMDL),NRHF(ISYMI))
                     CALL CCHO_SCVEC(WORK(KOFFK3),WORK(KOFFO),
     &                               NT1AM(ISYMDL),NRHF(ISYMI))
                  ENDDO
               ENDDO
C
C---------------------------------------------------------------
C              Calculate contributions in loop over available a.
C---------------------------------------------------------------
C
               DO A = 1,NUMIA
C
                  IA = IA1 + A - 1
C
C--------------------------------------------------
C                 Get original T(ck,i;a) subblocks.
C--------------------------------------------------
C
                  KOFT1 = KT21 + NCKI(ISYCKI)*(A - 1)
                  KOFT2 = KT22 + NCKI(ISYCKI)*(A - 1)
                  CALL DCOPY(NCKI(ISYCKI),WORK(KOFT1),1,WORK(KT23),1)
                  CALL DCOPY(NCKI(ISYCKI),WORK(KOFT2),1,WORK(KT24),1)
C
C-----------------------------------
C                 Loop i-symmetries.
C-----------------------------------
C
                  DO ISYMI = 1,NSYM
C
                     NI = NRHF(ISYMI)
                     IF (NI .LE. 0) GOTO 998
C
                     ISYMCK = MULD2H(ISYMI,ISYCKI)
                     ISYMB  = MULD2H(ISYMCK,ISYCKB)
                     ISYMDL = MULD2H(ISYMI,ISYDLI)
C
C--------------------------------------
C                    Calculate V and W.
C--------------------------------------
C
                     NDL   = NT1AM(ISYMDL)
                     NTODL = MAX(NDL,1)
                     NB    = NVIR(ISYMB)
                     NTOTB = MAX(NB,1)
                     NIJ   = NI*NUMIJ
C
                     KOFFI = NCKATR(ISYDLB)*(A - 1)
     &                     + ICKATR(ISYMDL,ISYMB) + 1
                     KOFK1 = KKIN1 + IOFF4(ISYMI)
C
                     CALL DGEMM('T','N',NB,NIJ,NDL,ONE,
     &                          XIINT(KOFFI),NTODL,WORK(KOFK1),NTODL,
     &                          ZERO,WORK(KVMAT),NTOTB)
C
                     KOFFJ = KOFFI
                     KOFK2 = KKIN2 + IOFF4(ISYMI)
                     KOFK3 = KKIN3 + IOFF4(ISYMI)
C
                     CALL DGEMM('T','N',NB,NIJ,NDL,XMONE,
     &                          XIINT(KOFFI),NTODL,WORK(KOFK2),NTODL,
     &                          ZERO,WORK(KWMAT),NTOTB)
                     CALL DGEMM('T','N',NB,NIJ,NDL,ONE,
     &                          XJINT(KOFFJ),NTODL,WORK(KOFK3),NTODL,
     &                          ONE,WORK(KWMAT),NTOTB)
C
C---------------------------------------------------------------
C                    Scale T23 and T24 blocks with full vectors.
C---------------------------------------------------------------
C
                     KOFFV = KCHOV + NTOT1*(A - 1) + IOFF2(ISYMCK)
                     KOFT3 = KT23  + ICKI(ISYMCK,ISYMI)
                     KOFT4 = KT24  + ICKI(ISYMCK,ISYMI)
                     CALL CCHO_SCVEC(WORK(KOFT3),WORK(KOFFV),
     &                               NT1AM(ISYMCK),NRHF(ISYMI))
                     CALL CCHO_SCVEC(WORK(KOFT4),WORK(KOFFV),
     &                               NT1AM(ISYMCK),NRHF(ISYMI))
C
C--------------------------------------
C                    Calculate Q and R.
C--------------------------------------
C
                     NCK   = NT1AM(ISYMCK)
                     NTOCK = MAX(NCK,1)
                     NBJ   = NVIR(ISYMB)*NUMIJ
                     NTOBJ = MAX(NBJ,1)
C
                     KOFT5 = KT25 + IOFF3(ISYMB)
C
                     CALL DGEMM('T','N',NBJ,NI,NCK,ONE,
     &                          WORK(KOFT5),NTOCK,WORK(KOFT3),NTOCK,
     &                          ZERO,WORK(KQMAT),NTOBJ)
C
                     CALL DCOPY(NBJ*NI,WORK(KQMAT),1,WORK(KRMAT),1)
C
                     CALL DGEMM('T','N',NBJ,NI,NCK,XMONE,
     &                          WORK(KOFT5),NTOCK,WORK(KOFT4),NTOCK,
     &                          TWO,WORK(KQMAT),NTOBJ)
C
                     KOFT6 = KT26 + IOFF3(ISYMB)
C
                     CALL DGEMM('T','N',NBJ,NI,NCK,ONE,
     &                          WORK(KOFT6),NTOCK,WORK(KOFT4),NTOCK,
     &                          ONE,WORK(KRMAT),NTOBJ)
C
                     DO I = 1,NI
C
C----------------------------------------------------------------------
C                       Get scaled occ. int. for S and U intermediates.
C----------------------------------------------------------------------
C
                        DO B = 1,NVIR(ISYMB)
C
                           IB = IVIR(ISYMB) + B - NRHFT
C
                           DO ISYMM = 1,NSYM
C
                              ISYMN = MULD2H(ISYMM,ISYMMN)
C
                              DO M = 1,NRHF(ISYMM)
                                 DO N = 1,NRHF(ISYMN)
C
                                    MN = IMATIJ(ISYMM,ISYMN)
     &                                 + NRHF(ISYMM)*(N - 1) + M
                                    NM = IMATIJ(ISYMN,ISYMM)
     &                                 + NRHF(ISYMN)*(M - 1) + N
C
                                    NMI = ISJIK(ISYMMN,ISYMI)
     &                                  + NMATIJ(ISYMMN)*(I - 1) + NM
C
                                    NMIB = ISJIKA(ISYMB,ISYMB)
     &                                   + NMAJIK(ISYMB)*(B - 1) + NMI
C
                                    KNMB = KCHOL
     &                                   + NMATIJ(ISYMMN)*(IB - 1)
     &                                   + NM - 1
                                    KOF4 = KKIN4
     &                                   + NMATIJ(ISYMMN)*(B - 1)
     &                                   + MN - 1
                                    KOF5 = KKIN5
     &                                   + NMATIJ(ISYMMN)*(B - 1)
     &                                   + NM - 1
C
                                    WORK(KOF4) = WORK(KNMB)*XKINT(NMIB)
                                    WORK(KOF5) = WORK(KOF4)
C
                                 ENDDO
C
                              ENDDO
C
                           ENDDO
C
                        ENDDO
C
C-----------------------------------------
C                       Calculate S and U.
C-----------------------------------------
C
                        NMN   = NMATIJ(ISYMMN)
                        NTOMN = MAX(NMN,1)
                        NB    = NVIR(ISYMB)
                        NTOTB = MAX(NB,1)
C
                        KOFFK = ISJIKA(ISYMA,ISYMA)
     &                        + NMAJIK(ISYMA)*(IA - 1)
     &                        + ISJIK(ISYMMN,ISYMJ)
     &                        + NMATIJ(ISYMMN)*(J1 - 1)
     &                        + 1
C
                        CALL DGEMM('T','N',NB,NUMIJ,NMN,XMONE,
     &                             WORK(KKIN4),NTOMN,XKINT(KOFFK),NTOMN,
     &                             ZERO,WORK(KSMAT),NTOTB)
C
                        CALL DGEMM('T','N',NB,NUMIJ,NMN,ONE,
     &                             WORK(KKIN5),NTOMN,XKINT(KOFFK),NTOMN,
     &                             ZERO,WORK(KUMAT),NTOTB)
C
                        DO J = 1,NUMIJ
C
C---------------------------------------------------------
C                          Calculate energy contributions.
C---------------------------------------------------------
C
                           KBJ   = NB*(J - 1)
                           KBJI  = NB*NUMIJ*(I - 1) + KBJ
                           KBIJ  = NB*NI*(J - 1) + NB*(I - 1)
C
                           KOFFQ = KQMAT + KBJI
                           KOFFV = KVMAT + KBIJ
                           KOFFS = KSMAT + KBJ
C
                           KOFFR = KRMAT + KBJI
                           KOFFW = KWMAT + KBIJ
                           KOFFU = KUMAT + KBJ
C
                           EJCOR = EJCOR
     &                        + TWO*DDOT(NB,WORK(KOFFQ),1,WORK(KOFFV),1)
     &                        + DDOT(NB,WORK(KOFFQ),1,WORK(KOFFS),1)
     &                        + TWO*DDOT(NB,WORK(KOFFR),1,WORK(KOFFW),1)
     &                        + DDOT(NB,WORK(KOFFR),1,WORK(KOFFU),1)
C
                        ENDDO
C
                     ENDDO
C
  998                CONTINUE
C
                  ENDDO
C
               ENDDO
C
               E4J = E4J + EJCOR
C 
             ENERGJ(ICHO) = ENERGJ(ICHO) + EJCOR
C 
C Decommented by Domenico
              TNOW = SECOND()
              DELTAT = TNOW - TLAST
              TLAST = TNOW
              SCNDSJ(ICHO) = SCNDSJ(ICHO) 
     &                      + DELTAT
C Decommented by Domenico
C
C
c              IF (PRINT) THEN
c                 WRITE(LUPRI,'(15X,A,I3,A,/,15X,A)')
c    &            'Status after Cholesky vector',ICHO,':',
c    &            '--------------------------------'
c                 IF (ABS(EJCOR) .LT. THRCHO) THEN
c                    WRITE(LUPRI,'(15X,A)') 'J term converged'
c                 ELSE
c                    WRITE(LUPRI,'(15X,A)') 'J term not converged'
c                 ENDIF
c                 TIM = SECOND() - TIMT
c                 WRITE(LUPRI,'(15X,A,F10.2,A,/)')
c    &            'Accumulated J-time: ',TIM,' seconds'
c              ENDIF
C
               IF (DABS(EJCOR) .LT. THRCHO) GOTO 999
C
            ENDDO
C
  999       CONTINUE
C
         ENDDO
C
 1000    CONTINUE
C
      ENDDO

c     lckai  = NCKI(ISYMA)*NUMIA
c     ldlba  = NCKATR(ISYMA)*NUMIA
c     ljild  = NTRAOC(1)
c     xinorm = dsqrt(ddot(ldlba,XIINT,1,XIINT,1))
c     xjnorm = dsqrt(ddot(ldlba,XJINT,1,XJINT,1))
c     xknorm = dsqrt(ddot(ljild,XKINT,1,XKINT,1))
c     xtnorm = dsqrt(ddot(NT2SQ(1),T2VO,1,T2VO,1))
c     write(LUPRI,*)
c     write(LUPRI,*) '   Exiting  JTERM:'
c     write(LUPRI,*) '   ==============='
c     write(LUPRI,*) '   IA1,NUMIA,ISYMA    : ',IA1,NUMIA,ISYMA
c     write(LUPRI,*) '   Norm of I(dl,b;#a) : ',xinorm
c     write(LUPRI,*) '   Norm of J(dl,b;#a) : ',xjnorm
c     write(LUPRI,*) '   Norm of K(ji,l;d)  : ',xknorm
c     write(LUPRI,*) '   Norm of T2VO       : ',xtnorm
c     write(LUPRI,*)

C
      RETURN
      END
C  /* Deck ccho_jxt2 */
      SUBROUTINE CCHO_JXT2(T2VO,T21,T22,IA1,NUMIA,ISYMA)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract subblocks of T2VO.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2VO(*), T21(*), T22(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYCKI = ISYMA
C
      DO ISYMKI = 1,NSYM
C
         ISYMCA = ISYMKI
         ISYMC  = MULD2H(ISYMKI,ISYCKI)
C
         IF (NVIR(ISYMC) .EQ. 0) GOTO 1000
C
         DO ISYMI = 1,NSYM
C
            ISYMK  = MULD2H(ISYMI,ISYMKI)
            ISYMCK = MULD2H(ISYMK,ISYMC)
C
            DO I = 1,NRHF(ISYMI)
               DO K = 1,NRHF(ISYMK)
C
                  KI = IMATIJ(ISYMK,ISYMI) + NRHF(ISYMK)*(I - 1) + K
                  IK = IMATIJ(ISYMI,ISYMK) + NRHF(ISYMI)*(K - 1) + I
C
                  DO A = 1,NUMIA
C
                     IA = IA1 + A - 1
C
                     KOFFS = NCKI(ISYCKI)*(A - 1) + ICKI(ISYMCK,ISYMI)
     &                     + NT1AM(ISYMCK)*(I - 1)
     &                     + IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1)
     &                     + 1
C
                     KOFFC = IT2VO(ISYMCA,ISYMKI)
     &                     + NMATAB(ISYMCA)*(KI - 1)
     &                     + IMATAB(ISYMC,ISYMA)
     &                     + NVIR(ISYMC)*(IA - 1) + 1
                     KOFFX = IT2VO(ISYMCA,ISYMKI)
     &                     + NMATAB(ISYMCA)*(IK - 1)
     &                     + IMATAB(ISYMC,ISYMA)
     &                     + NVIR(ISYMC)*(IA - 1) + 1
C
                     CALL DCOPY(NVIR(ISYMC),T2VO(KOFFC),1,T21(KOFFS),1)
                     CALL DCOPY(NVIR(ISYMC),T2VO(KOFFX),1,T22(KOFFS),1)
C
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDDO
C
 1000    CONTINUE
C
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     DO A = 1,NUMIA
c        DO ISYMI = 1,NSYM
c           ISYMCK = MULD2H(ISYMI,ISYMA)
c           DO I = 1,NRHF(ISYMI)
c              DO ISYMK = 1,NSYM
c                 ISYMC  = MULD2H(ISYMK,ISYMCK)
c                 ISYMCI = MULD2H(ISYMC,ISYMI)
c                 DO K = 1,NRHF(ISYMK)
c                    DO C = 1,NVIR(ISYMC)
c                       KCKIA = NCKI(ISYMA)*(A - 1)
c    &                        + ICKI(ISYMCK,ISYMI)
c    &                        + NT1AM(ISYMCK)*(I - 1)
c    &                        + IT1AM(ISYMC,ISYMK)
c    &                        + NVIR(ISYMC)*(K - 1)
c    &                        + C
c                       KCIKA = NCKI(ISYMA)*(A - 1)
c    &                        + ICKI(ISYMCI,ISYMK)
c    &                        + NT1AM(ISYMCI)*(K - 1)
c    &                        + IT1AM(ISYMC,ISYMI)
c    &                        + NVIR(ISYMC)*(I - 1)
c    &                        + C
c                       TEST  = T21(KCKIA) - T22(KCIKA)
c                       NTST  = NTST + 1
c                       IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*) '      JXT2: ISYMC,ISYMK,ISYMI,ISYMA: ',
c    &                        ISYMC,ISYMK,ISYMI,ISYMA
c     write(LUPRI,*) '      JXT2: C,K,I,A: ',C,K,I,A
c     write(LUPRI,*) '      JXT2: Diff.  : ',TEST
c     ICOUNE = ICOUNE + 1
c                       ENDIF
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     ndim = NCKI(ISYMA)*NUMIA
c     write(LUPRI,*) '      JXT2: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested (dim = ',ndim,')'
C
      RETURN
      END
C  /* Deck ccho_jxk */
      SUBROUTINE CCHO_JXK(XKINT,XKIN1,XKIN2,J1,NUMIJ,ISYMJ,IOFF4)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract occupied integrals.
C
#include <implicit.h>
#include <priunit.h>
      INTEGER DLIJ
      DIMENSION XKINT(*), XKIN1(*), XKIN2(*)
      INTEGER IOFF4(8)
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYDLI = ISYMJ
C
      DO ISYMI = 1,NSYM
C
         ISYMDL = MULD2H(ISYMI,ISYDLI)
         ISYMJI = MULD2H(ISYMI,ISYMJ)
C
         DO J = 1,NUMIJ
C
            JACT = J1 + J - 1
C
            DO I = 1,NRHF(ISYMI)
C
               JI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I - 1) + JACT
C
               DO ISYML = 1,NSYM
C
                  ISYMD  = MULD2H(ISYML,ISYMDL)
                  ISYJIL = MULD2H(ISYML,ISYMJI)
                  ISYMLI = MULD2H(ISYML,ISYMI)
C
                  DO L = 1,NRHF(ISYML)
C
                     LI  = IMATIJ(ISYML,ISYMI) + NRHF(ISYML)*(I - 1) + L
                     JIL = ISJIK(ISYMJI,ISYML)
     &                   + NMATIJ(ISYMJI)*(L - 1) + JI
                     LIJ = ISJIK(ISYMLI,ISYMJ)
     &                   + NMATIJ(ISYMLI)*(JACT - 1) + LI
C
                     DO D = 1,NVIR(ISYMD)
C
                        DLIJ = IOFF4(ISYMI)
     &                       + NT1AM(ISYMDL)*NRHF(ISYMI)*(J - 1)
     &                       + NT1AM(ISYMDL)*(I - 1)
     &                       + IT1AM(ISYMD,ISYML)
     &                       + NVIR(ISYMD)*(L - 1) + D
         
C
                        KOFD = ISJIKA(ISYJIL,ISYMD)
     &                       + NMAJIK(ISYJIL)*(D - 1)
C
                        JILD = KOFD + JIL
                        LIJD = KOFD + LIJ
C
                        XKIN1(DLIJ) = XKINT(JILD)
                        XKIN2(DLIJ) = XKINT(LIJD)
C
                     ENDDO
C
                  ENDDO
C
               ENDDO
C
            ENDDO
C
         ENDDO
C
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     ISYML  = ISYMJ
c     DO ISYMI = 1,NSYM
c        ISYMDL = MULD2H(ISYMI,ISYDLI)
c        ISYMD  = MULD2H(ISYMDL,ISYML)
c        ISYMDJ = MULD2H(ISYMD,ISYMJ)
c        DO J = 1,NUMIJ
c           JJ = J1 + J - 1
c           DO I = 1,NRHF(ISYMI)
c              DO L = 1,NUMIJ
c                 LL = J1 + L - 1
c                 DO D = 1,NVIR(ISYMD)
c                    KDLIJ = IOFF4(ISYMI)
c    &                     + NT1AM(ISYMDL)*NRHF(ISYMI)*(J - 1)
c    &                     + NT1AM(ISYMDL)*(I - 1)
c    &                     + IT1AM(ISYMD,ISYML)
c    &                     + NVIR(ISYMD)*(LL - 1) + D
c                    KDJIL = IOFF4(ISYMI)
c    &                     + NT1AM(ISYMDJ)*NRHF(ISYMI)*(L - 1)
c    &                     + NT1AM(ISYMDJ)*(I - 1)
c    &                     + IT1AM(ISYMD,ISYMJ)
c    &                     + NVIR(ISYMD)*(JJ - 1) + D
c                    TEST  = XKIN1(KDLIJ) - XKIN2(KDJIL)
c                    IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*) '      JXK: ISYMD,ISYML,ISYMI,ISYMJ: ',
c    &                       ISYMD,ISYML,ISYMI,ISYMJ
c     write(LUPRI,*) '      JXK: D,L,I,J: ',D,L,I,J
c     write(LUPRI,*) '      JXK: Diff.  : ',TEST
c     ICOUNE = ICOUNE + 1
c                    ENDIF
c                    NTST = NTST + 1
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     ndim = NCKI(ISYDLI)*NUMIJ
c     write(LUPRI,*) '      JXK: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested (dim = ',ndim,')'
C
      RETURN
      END
C  /* Deck ccho_jxt2cme */
      SUBROUTINE CCHO_JXT2CME(T2VO,T25,T26,J1,NUMIJ,ISYMJ,IOFF3)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract T2 amplitudes as 2CME and their transposed.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2VO(*), T25(*), T26(*)
      INTEGER IOFF3(8)
      INTEGER BJ
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (XMONE = -1.00D0, TWO = 2.00D0)
C
      ISYCKB = ISYMJ
C
      LENTOT = NCKATR(ISYCKB)*NUMIJ
      CALL DZERO(T25,LENTOT)
      CALL DZERO(T26,LENTOT)
C
      DO ISYMKJ = 1,NSYM
C
         ISYMCB = ISYMKJ
         ISYMK  = MULD2H(ISYMKJ,ISYMJ)
C
         DO J = 1,NUMIJ
C
            JACT = J1 + J - 1
C
            DO K = 1,NRHF(ISYMK)
C
               KJ = IMATIJ(ISYMK,ISYMJ) + NRHF(ISYMK)*(JACT - 1) + K
               JK = IMATIJ(ISYMJ,ISYMK) + NRHF(ISYMJ)*(K - 1) + JACT
C
               DO ISYMB = 1,NSYM
C
                  ISYMC  = MULD2H(ISYMB,ISYMCB)
                  ISYMCK = MULD2H(ISYMC,ISYMK)
C
                  IF (NVIR(ISYMC) .GT. 0) THEN
C
                     KOFCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1)
     &                     + 1
C
                     DO B = 1,NVIR(ISYMB)
C
                        KOFCB = IMATAB(ISYMC,ISYMB)
     &                        + NVIR(ISYMC)*(B - 1) + 1
                        KOFF1 = IT2VO(ISYMCB,ISYMKJ)
     &                        + NMATAB(ISYMCB)*(KJ - 1) + KOFCB
                        KOFF2 = IT2VO(ISYMCB,ISYMKJ)
     &                        + NMATAB(ISYMCB)*(JK - 1) + KOFCB
C
                        BJ    = NVIR(ISYMB)*(J - 1) + B
                        KOFFS = IOFF3(ISYMB) + NT1AM(ISYMCK)*(BJ - 1)
     &                        + KOFCK
C
                        CALL DAXPY(NVIR(ISYMC),TWO,T2VO(KOFF1),1,
     &                                             T25(KOFFS),1)
                        CALL DAXPY(NVIR(ISYMC),XMONE,T2VO(KOFF2),1,
     &                                               T25(KOFFS),1)
C
                        CALL DAXPY(NVIR(ISYMC),TWO,T2VO(KOFF2),1,
     &                                             T26(KOFFS),1)
                        CALL DAXPY(NVIR(ISYMC),XMONE,T2VO(KOFF1),1,
     &                                               T26(KOFFS),1)
C
                     ENDDO
C
                  ENDIF
C
               ENDDO
C
            ENDDO
C
         ENDDO
C
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     DO ISYMB = 1,NSYM
c        ISYMCK = MULD2H(ISYMB,ISYMJ)
c        DO J = 1,NUMIJ
c           DO B = 1,NVIR(ISYMB)
c              DO ISYMK = 1,NSYM
c                 ISYMC  = MULD2H(ISYMK,ISYMCK)
c                 ISYMBK = MULD2H(ISYMK,ISYMB)
c                 DO K = 1,NRHF(ISYMK)
c                    DO C = 1,NVIR(ISYMC)
c                       KCKBJ = IOFF3(ISYMB)
c    &                        + NT1AM(ISYMCK)*NVIR(ISYMB)*(J - 1)
c    &                        + NT1AM(ISYMCK)*(B - 1)
c    &                        + IT1AM(ISYMC,ISYMK)
c    &                        + NVIR(ISYMC)*(K - 1) + C
c                       KBKCJ = IOFF3(ISYMC)
c    &                        + NT1AM(ISYMBK)*NVIR(ISYMC)*(J - 1)
c    &                        + NT1AM(ISYMBK)*(C - 1)
c    &                        + IT1AM(ISYMB,ISYMK)
c    &                        + NVIR(ISYMB)*(K - 1) + B
c                       NTST  = NTST + 1
c                       TEST  = T25(KCKBJ) - T26(KBKCJ)
c                       IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*) '      JXT2CME: ISYMC,ISYMK,ISYMB,ISYMJ: ',
c    &                           ISYMC,ISYMK,ISYMB,ISYMJ
c     write(LUPRI,*) '      JXT2CME: C,K,B,J: ',C,K,B,J
c     write(LUPRI,*) '      JXT2CME: Diff.  : ',TEST
c     ICOUNE = ICOUNE + 1
c                       ENDIF
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     ndim = NCKATR(ISYMJ)*NUMIJ
c     write(LUPRI,*) '      JXT2CME: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested (dim = ',ndim,')'
C
      RETURN
      END
C  /* Deck ccho_jrsrt */
      SUBROUTINE CCHO_JRSRT(XJINT,XJ1,XJ2,NUMIB,ISYMB,IOFF2)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Resort virtual integrals and set up corresponding index array
C     IOFF2.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XJINT(*),XJ1(*),XJ2(*)
      INTEGER IOFF2(8)
      INTEGER DE,EI,DI,BI,DEBI,EIDB,DIEB
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYDEI = ISYMB
C
C------------------------
C     Set up index array.
C------------------------
C
      ICOUNT = 0
C
      DO ISYMI = 1,NSYM
C
         ISYMDE = MULD2H(ISYMI,ISYDEI)
C
         IOFF2(ISYMI) = ICOUNT
C
         ICOUNT = ICOUNT + NMATAB(ISYMDE)*NUMIB*NRHF(ISYMI)
C
      ENDDO

c     ndim = NCKATR(ISYDEI)*NUMIB
c     if (ICOUNT .NE. ndim) then
c        write(LUPRI,*) '   JRSRT: ICOUNT,ndim: ',ICOUNT,ndim
c        CALL QUIT(' dimension problem in JRSRT ')
c     endif
c     write(LUPRI,*) '   JRSRT: IOFF2: ',(IOFF2(ISYM),ISYM=1,NSYM)
C
C------------
C     Resort.
C------------
C
      DO ISYMI = 1,NSYM
C
         ISYMDE = MULD2H(ISYMI,ISYDEI)
C
         DO I = 1,NRHF(ISYMI)
            DO B = 1,NUMIB
C
               BI = NUMIB*(I - 1) + B
C
               DO ISYME = 1,NSYM
C
                  ISYMD  = MULD2H(ISYME,ISYMDE)
                  ISYMDI = MULD2H(ISYMD,ISYMI)
                  ISYMEI = MULD2H(ISYME,ISYMI)
C
                  DO E = 1,NVIR(ISYME)
C
                     EI = IT1AM(ISYME,ISYMI) + NVIR(ISYME)*(I - 1) + E
C
                     DO D = 1,NVIR(ISYMD)
C
                        DE = IMATAB(ISYMD,ISYME) + NVIR(ISYMD)*(E - 1)
     &                     + D
                        DI = IT1AM(ISYMD,ISYMI)  + NVIR(ISYMD)*(I - 1)
     &                     + D
C
                        DIEB = NCKATR(ISYDEI)*(B - 1)
     &                       + ICKATR(ISYMDI,ISYME)
     &                       + NT1AM(ISYMDI)*(E - 1) + DI
                        EIDB = NCKATR(ISYDEI)*(B - 1)
     &                       + ICKATR(ISYMEI,ISYMD)
     &                       + NT1AM(ISYMEI)*(D - 1) + EI
C
                        DEBI = IOFF2(ISYMI) + NMATAB(ISYMDE)*(BI - 1)
     &                       + DE
C
                        XJ1(DEBI) = XJINT(EIDB)
                        XJ2(DEBI) = XJINT(DIEB)
C
                     ENDDO
C
C
                  ENDDO
C
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     DO ISYMI = 1,NSYM
c        ISYMDE = MULD2H(ISYMI,ISYDEI)
c        DO I = 1,NRHF(ISYMI)
c           DO B = 1,NUMIB
c              NBI = NUMIB*(I - 1) + B
c              DO ISYME = 1,NSYM
c                 ISYMD = MULD2H(ISYME,ISYMDE)
c                 DO E = 1,NVIR(ISYME)
c                    DO D = 1,NVIR(ISYMD)
c                       NDE = IMATAB(ISYMD,ISYME)
c    &                      + NVIR(ISYMD)*(E - 1) + D
c                       NED = IMATAB(ISYME,ISYMD)
c    &                      + NVIR(ISYME)*(D - 1) + E
c                       KDEBI = IOFF2(ISYMI)
c    &                        + NMATAB(ISYMDE)*(NBI - 1)
c    &                        + NDE
c                       KEDBI = IOFF2(ISYMI)
c    &                        + NMATAB(ISYMDE)*(NBI - 1)
c    &                        + NED
c                       TEST  = XJ1(KDEBI) - XJ2(KEDBI)
c                       NTST  = NTST + 1
c                       IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*)
c    & '   JRSRT: ISYMD,ISYME,ISYMB,ISYMI: ',
c    & ISYMD,ISYME,ISYMB,ISYMI
c     write(LUPRI,*)
c    & '   JRSRT: D,E,B,I: ',D,E,IB1+B-1,I
c     write(LUPRI,*)
c    & '   JRSRT: Diff.  : ',TEST
c     ICOUNE = ICOUNE + 1
c                       ENDIF
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     write(LUPRI,*) '   JRSRT: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested'
c     write(LUPRI,*) '   JRSRT: (ICOUNT = ',ICOUNT,')'

c     ndim   = NCKATR(ISYDEI)*NUMIB
c     xjnorm = dsqrt(ddot(ndim,XJINT,1,XJINT,1))
c     xj1nrm = dsqrt(ddot(ndim,XJ1,1,XJ1,1))
c     xj2nrm = dsqrt(ddot(ndim,XJ2,1,XJ2,1))
c     write(LUPRI,*) '   JRSRT: Norm of  J(ei,d;#b): ',xjnorm
c     write(LUPRI,*) '   JRSRT: Norm of J1(de,#b,i): ',xj1nrm
c     write(LUPRI,*) '   JRSRT: Norm of J2(de,#b,i): ',xj2nrm
C
      RETURN
      END
C  /* Deck ccho_ix2cme */
      SUBROUTINE CCHO_IX2CME(T2VO,T2AM3,T2AM4,IB1,NUMIB,ISYMB,IOFF3)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Sort amplitudes and set up corresponding index array
C     IOFF3.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2VO(*),T2AM3(*),T2AM4(*)
      INTEGER IOFF3(8)
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (XMONE = -1.00D0, TWO = 2.00D0)
C
      ISYCKJ = ISYMB
C
C----------------
C     Initialize.
C----------------
C
      LENTOT = NCKI(ISYCKJ)*NUMIB
      CALL DZERO(T2AM3,LENTOT)
      CALL DZERO(T2AM4,LENTOT)
C
C------------------
C     Set up IOFF3.
C------------------
C
      ICOUNT = 0
C
      DO ISYMJ = 1,NSYM
C
         ISYMCK = MULD2H(ISYMJ,ISYCKJ)
C
         IOFF3(ISYMJ) = ICOUNT
C
         ICOUNT = ICOUNT + NT1AM(ISYMCK)*NUMIB*NRHF(ISYMJ)
C
      ENDDO

c     if (ICOUNT .NE. LENTOT) then
c        write(LUPRI,*) '   IX2CME: Dimension problem, ICOUNT,LENTOT: ',
c    &              ICOUNT,LENTOT
c        CALL QUIT(' dimension problem in IX2CME ')
c     endif
c     write(LUPRI,*) '   IX2CME: IOFF3: ',(IOFF3(ISYM),ISYM=1,NSYM)
C
C-------------
C     Extract.
C-------------
C
      DO ISYMJ = 1,NSYM
C
         ISYMCK = MULD2H(ISYMJ,ISYCKJ)
C
         DO J = 1,NRHF(ISYMJ)
            DO IB = 1,NUMIB
C
               B = IB1 + IB - 1
C
               IBJ = NUMIB*(J - 1) + IB
C
               DO ISYMK = 1,NSYM
C
                  ISYMC  = MULD2H(ISYMK,ISYMCK)
                  ISYMCB = MULD2H(ISYMC,ISYMB)
                  ISYMKJ = MULD2H(ISYMK,ISYMJ)
C
                  IF (NVIR(ISYMC) .GT. 0) THEN
C
                     DO K = 1,NRHF(ISYMK)
C
                        KJ = IMATIJ(ISYMK,ISYMJ) + NRHF(ISYMK)*(J - 1)
     &                     + K
                        JK = IMATIJ(ISYMJ,ISYMK) + NRHF(ISYMJ)*(K - 1)
     &                     + J
C
                        KOFFC = IT2VO(ISYMCB,ISYMKJ)
     &                        + NMATAB(ISYMCB)*(KJ - 1)
     &                        + IMATAB(ISYMC,ISYMB)
     &                        + NVIR(ISYMC)*(B - 1) + 1
                        KOFFX = IT2VO(ISYMCB,ISYMKJ)
     &                        + NMATAB(ISYMCB)*(JK - 1)
     &                        + IMATAB(ISYMC,ISYMB)
     &                        + NVIR(ISYMC)*(B - 1) + 1
C
                        KOFFT = IOFF3(ISYMJ) + NT1AM(ISYMCK)*(IBJ - 1)
     &                        + IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1)
     &                        + 1
C
                        CALL DAXPY(NVIR(ISYMC),TWO,T2VO(KOFFC),1,
     &                                             T2AM3(KOFFT),1)
                        CALL DAXPY(NVIR(ISYMC),XMONE,T2VO(KOFFX),1,
     &                                               T2AM3(KOFFT),1)
C
                        CALL DAXPY(NVIR(ISYMC),TWO,T2VO(KOFFX),1,
     &                                             T2AM4(KOFFT),1)
                        CALL DAXPY(NVIR(ISYMC),XMONE,T2VO(KOFFC),1,
     &                                               T2AM4(KOFFT),1)
C
                     ENDDO
C
                  ENDIF
C
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     DO ISYMJ = 1,NSYM
c        ISYMCK = MULD2H(ISYMJ,ISYCKJ)
c        DO J = 1,NRHF(ISYMJ)
c           DO B = 1,NUMIB
c              NBJ = NUMIB*(J - 1) + B
c              DO ISYMK = 1,NSYM
c                 ISYMC  = MULD2H(ISYMK,ISYMCK)
c                 ISYMCJ = MULD2H(ISYMK,ISYCKJ)
c                 DO K = 1,NRHF(ISYMK)
c                    NBK = NUMIB*(K - 1) + B
c                    DO C = 1,NVIR(ISYMC)
c                       KCKBJ = IOFF3(ISYMJ)
c    &                        + NT1AM(ISYMCK)*(NBJ - 1)
c    &                        + IT1AM(ISYMC,ISYMK)
c    &                        + NVIR(ISYMC)*(K - 1) + C
c                       KCJBK = IOFF3(ISYMK)
c    &                        + NT1AM(ISYMCJ)*(NBK - 1)
c    &                        + IT1AM(ISYMC,ISYMJ)
c    &                        + NVIR(ISYMC)*(J - 1) + C
c                       TEST  = T2AM3(KCKBJ) - T2AM4(KCJBK)
c                       NTST  = NTST + 1
c                       IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*)
c    & '   IX2CME: ISYMC,ISYMK,ISYMB,ISYMJ: ',
c    & ISYMC,ISYMK,ISYMB,ISYMJ
c     write(LUPRI,*)
c    & '   IX2CME: C,K,B,J: ',C,K,IB1+B-1,J
c     write(LUPRI,*)
c    & '   IX2CME: Diff.  : ',TEST
c     ICOUNE = ICOUNE + 1
c                       ENDIF
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     write(LUPRI,*) '   IX2CME: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested'
c     write(LUPRI,*) '   IX2CME: (ICOUNT = ',ICOUNT,')'
c     len    = NCKI(ISYMB)*NUMIB
c     t3norm = dsqrt(ddot(len,T2AM3,1,T2AM3,1))
c     t4norm = dsqrt(ddot(len,T2AM4,1,T2AM4,1))
c     write(LUPRI,*) '   IX2CME: Norm of T3: ',t3norm
c     write(LUPRI,*) '   IX2CME: Norm of T4: ',t4norm
C
      RETURN
      END
C  /* Deck ccho_ixktil */
      SUBROUTINE CCHO_IXKTIL(XKINT,XKTIL,IOFSET,IOFF4,NUMIA,ISYMA)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Sort virtual integrals for I term.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XKINT(*), XKTIL(*)
      INTEGER IOFSET(8),IOFF4(8)
      INTEGER AJ, EA
#include <ccorb.h>
#include <ccsdsym.h>
C
      DO ISYMJ = 1,NSYM
C
         ISYMDE = MULD2H(ISYMJ,ISYMA)
C
         DO J = 1,NRHF(ISYMJ)
            DO A = 1,NUMIA
C
               AJ = NUMIA*(J - 1) + A
C
               DO ISYME = 1,NSYM
C
                  ISYMD  = MULD2H(ISYME,ISYMDE)
                  ISYMDJ = MULD2H(ISYMD,ISYMJ)
C
                  IF (NVIR(ISYMD) .GT. 0) THEN
C
                     DO E = 1,NVIR(ISYME)
C
                        EA = NVIR(ISYME)*(A - 1) + E
C
                        KOFF1 = IOFSET(ISYME)
     &                        + NT1AM(ISYMDJ)*(EA - 1)
     &                        + IT1AM(ISYMD,ISYMJ)
     &                        + NVIR(ISYMD)*(J - 1) + 1
                        KOFF2 = IOFF4(ISYMJ)
     &                        + NMATAB(ISYMDE)*(AJ - 1)
     &                        + IMATAB(ISYMD,ISYME)
     &                        + NVIR(ISYMD)*(E - 1) + 1
C
                        CALL DCOPY(NVIR(ISYMD),XKINT(KOFF1),1,
     &                                         XKTIL(KOFF2),1)
C
                     ENDDO
C
                  ENDIF
C
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO

c     ndim   = NCKATR(ISYMA)*NUMIA
c     xknorm = dsqrt(ddot(ndim,XKINT,1,XKINT,1))
c     xktnrm = dsqrt(ddot(ndim,XKTIL,1,XKTIL,1))
c     test   = xknorm - xktnrm
c      if (dabs(test) .gt. 1.0d-15) then
c        write(LUPRI,*) '   IXKTIL: Norm of  nK(dj,e,#a): ',xknorm
c        write(LUPRI,*) '   IXKTIL: Norm of nKT(de,#a,j): ',xktnrm
c      endif
C
      RETURN
      END
C  /* Deck ccho_ixt2 */
      SUBROUTINE CCHO_IXT2(T2VO,T2AM1,T2AM2,IA1,NUMIA,ISYMA,IOFF5)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Sort amplitudes and set up corresponding index array
C     IOFF5.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2VO(*),T2AM1(*),T2AM2(*)
      INTEGER IOFF5(8)
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYCKI = ISYMA
C
C------------------
C     Set up IOFF5.
C------------------
C
      ICOUNT = 0
C
      DO ISYMI = 1,NSYM
C
         ISYMCK = MULD2H(ISYMI,ISYCKI)
C
         IOFF5(ISYMI) = ICOUNT
C
         ICOUNT = ICOUNT + NT1AM(ISYMCK)*NUMIA*NRHF(ISYMI)
C
      ENDDO

c     ndim = NCKI(ISYCKI)*NUMIA
c     if (ICOUNT .NE. ndim) then
c        write(LUPRI,*) '   IXT2: dim. error: ICOUNT,ndim: ',ICOUNT,ndim
c        CALL QUIT(' dimension problem in IXT2 ')
c     endif
c     write(LUPRI,*) '   IXT2: IOFF5: ',(IOFF5(ISYM),ISYM=1,NSYM)
C
C-------------
C     Extract.
C-------------
C
      DO ISYMI = 1,NSYM
C
         ISYMCK = MULD2H(ISYMI,ISYCKI)
C
         DO I = 1,NRHF(ISYMI)
            DO IA = 1,NUMIA
C
               A = IA1 + IA - 1
C
               IAI = NUMIA*(I - 1) + IA
C
               DO ISYMK = 1,NSYM
C
                  ISYMC  = MULD2H(ISYMK,ISYMCK)
                  ISYMCA = MULD2H(ISYMC,ISYMA)
                  ISYMKI = MULD2H(ISYMK,ISYMI)
C
                  IF (NVIR(ISYMC) .GT. 0) THEN
C
                     DO K = 1,NRHF(ISYMK)
C
                        KI = IMATIJ(ISYMK,ISYMI) + NRHF(ISYMK)*(I - 1)
     &                     + K
                        IK = IMATIJ(ISYMI,ISYMK) + NRHF(ISYMI)*(K - 1)
     &                     + I
C
                        KOFFC = IT2VO(ISYMCA,ISYMKI)
     &                        + NMATAB(ISYMCA)*(KI - 1)
     &                        + IMATAB(ISYMC,ISYMA)
     &                        + NVIR(ISYMC)*(A - 1) + 1
                        KOFFX = IT2VO(ISYMCA,ISYMKI)
     &                        + NMATAB(ISYMCA)*(IK - 1)
     &                        + IMATAB(ISYMC,ISYMA)
     &                        + NVIR(ISYMC)*(A - 1) + 1
C
                        KOFFT = IOFF5(ISYMI) + NT1AM(ISYMCK)*(IAI - 1)
     &                        + IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1)
     &                        + 1
C
                        CALL DCOPY(NVIR(ISYMC),T2VO(KOFFC),1,
     &                                         T2AM1(KOFFT),1)
                        CALL DCOPY(NVIR(ISYMC),T2VO(KOFFX),1,
     &                                         T2AM2(KOFFT),1)
C
                     ENDDO
C
                  ENDIF
C
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     DO ISYMI = 1,NSYM
c        ISYMCK = MULD2H(ISYMI,ISYCKI)
c        DO I = 1,NRHF(ISYMI)
c           DO A = 1,NUMIA
c              NAI = NUMIA*(I - 1) + A
c              DO ISYMK = 1,NSYM
c                 ISYMC  = MULD2H(ISYMK,ISYMCK)
c                 ISYMCI = MULD2H(ISYMK,ISYCKI)
c                 DO K = 1,NRHF(ISYMK)
c                    NAK = NUMIA*(K - 1) + A
c                    DO C = 1,NVIR(ISYMC)
c                       KCKAI = IOFF5(ISYMI)
c    &                        + NT1AM(ISYMCK)*(NAI - 1)
c    &                        + IT1AM(ISYMC,ISYMK)
c    &                        + NVIR(ISYMC)*(K - 1) + C
c                       KCIAK = IOFF5(ISYMK)
c    &                        + NT1AM(ISYMCI)*(NAK - 1)
c    &                        + IT1AM(ISYMC,ISYMI)
c    &                        + NVIR(ISYMC)*(I - 1) + C
c                       TEST  = T2AM1(KCKAI) - T2AM2(KCIAK)
c                       NTST  = NTST + 1
c                       IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*)
c    & '   IXT2: ISYMC,ISYMK,ISYMA,ISYMI: ',
c    & ISYMC,ISYMK,ISYMA,ISYMI
c     write(LUPRI,*)
c    & '   IXT2: C,K,A,I: ',C,K,IA1+A-1,I
c     write(LUPRI,*)
c    & '   IXT2: Diff.  : ',TEST
c     ICOUNE = ICOUNE + 1
c                       ENDIF
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     write(LUPRI,*) '   IXT2: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested'
c     write(LUPRI,*) '   IXT2: (ICOUNT = ',ICOUNT,')'
c     len    = NCKI(ISYMA)*NUMIA
c     t1norm = dsqrt(ddot(len,T2AM1,1,T2AM1,1))
c     t2norm = dsqrt(ddot(len,T2AM2,1,T2AM2,1))
c     write(LUPRI,*) '   IXT2: Norm of T1: ',t1norm
c     write(LUPRI,*) '   IXT2: Norm of T2: ',t2norm
C
      RETURN
      END
C  /* Deck ccho_isclt */
      SUBROUTINE CCHO_ISCLT(T2AM1,T2AM2,IOFF5,CHOVEC,NUMIA,ISYMA)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Scale amplitudes for I term.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2AM1(*),T2AM2(*),CHOVEC(*)
      INTEGER IOFF5(8)
      INTEGER AI,CK,CKAI
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYCKI = ISYMA

c     ndim   = NCKI(ISYCKI)*NUMIA
c     t1norm = dsqrt(ddot(ndim,T2AM1,1,T2AM1,1))
c     t2norm = dsqrt(ddot(ndim,T2AM2,1,T2AM2,1))
c     diff   = t2norm - t1norm
c      if (dabs(diff) .gt. 1.0d-15) then
c        write(LUPRI,*) '   ISCLT: Norm of 1T(ck,#ai)  [input]: ',t1norm
c        write(LUPRI,*) '   ISCLT: Norm of 2T(ck,#ai)  [input]: ',t2norm
c      endif
C
      DO ISYMI = 1,NSYM
C
         ISYMCK = MULD2H(ISYMI,ISYCKI)
C
         DO I = 1,NRHF(ISYMI)
            DO A = 1,NUMIA
C
               AI = NUMIA*(I - 1) + A
C
               DO ISYMK = 1,NSYM
C
                  ISYMC  = MULD2H(ISYMK,ISYMCK)
                  ISYMCI = MULD2H(ISYMK,ISYCKI)
C
                  IF (NVIR(ISYMC) .GT. 0) THEN
                     DO K = 1,NRHF(ISYMK)
C
                        KOFFT = IOFF5(ISYMI) + NT1AM(ISYMCK)*(AI - 1)
     &                        + IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1)
     &                        + 1
C
                        KOFFV = ICKI(ISYMCK,ISYMI)
     &                        + NT1AM(ISYMCK)*(I - 1)
     &                        + IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1)
     &                        + 1
                        CALL CCHO_SCVEC(T2AM1(KOFFT),CHOVEC(KOFFV),
     &                                  NVIR(ISYMC),1)
C
                        KOFFV = ICKI(ISYMCI,ISYMK)
     &                        + NT1AM(ISYMCI)*(K - 1)
     &                        + IT1AM(ISYMC,ISYMI) + NVIR(ISYMC)*(I - 1)
     &                        + 1
                        CALL CCHO_SCVEC(T2AM2(KOFFT),CHOVEC(KOFFV),
     &                                  NVIR(ISYMC),1)
C
                     ENDDO
                  ENDIF
C
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     DO ISYMI = 1,NSYM
c        ISYMCK = MULD2H(ISYMI,ISYCKI)
c        DO I = 1,NRHF(ISYMI)
c           DO A = 1,NUMIA
c              NAI = NUMIA*(I - 1) + A
c              DO ISYMK = 1,NSYM
c                 ISYMC  = MULD2H(ISYMK,ISYMCK)
c                 ISYMCI = MULD2H(ISYMK,ISYCKI)
c                 DO K = 1,NRHF(ISYMK)
c                    NAK = NUMIA*(K - 1) + A
c                    DO C = 1,NVIR(ISYMC)
c                       KCKAI = IOFF5(ISYMI)
c    &                        + NT1AM(ISYMCK)*(NAI - 1)
c    &                        + IT1AM(ISYMC,ISYMK)
c    &                        + NVIR(ISYMC)*(K - 1) + C
c                       KCIAK = IOFF5(ISYMK)
c    &                        + NT1AM(ISYMCI)*(NAK - 1)
c    &                        + IT1AM(ISYMC,ISYMI)
c    &                        + NVIR(ISYMC)*(I - 1) + C
c                       TEST  = T2AM1(KCKAI) - T2AM2(KCIAK)
c                       NTST  = NTST + 1
c                       IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*)
c    & '   ISCLT: ISYMC,ISYMK,ISYMA,ISYMI: ',
c    & ISYMC,ISYMK,ISYMA,ISYMI
c     write(LUPRI,*)
c    & '   ISCLT: C,K,A,I: ',C,K,IA1+A-1,I
c     write(LUPRI,*)
c    & '   ISCLT: Diff.  : ',TEST
c     ICOUNE = ICOUNE + 1
c                       ENDIF
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     ndim   = NCKI(ISYCKI)*NUMIA
c     t1norm = dsqrt(ddot(ndim,T2AM1,1,T2AM1,1))
c     t2norm = dsqrt(ddot(ndim,T2AM2,1,T2AM2,1))
c     diff   = t2norm - t1norm
c      if (dabs(diff) .gt. 1.0d-15) then
c        write(LUPRI,*) '   ISCLT: Norm of 1T(ck,#ai) [output]: ',t1norm
c        write(LUPRI,*) '   ISCLT: Norm of 2T(ck,#ai) [output]: ',t2norm
c      endif
c     write(LUPRI,*) '   ISCLT: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested'
c     write(LUPRI,*) '   ISCLT: (ndim = ',ndim,')'
C
      RETURN
      END
C  /* Deck ccho_iterm */
      SUBROUTINE CCHO_ITERM(XJ1,XJ2,T3,T4,XK,T1,T2,VMAT,WMAT,QMAT,RMAT,
     &                      CORI,NUMIA,NUMIB,ISYMA,ISYMB,IOFF2,IOFF3,
     &                      IOFF4,IOFF5)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Calculate contribution to I term from scaled quantities:
C
C       CORI = CORI + Q(ai,bj)*V(aj,bi) + R(ai,bj)*W(aj,bi)
C
C     where
C
C       Q(ai,bj) =   Sum(ck) d(cki) * [ 2 * t(ai,ck) - t(ak,ci) ]
C
C                                   * [ 2 * t(ck,bj) - t(cj,bk) ]
C
C       R(ai,bj) =   Sum(ck) d(cki) * t(ck,ai) * [ 2 * t(ck,bj) - t(cj,bk) ]
C
C                  + Sum(ck) d(cki) * t(ci,ak) * [ 2 * t(cj,bk) - t(ck,bj) ]
C
C       V(aj,bi) = - Sum(de) d(dej) (dj|ea) * (ie|bd)
C
C       W(aj,bi) =   Sum(de) d(dej) (dj|ea) * (id|be)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XJ1(*),XJ2(*),T3(*),T4(*),XK(*),T1(*),T2(*)
      DIMENSION VMAT(*),WMAT(*),QMAT(*),RMAT(*)
      INTEGER IOFF2(8),IOFF3(8),IOFF4(8),IOFF5(8)
      INTEGER BI,BJ
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (XMONE = -1.00D0, ZERO = 0.00D0)
      PARAMETER (ONE   =  1.00D0, TWO  = 2.00D0)
C
c     lenj1  = NCKATR(ISYMB)*NUMIB
c     lenj2  = lenj1
c     lenkt  = NCKATR(ISYMA)*NUMIA
c     lent1  = NCKI(ISYMA)*NUMIA
c     lent2  = lent1
c     lent3  = NCKI(ISYMB)*NUMIB
c     lent4  = lent3
c     xj1nrm = dsqrt(ddot(lenj1,XJ1,1,XJ1,1))
c     xj2nrm = dsqrt(ddot(lenj2,XJ2,1,XJ2,1))
c     xktnrm = dsqrt(ddot(lenkt,XK,1,XK,1))
c     xt1nrm = dsqrt(ddot(lent1,T1,1,T1,1))
c     xt2nrm = dsqrt(ddot(lent2,T2,1,T2,1))
c     xt3nrm = dsqrt(ddot(lent3,T3,1,T3,1))
c     xt4nrm = dsqrt(ddot(lent4,T4,1,T4,1))
c     write(LUPRI,*)
c     write(LUPRI,*) '   ITERM: NUMIB, ISYMB: ',NUMIB,ISYMB
c     write(LUPRI,*) '   ITERM: NUMIA, ISYMA: ',NUMIA,ISYMA
c     write(LUPRI,*) '   ITERM: IOFF2: ',(IOFF2(ISYM),ISYM=1,NSYM)
c     write(LUPRI,*) '   ITERM: IOFF3: ',(IOFF3(ISYM),ISYM=1,NSYM)
c     write(LUPRI,*) '   ITERM: IOFF4: ',(IOFF4(ISYM),ISYM=1,NSYM)
c     write(LUPRI,*) '   ITERM: IOFF5: ',(IOFF5(ISYM),ISYM=1,NSYM)
c     write(LUPRI,*) '   ITERM: Norm of  J1: ',xj1nrm
c     write(LUPRI,*) '   ITERM: Norm of  J2: ',xj2nrm
c     write(LUPRI,*) '   ITERM: Norm of nKT: ',xktnrm
c     write(LUPRI,*) '   ITERM: Norm of nT1: ',xt1nrm
c     write(LUPRI,*) '   ITERM: Norm of nT2: ',xt2nrm
c     write(LUPRI,*) '   ITERM: Norm of  T3: ',xt3nrm
c     write(LUPRI,*) '   ITERM: Norm of  T4: ',xt4nrm
c     write(LUPRI,*)
C
      ISYMIJ = MULD2H(ISYMA,ISYMB)
C
      DO ISYMJ = 1,NSYM
C
         ISYMI = MULD2H(ISYMJ,ISYMIJ)
C
         ISYMCK = MULD2H(ISYMI,ISYMA)
         ISYMDE = MULD2H(ISYMI,ISYMB)
C
         NI  = NRHF(ISYMI)
         NJ  = NRHF(ISYMJ)
         NCK = NT1AM(ISYMCK)
         NDE = NMATAB(ISYMDE)
C
         NTOTCK = MAX(NCK,1)
         NTOTDE = MAX(NDE,1)
C
         IF ((NI.GT.0) .AND. (NJ.GT.0)) THEN
C
c      write(LUPRI,*) '   ITERM: ISYMA,ISYMI,ISYMB,ISYMJ: ',
c     &                      ISYMA,ISYMI,ISYMB,ISYMJ
c      write(LUPRI,*) '   ITERM: Space for v,w,q,r: ',NI*NJ*NUMIA*NUMIB
C
            NAI = NUMIA*NI
            NAJ = NUMIA*NJ
            NBI = NUMIB*NI
            NBJ = NUMIB*NJ
C
C-----------------------------
C           Calculate V and W.
C-----------------------------
C
            KOFFK = IOFF4(ISYMJ) + 1
            KOFFJ = IOFF2(ISYMI) + 1
C
            CALL DGEMM('T','N',NAJ,NBI,NDE,
     &                 XMONE,XK(KOFFK),NTOTDE,XJ1(KOFFJ),NTOTDE,
     &                 ZERO,VMAT,NAJ)
C
            CALL DGEMM('T','N',NAJ,NBI,NDE,
     &                 ONE,XK(KOFFK),NTOTDE,XJ2(KOFFJ),NTOTDE,
     &                 ZERO,WMAT,NAJ)
C
C-----------------------------
C           Calculate Q and R.
C-----------------------------
C
            KOFT1 = IOFF5(ISYMI) + 1
            KOFT3 = IOFF3(ISYMJ) + 1
C
            CALL DGEMM('T','N',NAI,NBJ,NCK,
     &                 ONE,T1(KOFT1),NTOTCK,T3(KOFT3),NTOTCK,
     &                 ZERO,QMAT,NAI)
C
            CALL DCOPY(NAI*NBJ,QMAT,1,RMAT,1)
C
            KOFT2 = KOFT1
C
            CALL DGEMM('T','N',NAI,NBJ,NCK,
     &                 XMONE,T2(KOFT2),NTOTCK,T3(KOFT3),NTOTCK,
     &                 TWO,QMAT,NAI)
C
            KOFT4 = KOFT3
C
            CALL DGEMM('T','N',NAI,NBJ,NCK,
     &                 ONE,T2(KOFT2),NTOTCK,T4(KOFT4),NTOTCK,
     &                 ONE,RMAT,NAI)
C
C-----------------------------------------
C           Calculate energy contribution.
C-----------------------------------------
C
            DO J = 1,NJ
               DO B = 1,NUMIB
C
                  BJ = NUMIB*(J - 1) + B
C
                  DO I = 1,NI
C
                     BI = NUMIB*(I - 1) + B
C
                     KOFFQ = NAI*(BJ - 1) + NUMIA*(I - 1) + 1
                     KOFFV = NAJ*(BI - 1) + NUMIA*(J - 1) + 1
                     KOFFR = KOFFQ
                     KOFFW = KOFFV
C
                     EI1  = DDOT(NUMIA,QMAT(KOFFQ),1,VMAT(KOFFV),1)
                     EI2  = DDOT(NUMIA,RMAT(KOFFR),1,WMAT(KOFFW),1)
                     CORI = CORI + EI1 + EI2
C
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDIF
C
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_gterm */
      SUBROUTINE CCHO_GTERM(XIINT,XJINT,T2VO,T1AM,XIAJB,FOCKD,NUMCHO,
     &                      CHOELE,WORK,LWORK,EG,IB1,ISYMB,NUMIB,
     &                      FBATCH,PRINT,NONI)
c    &                      EYS,EXR,EWQ,EVP)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Calculate the G term:
C
C        EG = EG + Sum(aibj) [ V(aj,bi) * P(ai,bj) + W(aj,bi) * Q(ai,bj) ]
C
C                + Sum(bjik) X(bj,ik) * R(bi,kj)
C
C                + Sum(bi)   Y(bi) * S(bi)
C
C     where
C
C        V(aj,bi) =   Sum(d)   d(djb) * (dj|ba) * t1(di)
C
C        W(aj,bi) =   Sum(d)   d(djb) * (bj|da) * t1(di)
C
C        X(bj,ik) =   Sum(da)  d(dja) * (dj|ab) * L(iakd)
C
C        Y(bi)    =   Sum(dja) d(dja) * (dj|ab) * L(iajd)  =  Sum(j) X(bj,ij)
C
C        P(ai,bj) =   Sum(ck)  d(cki) * t(ck,ai) * L(kcjb)
C
C                   + Sum(ck)  d(cki) * t(ci,ak) * L(jckb)
C
C        Q(ai,bj) = - Sum(ck)  d(cki) * s(ck,ai) * L(kcjb)
C
C        R(bi,kj) =   Sum(c)   d(cki) * t(bi,ck) * t1(cj)
C
C        S(bi)    = - Sum(ck)  d(cki) * s(ck,bi) * t1(ck)
C
C     and s(ai,bj) = 2 * t(ai,bj) - t(aj,bi), d(aij) and d(aib) denote
C     the occupied and virtual parts of the Cholesky decomposition of the
C     orbital energy denominator, and L(iajb) = 2 * (ia|jb) - (ja|ib).
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XIINT(*),XJINT(*),T2VO(*),T1AM(*),XIAJB(*),FOCKD(*)
      DIMENSION CHOELE(*)
      DIMENSION WORK(LWORK)
      LOGICAL   FBATCH,PRINT
      INTEGER IOFF1(8),IOFF2(8),IOFF3(8),IOFF4(8),IOFF5(8),IOFF6(8)
      INTEGER IOFT1(8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.00D0, ZERO = 0.00D0)
      PARAMETER (ONE   =  1.00D0, TWO  = 2.00D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C 
      TLAST = SECOND()
C 
      IF (PRINT) THEN
         TIMT = SECOND()
         WRITE(LUPRI,'(6X,A,/,6X,A,/)')
     &   'Calculation of the G term:',
     &   '=========================='
      ENDIF

c     leni   = NCKATR(ISYMB)*NUMIB
c     lenj   = leni
c     lent2  = NT2SQ(1)
c     lent1  = NT1AM(1)
c     liajb  = NT2AM(1)
c     xinorm = dnorm2(leni,XIINT,1)
c     xjnorm = dnorm2(lenj,XJINT,1)
c     xt2nrm = dnorm2(lent2,T2VO,1)
c     xt1nrm = dnorm2(lent1,T1AM,1)
c     xniajb = dnorm2(liajb,XIAJB,1)
c     if (lwork .lt. nt2sq(1)) CALL QUIT(' crappy test in GTERM ')
c     call cc_t2sq(XIAJB,WORK,1)
c     xsiajb = dnorm2(nt2sq(1),WORK,1)
c     write(LUPRI,*)
c     write(LUPRI,*) '   Entering GTERM:'
c     write(LUPRI,*) '   ==============='
c     write(LUPRI,*) '   IB1,NUMIB,ISYMB   : ',IB1,NUMIB,ISYMB
c     write(LUPRI,*) '   Norm of I(dj,a;#b): ',xinorm
c     write(LUPRI,*) '   Norm of J(dj,a;#b): ',xjnorm
c     write(LUPRI,*) '   Norm of T2VO      : ',xt2nrm
c     write(LUPRI,*) '   Norm of T1AM      : ',xt1nrm
c     write(LUPRI,*) '   Norm of L(iajb)   : ',xniajb,' (packed)'
c     write(LUPRI,*) '   Norm of L(iajb)   : ',xsiajb,' (squared)'
C
C----------------
C     T1 handles.
C----------------
C
      ICOUNT = 0
      DO ISYM = 1,NSYM
         IOFT1(ISYM) = ICOUNT
         ICOUNT = ICOUNT + NT1AM(ISYM)
      ENDDO
      NTOT1 = ICOUNT

c     write(LUPRI,*) '   GTERM: NTOT1: ',NTOT1
c     write(LUPRI,*) '   GTERM: IOFT1: ',(IOFT1(ISYM),ISYM=1,NSYM)
C
C----------------
C     Allocation.
C----------------
C
      KCHOL = 1
      KCHOV = KCHOL + NCKATR(ISYMB)
      KM1   = KCHOV + NTOT1*NUMIB
      KM2   = KM1   + NCKI(ISYMB)*NUMIB
      KJ1   = KM2   + NCKI(ISYMB)*NUMIB
      KJ2   = KJ1   + NCKATR(ISYMB)*NUMIB
      KI1   = KJ2   + NCKATR(ISYMB)*NUMIB
      KEND1 = KI1   + NCKATR(ISYMB)*NUMIB
      LWRK1 = LWORK - KEND1 + 1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Insufficient memory in CCHO_GTERM'
         WRITE(LUPRI,*) 'Need (more than): ',KEND1-1
         WRITE(LUPRI,*) 'Available       : ',LWORK
         CALL QUIT(' Insufficient memory in CCHO_GTERM ')
      ENDIF
C
C--------------------------------------------------
C     Extract integrals:
C     M1(ai,#bj) = L(ia,jb), M2(ai,#bj) = L(ja,ib).
C     Order according to IOFF3.
C--------------------------------------------------
C
      CALL CCHO_GXM12(XIAJB,WORK(KM1),WORK(KM2),IB1,NUMIB,ISYMB,IOFF3)
C
      DO ISYMI = 1,NSYM
C
         IF (NRHF(ISYMI) .LE. 0) GOTO 1000
C
         ISYMBI = MULD2H(ISYMI,ISYMB)
C
C----------------------------------------
C        Find largest a,j symmetry block.
C        Find largest k,j symmetry block.
C----------------------------------------
C
         MAXAJ = -1
         MAXKJ = -1
         DO ISYMJ = 1,NSYM
            ISYMA = MULD2H(ISYMJ,ISYMBI)
            MAXAJ = MAX(MAXAJ,NVIR(ISYMA)*NRHF(ISYMJ))
            MAXKJ = MAX(MAXKJ,NRHF(ISYMA)*NRHF(ISYMJ))
         ENDDO

c     write(LUPRI,*) '   GTERM: ISYMI,MAXAJ,MAXKJ: ',ISYMI,MAXAJ,MAXKJ
C
C---------------------
C        Set up batch.
C---------------------
C
         MINMEM = NTOT1 + 4*NCKATR(ISYMI) + 3*NT1AM(ISYMBI)*NUMIB
     &          + 2*MAXAJ*NUMIB + 2*MAXKJ*NUMIB
         IF (ISYMI .EQ. ISYMB) MINMEM = MINMEM + 2*NUMIB
         IF (FBATCH) THEN
            NEFI  = MIN(NONI,NRHF(ISYMI))
            LEFF  = NEFI*MINMEM + 1
            LWRKB = MIN(LWRK1,LEFF)
         ELSE
            LWRKB = LWRK1
         ENDIF
         NUMI = MIN(LWRKB/MINMEM,NRHF(ISYMI))
         IF (NUMI .LE. 0) THEN
            WRITE(LUPRI,*) 'Insufficient memory in CCHO_GTERM'
            WRITE(LUPRI,*) 'Minimum memory required: ',KEND1+MINMEM-1
            WRITE(LUPRI,*) 'Available in total     : ',LWORK
            WRITE(LUPRI,*) 'IB1,NUMIB,ISYMB        : ',IB1,NUMIB,ISYMB
            WRITE(LUPRI,*) 'Number of occ. and sym.: ',NRHF(ISYMI),ISYMI
            WRITE(LUPRI,*) 'Min. mem. for batch    : ',MINMEM
            WRITE(LUPRI,*) 'Available for batch    : ',LWRKB
            CALL QUIT(' Insufficient memory in CCHO_GTERM ')
         ENDIF
C
         NBATCH = (NRHF(ISYMI) - 1)/NUMI + 1
C
         IF (PRINT) THEN
            WRITE(LUPRI,'(9X,A,I1,A,/,9X,A)')
     &      'Batch over I, symmetry ',ISYMI,':',
     &      '-------------------------'
            WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10)')
     &      'Minimum work space required   : ',MINMEM,
     &      'Work space available for batch: ',LWRK1
            WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10,/)')
     &      'Number of occupied orbitals   : ',NRHF(ISYMI),
     &      'Required number of I-batches  : ',NBATCH
         ENDIF
C
         DO IBATCH = 1,NBATCH
C
            NUMII = NUMI
            IF (IBATCH .EQ. NBATCH) THEN
               NUMII = NRHF(ISYMI) - NUMI*(NBATCH - 1)
            ENDIF
C
            II1 = NUMI*(IBATCH - 1) + 1
C
            IF (PRINT) THEN
               WRITE(LUPRI,'(12X,A,I10,A,/,12X,A)')
     &         'I-batch number ',IBATCH,':',
     &         '--------------------------'
               WRITE(LUPRI,'(12X,A,I10,1X,I10,/)')
     &         'First and last I: ',II1,II1+NUMII-1
            ENDIF
C
C----------------------
C           Allocation.
C----------------------
C
            KCHOO = KEND1
            KM3   = KCHOO + NTOT1*NUMII
            KT1   = KM3   + NCKATR(ISYMI)*NUMII
            KT2   = KT1   + NCKATR(ISYMI)*NUMII
            KT3   = KT2   + NCKATR(ISYMI)*NUMII
            KT4   = KT3   + NCKATR(ISYMI)*NUMII
            KVMAT = KT4   + NT1AM(ISYMBI)*NUMIB*NUMII
            KWMAT = KVMAT + NT1AM(ISYMBI)*NUMIB*NUMII
            KPMAT = KWMAT + NT1AM(ISYMBI)*NUMIB*NUMII
            KQMAT = KPMAT + MAXAJ*NUMIB*NUMII
            KXMAT = KQMAT + MAXAJ*NUMIB*NUMII
            KRMAT = KXMAT + MAXKJ*NUMIB*NUMII
            KEND2 = KRMAT + MAXKJ*NUMIB*NUMII
            IF (ISYMI .EQ. ISYMB) THEN
               KYMAT = KEND2
               KSMAT = KYMAT + NUMIB*NUMII
               KEND2 = KSMAT + NUMIB*NUMII
            ENDIF
            LWRK2 = LWORK - KEND2 + 1
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Batching bug in CCHO_GTERM'
               CALL QUIT(' Error in CCHO_GTERM ')
            ENDIF
C
C-----------------------------------------------------------
C           Extract:
C           M3(da,#ik) = L(ia,kd). Order according to IOFF4.
C-----------------------------------------------------------
C
            CALL CCHO_GXM3(XIAJB,WORK(KM3),II1,NUMII,ISYMI,IOFF4)
C
C---------------------------------------------------
C           Extract:
C           T1(ck,a#i)  = T(ck,ai), ordering: IOFF5.
C           T2(ck,a#i)  = T(ci,ak), ordering: IOFF5.
C---------------------------------------------------
C
            CALL CCHO_GXT12(T2VO,WORK(KT1),WORK(KT2),
     &                      II1,NUMII,ISYMI,IOFF5)
C
C---------------------------------------------------------------
C           Construct 2CME amplitudes:
C           T3(ck,a#i) = 2*T(ck,ai) - T(ci,ak), ordering: IOFF5.
C---------------------------------------------------------------
C
            LENTOT = NCKATR(ISYMI)*NUMII
            CALL DCOPY(LENTOT,WORK(KT1),1,WORK(KT3),1)
            CALL DSCAL(LENTOT,TWO,WORK(KT3),1)
            CALL DAXPY(LENTOT,XMONE,WORK(KT2),1,WORK(KT3),1)
C
C---------------------------------------------------
C           Extract:
C           T4(#b#i,kc) = T(ck,bi), ordering: IOFF6.
C---------------------------------------------------
C
            CALL CCHO_GXT4(T2VO,WORK(KT4),II1,NUMII,ISYMI,
     &                     IB1,NUMIB,ISYMB,IOFF6)

c     write(LUPRI,*)
c     write(LUPRI,*) '   GTERM: calling _gt to check unscaled amp.:'
c     call chk_gt(WORK(KT1),WORK(KT4),II1,NUMII,ISYMI,
c    &            IB1,NUMIB,ISYMB,IOFF5,IOFF6)
c     write(LUPRI,*)
C
C-------------------------------
C           Start Cholesky loop.
C-------------------------------
C
c     write(LUPRI,*) '   GTERM: ****NOTICE**** putting NUMCHO = 1'
c           do icho = 1,1
            DO ICHO = 1,NUMCHO
C
C--------------------------
C              Inititalize.
C--------------------------
C
               IF (ISYMI .EQ. ISYMB) THEN
                  CALL DZERO(WORK(KYMAT),NUMIB*NUMII)
               ENDIF
C
               GCOR = ZERO
C
C----------------------------------------------------------------------------
C              Get Cholesky vectors.
C              d(dj,a) : actual vector, in KCHOL; all dj,a of symmetry=ISYMB.
C              d(dj,#b): actual vector, in KCHOV; all dj for #b.
C              d(ck,#i): update vector, in KCHOO; all ck for #i.
C----------------------------------------------------------------------------
C
               ISYDJA = ISYMB
               CALL CCHO_DECHO7(FOCKD,CHOELE,NUMCHO,ICHO,WORK(KCHOL),
     &                          WORK(KCHOO),WORK(KCHOV),ISYDJA,
     &                          IB1,NUMIB,ISYMB,II1,NUMII,ISYMI)

c     write(LUPRI,*) '   GTERM: ****NOTICE**** putting vectors = 1'
c     do idum = 1,nckatr(isymb)
c        KOFF = KCHOL + idum - 1
c        WORK(KOFF) = 1.0D0
c     enddo
c     do idum = 1,NTOT1*NUMIB
c        KOFF = KCHOV + idum - 1
c        WORK(KOFF) = 1.0D0
c     enddo
c     do idum = 1,NTOT1*NUMII
c        KOFF = KCHOO + idum - 1
c        WORK(KOFF) = 1.0D0
c     enddo

c     call chk_vec2(FOCKD,CHOELE,NUMCHO,ICHO,WORK(KCHOL),
c    &              WORK(KCHOO),WORK(KCHOV),ISYDJA,
c    &              IB1,NUMIB,ISYMB,II1,NUMII,ISYMI,WORK(KEND2),LWRK2)
C
C----------------------------------------------------------------------
C              Extract virtual integrals and scale with actual vectors:
C              J1(aj,#bd) = d(dj,#b) * J(dj,a;#b), ordering: IOFF1
C              I1(aj,#bd) = d(dj,#b) * I(dj,a;#b), ordering: IOFF1
C              J2(da,#bj) = d(dj,a)  * J(dj,a;#b), ordering: IOFF2
C----------------------------------------------------------------------
C
               CALL CCHO_GXSJ12I1(XJINT,XIINT,
     &                            WORK(KJ1),WORK(KJ2),WORK(KI1),
     &                            WORK(KCHOV),WORK(KCHOL),
     &                            IB1,NUMIB,ISYMB,
     &                            IOFF1,IOFF2,IOFT1,NTOT1)
C
C-----------------------------------------------------------------------
C              Scale amplitudes with update vectors.
C              T1(ck,a#i)   <- d(ck,#i) * T1(ck,a#i),   ordering: IOFF5.
C              T2(ck,a#i)   <- d(ck,#i) * T2(ck,a#i),   ordering: IOFF5.
C              T3(ck,a#i)   <- d(ck,#i) * T3(ck,a#i),   ordering: IOFF5.
C              T4(#b#i,k,c) <- d(ck,#i) * T4(#b#i,k,c), ordering: IOFF6.
C-----------------------------------------------------------------------
C
               CALL CCHO_GST1234(WORK(KT1),WORK(KT2),WORK(KT3),
     &                           WORK(KT4),WORK(KCHOO),II1,NUMII,ISYMI,
     &                           IB1,NUMIB,ISYMB,IOFF5,IOFF6,IOFT1,
     &                           NTOT1)

c     lckai = NCKATR(ISYMI)*NUMII
c     if (lwrk2 .gt. lckai) then
c        call dcopy(lckai,WORK(KT1),1,WORK(KEND2),1)
c        call dscal(lckai,2.0D0,WORK(KEND2),1)
c        call daxpy(lckai,-1.0D0,WORK(KT2),1,WORK(KEND2),1)
c        call daxpy(lckai,-1.0D0,WORK(KT3),1,WORK(KEND2),1)
c        diff = dnorm2(lckai,WORK(KEND2),1)
c        write(LUPRI,*) '   GTERM: (2*T1-T2)-T3: ',diff
c     endif

c     write(LUPRI,*)
c     write(LUPRI,*) '   GTERM: calling _gt to check scaled amp., ICHO = ',
c    &           ICHO,':'
c     call chk_gt(WORK(KT1),WORK(KT4),II1,NUMII,ISYMI,
c    &            IB1,NUMIB,ISYMB,IOFF5,IOFF6)
c     write(LUPRI,*)

c     write(LUPRI,*) '   GTERM: IOFF1: ',(IOFF1(JSYM),JSYM=1,NSYM)
c     write(LUPRI,*) '   GTERM: IOFF2: ',(IOFF2(JSYM),JSYM=1,NSYM)
c     write(LUPRI,*) '   GTERM: IOFF3: ',(IOFF3(JSYM),JSYM=1,NSYM)
c     write(LUPRI,*) '   GTERM: IOFF4: ',(IOFF4(JSYM),JSYM=1,NSYM)
c     write(LUPRI,*) '   GTERM: IOFF5: ',(IOFF5(JSYM),JSYM=1,NSYM)
c     write(LUPRI,*) '   GTERM: IOFF6: ',(IOFF6(JSYM),JSYM=1,NSYM)
c     write(LUPRI,*) '   GTERM: IOFT1: ',(IOFT1(JSYM),JSYM=1,NSYM)
C
C--------------------------------
C              Calculate V and W.
C--------------------------------
C
               ISYMD  = ISYMI
               ISYMAJ = MULD2H(ISYMD,ISYMB)
C
               NAJB  = NT1AM(ISYMAJ)*NUMIB
               NTAJB = MAX(NAJB,1)
               ND    = NVIR(ISYMD)
               NTOTD = MAX(ND,1)
C
               KOFFJ = KJ1 + IOFF1(ISYMD)
               KOFFT = IT1AM(ISYMD,ISYMI) + NVIR(ISYMD)*(II1 - 1) + 1
C
               CALL DGEMM('N','N',NAJB,NUMII,ND,
     &                    ONE,WORK(KOFFJ),NTAJB,T1AM(KOFFT),NTOTD,
     &                    ZERO,WORK(KVMAT),NTAJB)
C
               KOFFI = KI1 + IOFF1(ISYMD)
C
               CALL DGEMM('N','N',NAJB,NUMII,ND,
     &                    ONE,WORK(KOFFI),NTAJB,T1AM(KOFFT),NTOTD,
     &                    ZERO,WORK(KWMAT),NTAJB)
C
C-------------------------------------
C              Loop over j symmetries.
C-------------------------------------
C
               DO ISYMJ = 1,NSYM
C
                  IF (NRHF(ISYMJ) .LE. 0) GOTO 998
C
C-----------------------------------
C                 Calculate P and Q.
C-----------------------------------
C
                  ISYMCK = MULD2H(ISYMJ,ISYMB)
                  ISYMAI = ISYMCK
                  ISYMA  = MULD2H(ISYMAI,ISYMI)
C
                  NCK   = NT1AM(ISYMCK)
                  NTOCK = MAX(NCK,1)
                  NAI   = NVIR(ISYMA)*NUMII
                  NTOAI = MAX(NAI,1)
                  NBJ   = NUMIB*NRHF(ISYMJ)
C
                  KOFT1 = KT1 + IOFF5(ISYMA)
                  KOFT2 = KT2 + IOFF5(ISYMA)
                  KOFT3 = KT3 + IOFF5(ISYMA)
                  KOFM1 = KM1 + IOFF3(ISYMJ)
                  KOFM2 = KM2 + IOFF3(ISYMJ)
C
                  CALL DGEMM('T','N',NAI,NBJ,NCK,
     &                       ONE,WORK(KOFT1),NTOCK,WORK(KOFM1),NTOCK,
     &                       ZERO,WORK(KPMAT),NTOAI)
C
                  CALL DGEMM('T','N',NAI,NBJ,NCK,
     &                       XMONE,WORK(KOFT3),NTOCK,WORK(KOFM1),NTOCK,
     &                       ZERO,WORK(KQMAT),NTOAI)
C
                  CALL DGEMM('T','N',NAI,NBJ,NCK,
     &                       ONE,WORK(KOFT2),NTOCK,WORK(KOFM2),NTOCK,
     &                       ONE,WORK(KPMAT),NTOAI)
C
C-----------------------------------
C                 Calculate X and R.
C-----------------------------------
C
                  ISYMDA = MULD2H(ISYMJ,ISYMB)
                  ISYMIK = ISYMDA
                  ISYMK  = MULD2H(ISYMIK,ISYMI)
                  ISYMC  = ISYMJ
C
                  NDA   = NMATAB(ISYMDA)
                  NTODA = MAX(NDA,1)
                  NBJ   = NUMIB*NRHF(ISYMJ)
                  NTOBJ = MAX(NBJ,1)
                  NIK   = NUMII*NRHF(ISYMK)
C
                  NBIK  = NUMIB*NUMII*NRHF(ISYMK)
                  NTBIK = MAX(NBIK,1)
                  NC    = NVIR(ISYMC)
                  NTOTC = MAX(NC,1)
                  NJ    = NRHF(ISYMJ)
C
                  KOFJ2 = KJ2 + IOFF2(ISYMJ)
                  KOFM3 = KM3 + IOFF4(ISYMK)
                  KOFT4 = KT4 + IOFF6(ISYMC)
                  KOFFT = IT1AM(ISYMC,ISYMJ) + 1
C
                  CALL DGEMM('T','N',NBJ,NIK,NDA,
     &                       ONE,WORK(KOFJ2),NTODA,WORK(KOFM3),NTODA,
     &                       ZERO,WORK(KXMAT),NTOBJ)
C
                  CALL DGEMM('N','N',NBIK,NJ,NC,
     &                       ONE,WORK(KOFT4),NTBIK,T1AM(KOFFT),NTOTC,
     &                       ZERO,WORK(KRMAT),NTBIK)

c     lenx = NBJ*NIK
c     lenr = NBIK*NJ
c     lmax = NUMIB*NUMII*MAXKJ
c     if (lenx .NE. lenr) then
c        write(LUPRI,*) '   GTERM: lenx,lenr: ',lenx,lenr
c        CALL QUIT(' lenx .NE. lenr ')
c     endif
c     if (lenx .GT. lmax) then
c        write(LUPRI,*) '   GTERM: lenx,lmax: ',lenx,lmax
c        CALL QUIT(' lenx .GT. lmax ')
c     endif
c     xnorm = dnorm2(lenx,WORK(KXMAT),1)
c     rnorm = dnorm2(lenr,WORK(KRMAT),1)
c     write(LUPRI,*)
c    & '   GTERM: After XR cal.: ISYMB,ISYMJ,ISYMI,ISYMK,ISYMC: ',
c    &                           ISYMB,ISYMJ,ISYMI,ISYMK,ISYMC
c     write(LUPRI,*) '   GTERM: After XR cal.: Norm of X(#bj,#ik): ',xnorm
c     write(LUPRI,*) '   GTERM: After XR cal.: Norm of R(#b#i,kj): ',rnorm
c     write(LUPRI,*) '   GTERM: XMAT(bj,ik):'
c     call output(WORK(KXMAT),1,NBJ,1,NIK,NBJ,NIK,1,6)
c     write(LUPRI,*) '   GTERM: RMAT(bik,j):'
c     call output(WORK(KRMAT),1,NBIK,1,NJ,NBIK,NJ,1,6)
C
C-----------------------------
C                 Calculate Y.
C-----------------------------
C
                  IF (ISYMI .EQ. ISYMB) THEN
C
                     NBJ = NUMIB*NRHF(ISYMJ)
C
                     DO J = 1,NRHF(ISYMJ)
                        DO I = 1,NUMII
C
                           IJ = NUMII*(J - 1) + I
C
                           DO B = 1,NUMIB
C
                              KBI = NUMIB*(I - 1) + B
                              KBJ = NUMIB*(J - 1) + B
C
                              KOFFY = KYMAT + KBI - 1
                              KOFFX = KXMAT + NBJ*(IJ - 1) + KBJ - 1
C
                              WORK(KOFFY) = WORK(KOFFY) + WORK(KOFFX)
C
                           ENDDO
C
                        ENDDO
                     ENDDO
C
                  ENDIF
C
C------------------------------------------------
C                 Calculate energy contributions:
C                 GCOR = GCOR + V*P + W*Q
C------------------------------------------------
C
                  ISYMAJ = MULD2H(ISYMI,ISYMB)
                  ISYMA  = MULD2H(ISYMJ,ISYMAJ)
C
                  NAJ = NT1AM(ISYMAJ)
                  NAI = NVIR(ISYMA)*NUMII
C
                  DO I = 1,NUMII
C
                     IOFAI = NVIR(ISYMA)*(I - 1)
C
                     DO B = 1,NUMIB
C
                        KBI = NUMIB*(I - 1) + B
C
                        DO J = 1,NRHF(ISYMJ)
C
                           IOFAJ = IT1AM(ISYMA,ISYMJ)
     &                           + NVIR(ISYMA)*(J - 1)
                           KBJ   = NUMIB*(J - 1) + B
C
                           KOFFV = KVMAT + NAJ*(KBI - 1) + IOFAJ
                           KOFFW = KWMAT + NAJ*(KBI - 1) + IOFAJ
                           KOFFP = KPMAT + NAI*(KBJ - 1) + IOFAI
                           KOFFQ = KQMAT + NAI*(KBJ - 1) + IOFAI
C
                           GCOR = GCOR
     &                          + DDOT(NVIR(ISYMA),WORK(KOFFV),1,
     &                                             WORK(KOFFP),1)
     &                          + DDOT(NVIR(ISYMA),WORK(KOFFW),1,
     &                                             WORK(KOFFQ),1)
C
      EVP = EVP + DDOT(NVIR(ISYMA),WORK(KOFFV),1,WORK(KOFFP),1)
      EWQ = EWQ + DDOT(NVIR(ISYMA),WORK(KOFFW),1,WORK(KOFFQ),1)
C
                        ENDDO
C
                     ENDDO
C
                  ENDDO
C
C----------------------------------------------
C                 Calculate energy contribution:
C                 GCOR = GCOR + X*R
C----------------------------------------------
C
                  ISYMKJ = MULD2H(ISYMI,ISYMB)
                  ISYMK  = MULD2H(ISYMJ,ISYMKJ)
C
                  NBJ = NUMIB*NRHF(ISYMJ)
                  NBI = NUMIB*NUMII

c     lenx1 = NBJ*NUMII*NRHF(ISYMK)
c     lenr1 = NBI*NRHF(ISYMK)*NRHF(ISYMJ)
c     if (lenx1 .NE. lenr1) then
c        write(LUPRI,*) '   GTERM: lenx1,lenr1: ',lenx1,lenr1
c        CALL QUIT(' lenx1 .NE. lenr1 ')
c     endif
c     xnorm = dnorm2(lenx1,WORK(KXMAT),1)
c     rnorm = dnorm2(lenr1,WORK(KRMAT),1)
c     write(LUPRI,*) '   GTERM: Bef. en. cal.: ISYMB,ISYMJ,ISYMI,ISYMK: ',
c    &                                     ISYMB,ISYMJ,ISYMI,ISYMK
c     write(LUPRI,*) '   GTERM: Bef. en. cal.: Norm of X(#bj,#ik): ',xnorm
c     write(LUPRI,*) '   GTERM: Bef. en. cal.: Norm of R(#b#i,kj): ',rnorm
c     write(LUPRI,*) '   GTERM: XMAT(bj,ik):'
c     NIK = NUMII*NRHF(ISYMK)
c     call output(WORK(KXMAT),1,NBJ,1,NIK,NBJ,NIK,1,6)
c     write(LUPRI,*) '   GTERM: RMAT(bik,j):'
c     NBIK = NBI*NRHF(ISYMK)
c     NJ   = NRHF(ISYMJ)
c     call output(WORK(KRMAT),1,NBIK,1,NJ,NBIK,NJ,1,6)
c     write(LUPRI,*) '   GTERM: Before calc. EXR = ',EXR
C
                  DO K = 1,NRHF(ISYMK)
                     DO I = 1,NUMII
                        DO J = 1,NRHF(ISYMJ)
                           DO B = 1,NUMIB
C
                              KOFFX = KXMAT
     &                              + NUMIB*NRHF(ISYMJ)*NUMII*(K - 1)
     &                              + NUMIB*NRHF(ISYMJ)*(I - 1)
     &                              + NUMIB*(J - 1)
     &                              + B - 1
                              KOFFR = KRMAT
     &                              + NUMIB*NUMII*NRHF(ISYMK)*(J - 1)
     &                              + NUMIB*NUMII*(K - 1)
     &                              + NUMIB*(I - 1)
     &                              + B - 1

                              GCOR = GCOR + WORK(KOFFX)*WORK(KOFFR)

      EXR = EXR + WORK(KOFFX)*WORK(KOFFR)

                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
 
c     write(LUPRI,*) '   GTERM: After ISYMJ = ',ISYMJ,': EXR = ',EXR
C
  998             CONTINUE
C
               ENDDO
C
C---------------------------------------------------------------------
C              Calculate S and energy contribution: GCOR = GCOR + Y*S.
C---------------------------------------------------------------------
C
               IF (ISYMI .EQ. ISYMB) THEN
C
                  NCK   = NT1AM(1)
                  NTOCK = MAX(NCK,1)
C
                  DO I = 1,NUMII
C
                     IOFBI = NVIR(ISYMB)*(I - 1) + IB1
                     KOFT3 = KT3 + IOFF5(ISYMB) + NCK*(IOFBI - 1)
                     KOFFS = KSMAT + NUMIB*(I - 1)
C
                     CALL DGEMV('T',NCK,NUMIB,XMONE,WORK(KOFT3),NTOCK,
     &                          T1AM,1,ZERO,WORK(KOFFS),1)
C
                  ENDDO
C
                  GCOR = GCOR
     &                 + DDOT(NUMIB*NUMII,WORK(KYMAT),1,WORK(KSMAT),1)
C
      EYS = EYS + DDOT(NUMIB*NUMII,WORK(KYMAT),1,WORK(KSMAT),1)
C
               ENDIF
C
C--------------------------------------
C              Add energy contribution.
C--------------------------------------
C
               EG = EG + GCOR
C 
              ENERGG(ICHO) = ENERGG(ICHO) + GCOR
C 
C Decommented by Domenico
                    TNOW = SECOND()
                    DELTAT = TNOW - TLAST
                    TLAST = TNOW
                    SCNDSG(ICHO) = SCNDSG(ICHO) 
     &                          + DELTAT
C Decommented by Domenico
C 
C
c              IF (PRINT) THEN
c                 WRITE(LUPRI,'(15X,A,I3,A,/,15X,A)')
c    &            'Status after Cholesky vector',ICHO,':',
c    &            '--------------------------------'
c                 IF (ABS(GCOR) .LT. THRCHO) THEN
c                    WRITE(LUPRI,'(15X,A)') 'G term converged'
c                 ELSE
c                    WRITE(LUPRI,'(15X,A)') 'G term not converged'
c                 ENDIF
c                 TIM = SECOND() - TIMT
c                 WRITE(LUPRI,'(15X,A,F10.2,A,/)')
c    &            'Accumulated G-time: ',TIM,' seconds'
c              ENDIF

c     lckai = NCKATR(ISYMI)*NUMII
c     if (lwrk2 .gt. lckai) then
c        call dcopy(lckai,WORK(KT1),1,WORK(KEND2),1)
c        call dscal(lckai,2.0D0,WORK(KEND2),1)
c        call daxpy(lckai,-1.0D0,WORK(KT2),1,WORK(KEND2),1)
c        call daxpy(lckai,-1.0D0,WORK(KT3),1,WORK(KEND2),1)
c        diff = dnorm2(lckai,WORK(KEND2),1)
c        write(LUPRI,*) '   GTERM: end of cho loop:(2*T1-T2)-T3: ',diff
c     endif
C
C--------------------------------------------------
C              Skip remaining vectors if converged.
C--------------------------------------------------
C
               IF (DABS(GCOR) .LT. THRCHO) GOTO 999
C
            ENDDO
C
  999       CONTINUE
C
         ENDDO
C
 1000    CONTINUE
C
      ENDDO

c     leni   = NCKATR(ISYMB)*NUMIB
c     lenj   = leni
c     lent2  = NT2SQ(1)
c     lent1  = NT1AM(1)
c     liajb  = NT2AM(1)
c     xinorm = dnorm2(leni,XIINT,1)
c     xjnorm = dnorm2(lenj,XJINT,1)
c     xt2nrm = dnorm2(lent2,T2VO,1)
c     xt1nrm = dnorm2(lent1,T1AM,1)
c     xniajb = dnorm2(liajb,XIAJB,1)
c     if (lwork .lt. nt2sq(1)) CALL QUIT(' crappy test in GTERM ')
c     call cc_t2sq(XIAJB,WORK,1)
c     xsiajb = dnorm2(nt2sq(1),WORK,1)
c     write(LUPRI,*)
c     write(LUPRI,*) '   Exiting GTERM:'
c     write(LUPRI,*) '   =============='
c     write(LUPRI,*) '   IB1,NUMIB,ISYMB   : ',IB1,NUMIB,ISYMB
c     write(LUPRI,*) '   Norm of I(dj,a;#b): ',xinorm
c     write(LUPRI,*) '   Norm of J(dj,a;#b): ',xjnorm
c     write(LUPRI,*) '   Norm of T2VO      : ',xt2nrm
c     write(LUPRI,*) '   Norm of T1AM      : ',xt1nrm
c     write(LUPRI,*) '   Norm of L(iajb)   : ',xniajb,' (packed)'
c     write(LUPRI,*) '   Norm of L(iajb)   : ',xsiajb,' (squared)'
c     write(LUPRI,*)
C
      RETURN
      END
C  /* Deck ccho_gxm12 */
      SUBROUTINE CCHO_GXM12(XIAJB,XM1,XM2,IB1,NUMIB,ISYMB,IOFF3)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract integrals:
C
C        M1(ai,#bj) = L(ia,jb)
C        M2(ai,#bj) = L(ja,ib)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XIAJB(*),XM1(*),XM2(*)
      INTEGER IOFF3(8)
      INTEGER AI,BJ,AJ,BI,AIBJ,AJBI
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------------
C     Set up index array.
C------------------------
C
      ICOUNT = 0
C
      DO ISYMJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMJ,ISYMB)
C
         IOFF3(ISYMJ) = ICOUNT
C
         ICOUNT = ICOUNT + NT1AM(ISYMAI)*NUMIB*NRHF(ISYMJ)
C
      ENDDO
C
C-----------------------
C     Extract M1 and M2.
C-----------------------
C
      DO ISYMBJ = 1,NSYM
C
         ISYMAI = ISYMBJ
         ISYMJ  = MULD2H(ISYMBJ,ISYMB)
C
         DO J = 1,NRHF(ISYMJ)
            DO IB = 1,NUMIB
C
               B = IB1 + IB - 1
C
               IBJ = NUMIB*(J - 1) + IB
               BJ  = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
               DO ISYMI = 1,NSYM
C
                  ISYMA  = MULD2H(ISYMI,ISYMAI)
                  ISYMAJ = MULD2H(ISYMA,ISYMJ)
                  ISYMBI = MULD2H(ISYMI,ISYMB)
C
                  DO I = 1,NRHF(ISYMI)
C
                     BI = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I - 1) + B
C
                     DO A = 1,NVIR(ISYMA)
C
                        AI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1)
     &                     + A
                        AJ = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J - 1)
     &                     + A
C
                        AIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(AI,BJ)
                        AJBI = IT2AM(ISYMAJ,ISYMBI) + INDEX(AJ,BI)
C
                        KOFM = IOFF3(ISYMJ) + NT1AM(ISYMAI)*(IBJ - 1)
     &                       + AI
C
                        XM1(KOFM) = XIAJB(AIBJ)
                        XM2(KOFM) = XIAJB(AJBI)
C
                     ENDDO
C
                  ENDDO
C
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO
C
C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     DO ISYMJ = 1,NSYM
c        ISYMBJ = MULD2H(ISYMJ,ISYMB)
c        ISYMAI = ISYMBJ
c        DO J = 1,NRHF(ISYMJ)
c           DO B = 1,NUMIB
c              NBJ = NUMIB*(J - 1) + B
c              DO ISYMI = 1,NSYM
c                 ISYMA = MULD2H(ISYMI,ISYMAI)
c                 ISYMAJ = MULD2H(ISYMA,ISYMJ)
c                 ISYMBI = MULD2H(ISYMI,ISYMB)
c                 DO I = 1,NRHF(ISYMI)
c                    NBI = NUMIB*(I - 1) + B
c                    DO A = 1,NVIR(ISYMA)
c                       NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1)
c    &                      + A
c                       NAJ = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J - 1)
c    &                      + A
c                       KAIBJ = IOFF3(ISYMJ) + NT1AM(ISYMAI)*(NBJ - 1)
c    &                        + NAI
c                       KAJBI = IOFF3(ISYMI) + NT1AM(ISYMAJ)*(NBI - 1)
c    &                        + NAJ
c                       TEST  = XM1(KAIBJ) - XM2(KAJBI)
c                       NTST  = NTST + 1
c                       IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*) '      GXM12: ISYMA,ISYMI,ISYMB,ISYMJ: ',
c    &                         ISYMA,ISYMI,ISYMB,ISYMJ
c     write(LUPRI,*) '      GXM12: A,I,B,J: ',A,I,B,J
c     write(LUPRI,*) '      GXM12: Diff.  : ',TEST
c     ICOUNE = ICOUNE + 1
c                       ENDIF
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     ndim = NCKI(ISYMB)*NUMIB
c     write(LUPRI,*) '      GXM12: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested (dim = ',ndim,')'
c     xm1nrm = dnorm2(ndim,XM1,1)
c     xm2nrm = dnorm2(ndim,XM2,1)
c     write(LUPRI,*) '      GXM12: IB1,NUMIB,ISYMB   : ',IB1,NUMIB,ISYMB
c     write(LUPRI,*) '      GXM12: Norm of M1(ai,#bj): ',xm1nrm
c     write(LUPRI,*) '      GXM12: Norm of M2(ai,#bj): ',xm2nrm
C
      RETURN
      END
C  /* Deck ccho_gxm3 */
      SUBROUTINE CCHO_GXM3(XIAJB,XM3,II1,NUMII,ISYMI,IOFF4)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract integrals:
C
C        M3(da,#ik) = L(ia,kd)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XIAJB(*),XM3(*)
      INTEGER IOFF4(8)
      INTEGER AI,DK,DKAI,DA
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------------
C     Set up index array.
C------------------------
C
      ICOUNT = 0
C
      DO ISYMK = 1,NSYM
C
         ISYMDA = MULD2H(ISYMK,ISYMI)
C
         IOFF4(ISYMK) = ICOUNT
C
         ICOUNT = ICOUNT + NMATAB(ISYMDA)*NUMII*NRHF(ISYMK)
C
      ENDDO
C
C----------------
C     Extract M3.
C----------------
C
      DO ISYMAI = 1,NSYM
C
         ISYMA  = MULD2H(ISYMAI,ISYMI)
         ISYMDK = ISYMAI
C
         DO II = 1,NUMII
C
            I = II1 + II - 1
C
            DO A = 1,NVIR(ISYMA)
C
               AI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
               DO ISYMK = 1,NSYM
C
                  ISYMD  = MULD2H(ISYMK,ISYMDK)
                  ISYMDA = MULD2H(ISYMD,ISYMA)
C
                  DO K = 1,NRHF(ISYMK)
C
                     IK = NUMII*(K - 1) + II
C
                     DO D = 1,NVIR(ISYMD)
C
                        DK = IT1AM(ISYMD,ISYMK) + NVIR(ISYMD)*(K - 1)
     &                     + D
                        DA = IMATAB(ISYMD,ISYMA) + NVIR(ISYMD)*(A - 1)
     &                     + D
C
                        DKAI = IT2AM(ISYMDK,ISYMAI) + INDEX(DK,AI)
C
                        KOFM = IOFF4(ISYMK) + NMATAB(ISYMDA)*(IK - 1)
     &                       + DA
C
                        XM3(KOFM) = XIAJB(DKAI)
C
                     ENDDO
C
                  ENDDO
C
               ENDDO
C
            ENDDO
C
         ENDDO
C
      ENDDO
C
C DEBUG SECTION:
C***************
c     ndim   = NCKATR(ISYMI)*NUMII
c     xm3nrm = dnorm2(ndim,XM3,1)
c     write(LUPRI,*) '      GXM3: II1,NUMII,ISYMI   : ',II1,NUMII,ISYMI
c     write(LUPRI,*) '      GXM3: Norm of M3(da,#ik): ',xm3nrm
C
      RETURN
      END
C  /* Deck ccho_gxt12 */
      SUBROUTINE CCHO_GXT12(T2VO,T1,T2,II1,NUMII,ISYMI,IOFF5)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract amplitudes:
C
C        T1(ck,a#i)  = T(ck,ai), ordering: IOFF5.
C        T2(ck,a#i)  = T(ci,ak), ordering: IOFF5.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2VO(*),T1(*),T2(*)
      INTEGER IOFF5(8)
#include <ccorb.h>
#include <ccsdsym.h>
C
C------------------
C     Set up IOFF5.
C------------------
C
      ICOUNT = 0
C
      DO ISYMA = 1,NSYM
C
         ISYMCK = MULD2H(ISYMA,ISYMI)
C
         IOFF5(ISYMA) = ICOUNT
C
         ICOUNT = ICOUNT + NT1AM(ISYMCK)*NVIR(ISYMA)*NUMII
C
      ENDDO

c     if (ICOUNT .NE. NCKATR(ISYMI)*NUMII) then
c        write(LUPRI,*) 'WARNING: GXT12: ICOUNT, actual dim.: ',
c    &              ICOUNT,NCKATR(ISYMI)*NUMII
c     endif
c     write(LUPRI,*) '      GXT12: II1,NUMII,ISYMI: ',II1,NUMII,ISYMI
c     write(LUPRI,*) '      GXT12: IOFF5: ',(IOFF5(ISYM),ISYM=1,NSYM)
c     write(LUPRI,*) '      GXT12: Norm of T2VO: ',dnorm2(nt2sq(1),T2VO,1)
C
C-------------
C     Extract.
C-------------
C
      DO ISYMKI = 1,NSYM
C
         ISYMCA = ISYMKI
         ISYMK  = MULD2H(ISYMKI,ISYMI)
C
         DO II = 1,NUMII
C
            I = II1 + II - 1
C
            DO K = 1,NRHF(ISYMK)
C
               KI = IMATIJ(ISYMK,ISYMI) + NRHF(ISYMK)*(I - 1) + K
               IK = IMATIJ(ISYMI,ISYMK) + NRHF(ISYMI)*(K - 1) + I
C
               DO ISYMA = 1,NSYM
C
                  ISYMC  = MULD2H(ISYMA,ISYMCA)
                  ISYMCK = MULD2H(ISYMC,ISYMK)
C
                  NC = NVIR(ISYMC)
C
                  IF (NC .GT. 0) THEN
C
                     IOFCK = IT1AM(ISYMC,ISYMK) + NC*(K - 1) + 1
C
                     DO A = 1,NVIR(ISYMA)
C
                        AI = NVIR(ISYMA)*(II - 1) + A
C
                        IOFCA = IMATAB(ISYMC,ISYMA) + NC*(A - 1) + 1
C
                        KCAKI = IT2VO(ISYMCA,ISYMKI)
     &                        + NMATAB(ISYMCA)*(KI - 1) + IOFCA
                        KCAIK = IT2VO(ISYMCA,ISYMKI)
     &                        + NMATAB(ISYMCA)*(IK - 1) + IOFCA
C
                        KOFFT = IOFF5(ISYMA)
     &                        + NT1AM(ISYMCK)*(AI - 1) + IOFCK
C
                        CALL DCOPY(NC,T2VO(KCAKI),1,T1(KOFFT),1)
                        CALL DCOPY(NC,T2VO(KCAIK),1,T2(KOFFT),1)
C
                     ENDDO
C
                  ENDIF
C
               ENDDO
C
            ENDDO
C
         ENDDO
C
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     DO ISYMA = 1,NSYM
c        ISYMCK = MULD2H(ISYMA,ISYMI)
c        DO I = 1,NUMII
c           DO A = 1,NVIR(ISYMA)
c              NAI = NVIR(ISYMA)*(I - 1) + A
c              DO ISYMK = 1,NSYM
c                 ISYMAK = MULD2H(ISYMK,ISYMA)
c                 ISYMC  = MULD2H(ISYMK,ISYMCK)
c                 DO K = 1,NRHF(ISYMK)
c                    NAK = IT1AM(ISYMA,ISYMK) + NVIR(ISYMA)*(K - 1) + A
c                    DO C = 1,NVIR(ISYMC)
c                       NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1)
c    &                      + C
c                       NCI = NVIR(ISYMC)*(I - 1) + C
c                       KCKAI = IOFF5(ISYMA) + NT1AM(ISYMCK)*(NAI - 1)
c    &                        + NCK
c                       KAKCI = IOFF5(ISYMC) + NT1AM(ISYMAK)*(NCI - 1)
c    &                        + NAK
c                       TEST  = T1(KCKAI) - T2(KAKCI)
c                       NTST  = NTST + 1
c                       IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*) '      GXT12: ISYMC,ISYMK,ISYMA,ISYMI: ',
c    &                         ISYMC,ISYMK,ISYMA,ISYMI
c     write(LUPRI,*) '      GXT12: C,K,A,I: ',C,K,A,I
c     write(LUPRI,*) '      GXT12: Diff.  : ',TEST
c     ICOUNE = ICOUNE + 1
c                       ENDIF
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     ndim   = NCKATR(ISYMI)*NUMII
c     xt1nrm = dnorm2(ndim,T1,1)
c     xt2nrm = dnorm2(ndim,T2,1)
c     write(LUPRI,*) '      GXT12: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested (dim = ',ndim,')'
c     write(LUPRI,*) '      GXT12: II1,NUMII,ISYMI   : ',II1,NUMII,ISYMI
c     write(LUPRI,*) '      GXT12: Norm of T1(ck,a#i): ',xt1nrm
c     write(LUPRI,*) '      GXT12: Norm of T2(ck,a#i): ',xt2nrm
C
      RETURN
      END
C  /* Deck ccho_gxt4 */
      SUBROUTINE CCHO_GXT4(T2VO,T4,II1,NUMII,ISYMI,
     &                     IB1,NUMIB,ISYMB,IOFF6)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract amplitudes:
C
C        T4(ck,#b#i)  = T(ck,bi), ordering: IOFF6.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2VO(*),T4(*)
      INTEGER IOFF6(8)
#include <ccorb.h>
#include <ccsdsym.h>
C
C------------------
C     Set up IOFF6.
C------------------
C
      ISYMKC = MULD2H(ISYMI,ISYMB)
      LENBI  = NUMIB*NUMII
      ICOUNT = 0
C
      DO ISYMC = 1,NSYM
C
         ISYMK = MULD2H(ISYMC,ISYMKC)
C
         IOFF6(ISYMC) = ICOUNT
C
         ICOUNT = ICOUNT + LENBI*NRHF(ISYMK)*NVIR(ISYMC)
C
      ENDDO
C
C-------------
C     Extract.
C-------------
C
      DO ISYMIK = 1,NSYM
C
         ISYMBC = ISYMIK
         ISYMK  = MULD2H(ISYMIK,ISYMI)
         ISYMC  = MULD2H(ISYMBC,ISYMB)         
C
         DO K = 1,NRHF(ISYMK)
            DO II = 1,NUMII
C
               IOFBI = NUMIB*(II - 1) + 1
C
               I  = II1 + II - 1
               IK = IMATIJ(ISYMI,ISYMK) + NRHF(ISYMI)*(K - 1) + I
C
               DO C = 1,NVIR(ISYMC)
C
                  IOFBC = IMATAB(ISYMB,ISYMC) + NVIR(ISYMB)*(C - 1)
     &                  + IB1
                  KC    = NRHF(ISYMK)*(C - 1) + K
C
                  KBCIK = IT2VO(ISYMBC,ISYMIK)
     &                  + NMATAB(ISYMBC)*(IK - 1) + IOFBC
C
                  KOFFT = IOFF6(ISYMC) + LENBI*(KC - 1) + IOFBI
C
                  CALL DCOPY(NUMIB,T2VO(KBCIK),1,T4(KOFFT),1)
C
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO

C DEBUG SECTION:
C***************
c     ndim   = LENBI*NT1AM(ISYMKC)
c     xt4nrm = dnorm2(ndim,T4,1)
c     write(LUPRI,*) '      GXT4: IB1,NUMIB,ISYMB      : ',IB1,NUMIB,ISYMB
c     write(LUPRI,*) '      GXT4: II1,NUMII,ISYMI      : ',II1,NUMII,ISYMI
c     write(LUPRI,*) '      GXT4: Norm of T4(#b,#i,k,c): ',xt4nrm
C
      RETURN
      END
C  /* Deck ccho_gxsj12i1 */
      SUBROUTINE CCHO_GXSJ12I1(XJINT,XIINT,XJ1,XJ2,XI1,CHOV,CHOL,
     &                         IB1,NUMIB,ISYMB,IOFF1,IOFF2,IOFT1,NTOT1)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract and scale:
C     J1(aj,#bd) = d(dj,#b) * J(dj,a;#b), ordering: IOFF1
C     I1(aj,#bd) = d(dj,#b) * I(dj,a;#b), ordering: IOFF1
C     J2(da,#bj) = d(dj,a)  * J(dj,a;#b), ordering: IOFF2
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XJINT(*),XIINT(*),XJ1(*),XJ2(*),XI1(*),CHOV(*),CHOL(*)
      INTEGER IOFF1(8),IOFF2(8),IOFT1(8)
      INTEGER AJ,BJ,DJ,DA,BD,DJB,DJA,AJBD,DABJ,DJAB
#include <ccorb.h>
#include <ccsdsym.h>
C
C----------------------------
C     Set up IOFF1 and IOFF2.
C----------------------------
C
      ICOUN1 = 0
      ICOUN2 = 0
C
      DO ISYM = 1,NSYM
C
         ISYMAJ = MULD2H(ISYM,ISYMB)
         ISYMDA = ISYMAJ
C
         IOFF1(ISYM) = ICOUN1
         IOFF2(ISYM) = ICOUN2
C
         ICOUN1 = ICOUN1 + NT1AM(ISYMAJ)*NUMIB*NVIR(ISYM)
         ICOUN2 = ICOUN2 + NMATAB(ISYMDA)*NUMIB*NRHF(ISYM)
C
      ENDDO

c     write(LUPRI,*) '      GXSJ12I1: IB1,NUMIB,ISYMB: ',IB1,NUMIB,ISYMB
c     write(LUPRI,*) '      GXSJ12I1: IOFF1: ',(IOFF1(ISYM),ISYM=1,NSYM)
c     write(LUPRI,*) '      GXSJ12I1: IOFF2: ',(IOFF2(ISYM),ISYM=1,NSYM)
c     write(LUPRI,*) '      GXSJ12I1: IOFT1: ',(IOFT1(ISYM),ISYM=1,NSYM)
C
C-----------------------
C     Extract and scale.
C-----------------------
C
      DO B = 1,NUMIB
         DO ISYMA = 1,NSYM
C
            ISYMDJ = MULD2H(ISYMA,ISYMB)
C
            DO A = 1,NVIR(ISYMA)
               DO ISYMJ = 1,NSYM
C
                  ISYMD  = MULD2H(ISYMJ,ISYMDJ)
                  ISYMAJ = MULD2H(ISYMJ,ISYMA)
                  ISYMDA = MULD2H(ISYMD,ISYMA)
C
                  DO J = 1,NRHF(ISYMJ)
C
                     AJ = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J - 1) + A
                     BJ = NUMIB*(J - 1) + B
C
                     DO D = 1,NVIR(ISYMD)
C
                        DJ = IT1AM(ISYMD,ISYMJ) + NVIR(ISYMD)*(J - 1)
     &                     + D
                        DA = IMATAB(ISYMD,ISYMA) + NVIR(ISYMD)*(A - 1)
     &                     + D
                        BD = NUMIB*(D - 1) + B
C
                        DJB = NTOT1*(B - 1) + IOFT1(ISYMDJ) + DJ
                        DJA = ICKATR(ISYMDJ,ISYMA)
     &                      + NT1AM(ISYMDJ)*(A - 1) + DJ
C
                        AJBD = IOFF1(ISYMD) + NT1AM(ISYMAJ)*(BD - 1)
     &                       + AJ
                        DABJ = IOFF2(ISYMJ) + NMATAB(ISYMDA)*(BJ - 1)
     &                       + DA
C
                        DJAB = NCKATR(ISYMB)*(B - 1) + DJA
C
                        XJ1(AJBD) = CHOV(DJB)*XJINT(DJAB)
                        XI1(AJBD) = CHOV(DJB)*XIINT(DJAB)
                        XJ2(DABJ) = CHOL(DJA)*XJINT(DJAB)

c     write(LUPRI,*)
c     write(LUPRI,*) '      GXSJ12I1: AJBD,DABJ,DJAB: ',AJBD,DABJ,DJAB
c     write(LUPRI,*) '      GXSJ12I1: DJB,DJA       : ',DJB,DJA
c     write(LUPRI,*) '      GXSJ12I1: CHOV,CHOL     : ',CHOV(DJB),CHOL(DJA)
c     write(LUPRI,*) '      GXSJ12I1: XJ1,XI1,XJ2   : ',XJ1(AJBD),
c    &                                              XI1(AJBD),XJ2(DABJ)
c     write(LUPRI,*)
C
                     ENDDO
C
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     ISYMD  = ISYMB
c     NUMID  = NUMIB
c     ID1    = IB1
c     ISYMAJ = 1
c     DO LD = 1,NUMID
c        D = ID1 + LD - 1
c        DO LB = 1,NUMIB
c           B = IB1 + LB - 1
c           NBD = NUMIB*(D - 1) + LB
c           NDB = NUMID*(B - 1) + LD
c           DO ISYMJ = 1,NSYM
c              ISYMA = ISYMJ
c              DO J = 1,NRHF(ISYMJ)
c                 DO A = 1,NVIR(ISYMA)
c                    NAJ = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J - 1) + A
c                    KAJBD = IOFF1(ISYMD)
c    &                     + NT1AM(ISYMAJ)*(NBD - 1)
c    &                     + NAJ
c                    KAJDB = IOFF1(ISYMB)
c    &                     + NT1AM(ISYMAJ)*(NDB - 1)
c    &                     + NAJ
c                    TEST  = XJ1(KAJBD) - XI1(KAJDB)
c                    NTST  = NTST + 1
c                    IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*) '      GXSJ12I1: ISYMA,ISYMJ,ISYMB,ISYMD: ',
c    &                            ISYMA,ISYMJ,ISYMB,ISYMD
c     write(LUPRI,*) '      GXSJ12I1: A,J,B,D: ',A,J,B,D
c     write(LUPRI,*) '      GXSJ12I1: Diff.  : ',TEST
c                       ICOUNE = ICOUNE + 1
c                    ENDIF
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     ndim = NCKATR(ISYMB)*NUMIB
c     WRITE(LUPRI,*) '      GXSJ12I1: ',ICOUNE,' interchange errors ',
c    &           '(',NTST,' tested, dim = ',ndim,')'
c     xj1nrm = dnorm2(ndim,XJ1,1)
c     xj2nrm = dnorm2(ndim,XJ2,1)
c     xi1nrm = dnorm2(ndim,XI1,1)
c     write(LUPRI,*) '      GXSJ12I1: Norm of scaled J1(aj,#bd): ',xj1nrm
c     write(LUPRI,*) '      GXSJ12I1: Norm of scaled I1(aj,#bd): ',xi1nrm
c     write(LUPRI,*) '      GXSJ12I1: Norm of scaled J2(da,#bj): ',xj2nrm
C
      RETURN
      END
C  /* Deck ccho_gst1234 */
      SUBROUTINE CCHO_GST1234(T1,T2,T3,T4,CHOO,II1,NUMII,ISYMI,
     &                        IB1,NUMIB,ISYMB,
     &                        IOFF5,IOFF6,IOFT1,NTOT1)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Scale amplitudes with update vectors.
C
C        T1(ck,a#i)   <- d(ck,#i) * T1(ck,a#i),   ordering: IOFF5.
C        T2(ck,a#i)   <- d(ck,#i) * T2(ck,a#i),   ordering: IOFF5.
C        T3(ck,a#i)   <- d(ck,#i) * T3(ck,a#i),   ordering: IOFF5.
C        T4(#b#i,k,c) <- d(ck,#i) * T4(#b#i,k,c), ordering: IOFF6.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T1(*),T2(*),T3(*),T4(*),CHOO(*)
      INTEGER IOFF5(8),IOFF6(8),IOFT1(8)
      INTEGER AI,CK,BI
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYMBI = MULD2H(ISYMI,ISYMB)
      LENBI  = NUMIB*NUMII
C
      DO ISYM = 1,NSYM
C
C-----------------------
C        T1, T2, and T3.
C-----------------------
C
         ISYMA  = ISYM
         ISYMCK = MULD2H(ISYMA,ISYMI)
C
         DO I = 1,NUMII
            DO A = 1,NVIR(ISYMA)
C
               AI = NVIR(ISYMA)*(I - 1) + A
C
               DO CK = 1,NT1AM(ISYMCK)
C
                  KOFFT = IOFF5(ISYMA)  + NT1AM(ISYMCK)*(AI - 1) + CK
                  KOFFC = NTOT1*(I - 1) + IOFT1(ISYMCK) + CK
C
                  T1(KOFFT) = CHOO(KOFFC)*T1(KOFFT)
                  T2(KOFFT) = CHOO(KOFFC)*T2(KOFFT)
                  T3(KOFFT) = CHOO(KOFFC)*T3(KOFFT)
C
               ENDDO
C
            ENDDO
         ENDDO
C
C-----------
C        T4.
C-----------
C
         ISYMC  = ISYM
         ISYMCK = ISYMBI
         ISYMK  = MULD2H(ISYMC,ISYMCK)
C
         DO C = 1,NVIR(ISYMC)
            DO K = 1,NRHF(ISYMK)
C
               KC = NRHF(ISYMK)*(C - 1) + K
               CK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
               DO I = 1,NUMII
C
                  KOFFC = NTOT1*(I - 1) + IOFT1(ISYMCK) + CK
C
                  DO B = 1,NUMIB
C
                     BI = NUMIB*(I - 1) + B
C
                     KOFFT = IOFF6(ISYMC) + LENBI*(KC - 1) + BI
C
                     T4(KOFFT) = CHOO(KOFFC)*T4(KOFFT)
C
                  ENDDO
C
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO

C DEBUG SECTION
C**************
c     ICOUNE = 0
c     NTST   = 0
c     ISYMCK = MULD2H(ISYMI,ISYMB)
c     DO I = 1,NUMII
c        DO LB = 1,NUMIB
c           B   = IB1 + LB - 1
c           LBI = NUMIB*(I - 1) + LB
c           NBI = NVIR(ISYMB)*(I - 1) + B
c           DO ISYMK = 1,NSYM
c              ISYMC = MULD2H(ISYMK,ISYMCK)
c              DO K = 1,NRHF(ISYMK)
c                 DO C = 1,NVIR(ISYMC)
c                    NCK = IT1AM(ISYMC,ISYMK)
c    &                   + NVIR(ISYMC)*(K - 1) + C
c                    NKC = NRHF(ISYMK)*(C - 1) + K
c                    KCKBI = IOFF5(ISYMB)
c    &                     + NT1AM(ISYMCK)*(NBI - 1) + NCK
c                    KBIKC = IOFF6(ISYMC)
c    &                     + NUMIB*NUMII*(NKC - 1) + LBI
c                    TEST  = T1(KCKBI) - T4(KBIKC)
c                    NTST  = NTST + 1
c                    IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*) '      GST1234: ISYMC,ISYMK,ISYMB,ISYMI: ',
c    &                           ISYMC,ISYMK,ISYMB,ISYMI
c     write(LUPRI,*) '      GST1234: C,K,B,I: ',C,K,B,II1+I-1
c     write(LUPRI,*) '      GST1234: Diff.  : ',TEST
c     ICOUNE = ICOUNE + 1
c                    ENDIF
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     lent1  = NCKATR(ISYMI)*NUMII
c     lent2  = lent1
c     lent3  = lent1
c     lent4  = LENBI*NT1AM(ISYMBI) 
c     write(LUPRI,*) '      GST1234: ',ICOUNE,' T1/T4 errors out of ',
c    &           NTST,' tested (dim T4 = ',lent4,')'
c     t1norm = dnorm2(lent1,T1,1)
c     t2norm = dnorm2(lent2,T2,1)
c     t3norm = dnorm2(lent3,T3,1)
c     t4norm = dnorm2(lent4,T4,1)
c     write(LUPRI,*) '      GST1234: Norm of scaled T1(ck,a#i)  : ',t1norm
c     write(LUPRI,*) '      GST1234: Norm of scaled T2(ck,a#i)  : ',t2norm
c     write(LUPRI,*) '      GST1234: Norm of scaled T3(ck,a#i)  : ',t3norm
c     write(LUPRI,*) '      GST1234: Norm of scaled T4(#b#i,k,c): ',t4norm
C
      RETURN
      END
C  /* Deck ccho_g1term */
      SUBROUTINE CCHO_G1TERM(XIINT,XJINT,T2VO,T1AM,XIAJB,FOCKD,NUMCHO,
     &                       CHOELE,WORK,LWORK,EG1,IB1,ISYMB,NUMIB,
     &                       FBATCH,PRINT,NONA)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Calculate the G1 term:
C
C        EG1 = EG1 + Sum(ba) V(ba) * P(ba)
C
C     where
C
C        V(ba) =   Sum(dj)  d(djb) * L(djba) * t1(dj)
C
C        P(ba) = - Sum(cki) d(cki) * L(ibkc) * t(ai,ck)
C
C     and d(djb) and d(cki) denote
C     the virtual and occupied parts of the Cholesky decomposition of the
C     orbital energy denominator, L(iajb) = 2 * (ia|jb) - (ja|ib), and
C     L(djba) = 2 * (dj|ba) - (bj|da).
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XIINT(*),XJINT(*),T2VO(*),T1AM(*),XIAJB(*),FOCKD(*)
      DIMENSION CHOELE(*)
      DIMENSION WORK(LWORK)
      LOGICAL   FBATCH,PRINT
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.00D0, ZERO = 0.00D0)
      PARAMETER (ONE   =  1.00D0, TWO  = 2.00D0)
C 
      TLAST = SECOND()
C 
C
      IF (PRINT) TIMT = SECOND()

c     lenj  = NCKATR(ISYMB)*NUMIB
c     leni  = lenj
c     liajb = NT2AM(1)
c     lent2 = NT2SQ(1)
c     lent1 = NT1AM(1)
c     xjnorm = dnorm2(lenj,XJINT,1)
c     xinorm = dnorm2(leni,XIINT,1)
c     xniajb = dnorm2(liajb,XIAJB,1)
c     t2norm = dnorm2(lent2,T2VO,1)
c     t1norm = dnorm2(lent1,T1AM,1)
c     if (lwork .ge. lent2) then
c        call cc_t2sq(XIAJB,WORK,1)
c        xsiajb = dnorm2(lent2,WORK,1)
c     else
c        write(LUPRI,*) '   G1TERM: insuf. core for squaring L(iajb)'
c        write(LUPRI,*) '   G1TERM: LWORK = ',LWORK
c        xsiajb = 1.0D9
c     endif
c     if (lwork .ge. lenj) then
c        call dzero(WORK,lenj)
c        call daxpy(lenj,2.0D0,XJINT,1,WORK,1)
c        call daxpy(lenj,-1.0D0,XIINT,1,WORK,1)
c        xn2jmi = dnorm2(lenj,WORK,1)
c     else
c        write(LUPRI,*) '   G1TERM: insuf. core for 2J-I'
c        write(LUPRI,*) '   G1TERM: LWORK = ',LWORK
c        xn2jmi = 1.0D9
c     endif
c     write(LUPRI,*)
c     write(LUPRI,*) '   Entering G1TERM:'
c     write(LUPRI,*) '   ================'
c     write(LUPRI,*) '   LWORK = ',LWORK
c     write(LUPRI,*) '   IB1,NUMIB,ISYMB: ',IB1,NUMIB,ISYMB
c     write(LUPRI,*) '   Norm of J(dj,a;#b): ',xjnorm
c     write(LUPRI,*) '   Norm of I(dj,a;#b): ',xinorm
c     write(LUPRI,*) '   Norm of K=2*J - I : ',xn2jmi
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xniajb,' (packed)'
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xsiajb,' (squared)'
c     write(LUPRI,*) '   Norm of T1AM      : ',t1norm
c     write(LUPRI,*) '   Norm of T2VO      : ',t2norm
c     write(LUPRI,*)
C
C-----------------------
C     Assign symmetries.
C-----------------------
C
      ISYMA  = ISYMB
      ISYMDJ = 1 
      ISYIKC = ISYMB
C
      LENDJ  = NT1AM(ISYMDJ)
      LENDJB = LENDJ*NUMIB
      LENIKC = NMAIJA(ISYIKC)
C
      IF ((LENDJB.LE.0) .OR. (LENIKC.LE.0) .OR. (NUMIB.LE.0)
     &    .OR. (NVIR(ISYMA).LE.0)) RETURN
C
      IF (PRINT) THEN
         WRITE(LUPRI,'(6X,A,/,6X,A,/)')
     &   'Calculation of the G1 term:',
     &   '==========================='
      ENDIF
C
C----------------
C     Allocation.
C----------------
C
      KCHOO = 1
      KCHOV = KCHOO + LENIKC
      KMINT = KCHOV + LENDJB
      KEND1 = KMINT + NUMIB*LENIKC
      LWRK1 = LWORK - KEND1 + 1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Insufficient memory in CCHO_G1TERM'
         WRITE(LUPRI,*) 'Need (more than): ',KEND1-1
         WRITE(LUPRI,*) 'Available       : ',LWORK
         CALL QUIT(' Insufficient memory in CCHO_G1TERM ')
      ENDIF
C
C--------------------------
C     Extract integrals:
C     M(#b,ikc) = L(ib,kc).
C--------------------------
C
      CALL CCHO_G1XM(XIAJB,WORK(KMINT),IB1,NUMIB,ISYMB)

c     lmint = NUMIB*LENIKC
c     xmnrm = dnorm2(lmint,WORK(KMINT),1)
c     write(LUPRI,*) '   G1TERM: Norm of M(#b,ikc): ',xmnrm
c     write(LUPRI,*) '   G1TERM: memory available for a-batch: ',LWRK1
c     write(LUPRI,*) '   G1TERM: LENDJB,LENIKC,NUMIB: ',LENDJB,LENIKC,NUMIB
c     call flshfo(LUPRI)
C
C-------------------------
C     Set up batch over a.
C-------------------------
C
      MINMEM = LENDJB + LENIKC + 2*NUMIB
      IF (FBATCH) THEN
         NEFA  = MIN(NONA,NVIR(ISYMA))
         LEFF  = NEFA*MINMEM + 1
         LWRK1 = MIN(LWRK1,LEFF)
      ENDIF
      NUMA   = MIN(LWRK1/MINMEM,NVIR(ISYMA))
C
      IF (NUMA .LE. 0) THEN
         WRITE(LUPRI,*) 'Insufficient memory for batch in CCHO_G1TERM'
         WRITE(LUPRI,*) 'Memory needed for batch   : ',MINMEM
         WRITE(LUPRI,*) 'Memory available for batch: ',LWRK1
         WRITE(LUPRI,*) 'Memory available in total : ',LWORK
         WRITE(LUPRI,*) 'NUMA, NVIR(ISYMA), ISYMA  : ',NUMA,NVIR(ISYMA),
     &              ISYMA
         CALL QUIT(' Insufficient memory in CCHO_G1TERM ')
      ENDIF
C
      NBATCH = (NVIR(ISYMA) - 1)/NUMA + 1
C
      IF (PRINT) THEN
         WRITE(LUPRI,'(6X,A,I1,A,/,6X,A)')
     &   'Batch over A, symmetry ',ISYMA,':',
     &   '-------------------------'
         WRITE(LUPRI,'(6X,A,I10,/,6X,A,I10)')
     &   'Minimum work space required   : ',MINMEM,
     &   'Work space available for batch: ',LWRK1
         WRITE(LUPRI,'(6X,A,I10,/,6X,A,I10,/)')
     &   'Number of virtual orbitals    : ',NVIR(ISYMA),
     &   'Required number of A-batches  : ',NBATCH
      ENDIF

c     write(LUPRI,*) '   G1TERM: MINMEM,NUMA,NBATCH: ',
c    &               MINMEM,NUMA,NBATCH
c     call flshfo(LUPRI)
C
      DO IBATCH = 1,NBATCH
C
         NUMIA = NUMA
         IF (IBATCH .EQ. NBATCH) THEN
            NUMIA = NVIR(ISYMA) - NUMA*(NBATCH - 1)
         ENDIF
C
         IA1 = NUMA*(IBATCH - 1) + 1
C
         IF (PRINT) THEN
            WRITE(LUPRI,'(9X,A,I10,A,/,9X,A)')
     &      'A-batch number ',IBATCH,':',
     &      '--------------------------'
            WRITE(LUPRI,'(9X,A,I10,1X,I10,/)')
     &      'First and last A: ',IA1,IA1+NUMIA-1
         ENDIF

c     write(LUPRI,*) '   GTERM1: a-batch: IA1,NUMIA,ISYMA: ',
c    &                                IA1,NUMIA,ISYMA
c     write(LUPRI,*) '   GTERM1: a-batch: IBATCH,NBATCH  : ',
c    &                                IBATCH,NBATCH
C
C-------------------
C        Allocation.
C-------------------
C
         LENBA = NUMIB*NUMIA
C
         KKINT = KEND1
         KTG   = KKINT + LENDJB*NUMIA
         KVMAT = KTG   + LENIKC*NUMIA
         KPMAT = KVMAT + LENBA
         KEND2 = KPMAT + LENBA
         LWRK2 = LWORK - KEND2 + 1
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Batching bug in CCHO_G1TERM'
            CALL QUIT(' Error in CCHO_G1TERM ')
         ENDIF
C
C----------------------------------------------
C        Extract integrals:
C        K(dj,#b,#a) = 2*J(dj,a;b) - I(dj,a;b).
C----------------------------------------------
C
         CALL DZERO(WORK(KKINT),LENDJB*NUMIA)
C
         DO B = 1,NUMIB
            DO A = 1,NUMIA
C
               IA  = IA1 + A - 1
               NBA = NUMIB*(A - 1) + B
C
               KOFFJ = NCKATR(ISYMB)*(B - 1) + ICKATR(ISYMDJ,ISYMA)
     &               + NT1AM(ISYMDJ)*(IA - 1) + 1
               KOFFI = KOFFJ
               KOFFK = KKINT + NT1AM(ISYMDJ)*(NBA - 1)
C
               CALL DAXPY(NT1AM(ISYMDJ),TWO,XJINT(KOFFJ),1,
     &                                      WORK(KOFFK),1)
               CALL DAXPY(NT1AM(ISYMDJ),XMONE,XIINT(KOFFI),1,
     &                                        WORK(KOFFK),1)
C
            ENDDO
         ENDDO

c     lenk   = LENDJB*NUMIA
c     xknorm = dnorm2(lenk,WORK(KKINT),1)
c     write(LUPRI,*) '   G1TERM: Norm of K(dj,#b,#a): ',xknorm
C
C--------------------------------------------
C        Extract amplitudes:
C        TG(ikc,#a) = T(ck,ai) = T2VO(ac,ik).
C--------------------------------------------
C
         CALL CCHO_G1XT(T2VO,WORK(KTG),NUMIA,IA1,ISYMA,ISYIKC,
     &                  LENIKC)

c     lentg  = LENIKC*NUMIA
c     tgnorm = dnorm2(lentg,WORK(KTG),1)
c     write(LUPRI,*) '   G1TERM: Norm of TG(ikc,#a): ',tgnorm
C
C----------------------------
C        Start Cholesky loop.
C----------------------------
C
         DO ICHO = 1,NUMCHO
C
            G1COR = ZERO
C
C---------------------------------------
C           Get Cholesky update vectors:
C           d(ikc) and d(dj,#b).
C---------------------------------------
C
            CALL CCHO_DECHO6(FOCKD,CHOELE,NUMCHO,ICHO,WORK(KCHOO),
     &                       WORK(KCHOV),ISYIKC,ISYMDJ,IB1,NUMIB,ISYMB)

c     write(LUPRI,*) '   G1TERM: decho6 completed for ICHO = ',ICHO
C
C--------------------------------------------------
C           Scale virtual integrals and amplitudes:
C           K(dj#b,#a) <- d(dj#b)*K(dj#b,#a)
C           TG(ikc,#a) <- d(ikc)*TG(ikc,#a)
C--------------------------------------------------
C
            CALL CCHO_SCVEC(WORK(KKINT),WORK(KCHOV),LENDJB,NUMIA)
            CALL CCHO_SCVEC(WORK(KTG),WORK(KCHOO),LENIKC,NUMIA)

c     write(LUPRI,*) '   G1TERM: scaling completed'
C
C-----------------------------
C           Calculate V and P.
C-----------------------------
C
            CALL DGEMV('T',LENDJ,LENBA,ONE,WORK(KKINT),LENDJ,
     &                 T1AM,1,ZERO,WORK(KVMAT),1)
C
            CALL DGEMM('N','N',NUMIB,NUMIA,LENIKC,
     &                 XMONE,WORK(KMINT),NUMIB,WORK(KTG),LENIKC,
     &                 ZERO,WORK(KPMAT),NUMIB)
C
C-----------------------------------------
C           Calculate energy contribution.
C-----------------------------------------
C
            G1COR = G1COR + DDOT(LENBA,WORK(KVMAT),1,WORK(KPMAT),1)
C
            EG1 = EG1 + G1COR
C 
              ENERGG(ICHO) = ENERGG(ICHO) + G1COR
C 
C Decommented by Domenico
           TNOW = SECOND()
           DELTAT = TNOW - TLAST
           TLAST = TNOW
           SCNDSG(ICHO) = SCNDSG(ICHO) 
     &                   + DELTAT
C Decommented by Domenico
C 
C
c           IF (PRINT) THEN
c              WRITE(LUPRI,'(12X,A,I3,A,/,12X,A)')
c    &         'Status after Cholesky vector',ICHO,':',
c    &         '--------------------------------'
c              IF (ABS(G1COR) .LT. THRCHO) THEN
c                 WRITE(LUPRI,'(12X,A)') 'G1 term converged'
c              ELSE
c                 WRITE(LUPRI,'(12X,A)') 'G1 term not converged'
c              ENDIF
c              TIM = SECOND() - TIMT
c              WRITE(LUPRI,'(12X,A,F10.2,A,/)')
c    &         'Accumulated G1-time: ',TIM,' seconds'
c           ENDIF
C
            IF (DABS(G1COR) .LT. THRCHO) GOTO 1000
C
         ENDDO
C
 1000    CONTINUE
C
      ENDDO

c     lenj  = NCKATR(ISYMB)*NUMIB
c     leni  = lenj
c     liajb = NT2AM(1)
c     lent2 = NT2SQ(1)
c     lent1 = NT1AM(1)
c     xjnorm = dnorm2(lenj,XJINT,1)
c     xinorm = dnorm2(leni,XIINT,1)
c     xniajb = dnorm2(liajb,XIAJB,1)
c     t2norm = dnorm2(lent2,T2VO,1)
c     t1norm = dnorm2(lent1,T1AM,1)
c     if (lwork .ge. lent2) then
c        call cc_t2sq(XIAJB,WORK,1)
c        xsiajb = dnorm2(lent2,WORK,1)
c     else
c        write(LUPRI,*) '   G1TERM: insuf. core for squaring L(iajb)'
c        write(LUPRI,*) '   LWORK = ',LWORK
c        xsiajb = 1.0D9
c     endif
c     if (lwork .ge. lenj) then
c        call dzero(WORK,lenj)
c        call daxpy(lenj,2.0D0,XJINT,1,WORK,1)
c        call daxpy(lenj,-1.0D0,XIINT,1,WORK,1)
c        xn2jmi = dnorm2(lenj,WORK,1)
c     else
c        write(LUPRI,*) '   G1TERM: insuf. core for 2J-I'
c        write(LUPRI,*) '   G1TERM: LWORK = ',LWORK
c        xn2jmi = 1.0D9
c     endif
c     write(LUPRI,*)
c     write(LUPRI,*) '   Exiting G1TERM:'
c     write(LUPRI,*) '   ==============='
c     write(LUPRI,*) '   LWORK = ',LWORK
c     write(LUPRI,*) '   IB1,NUMIB,ISYMB: ',IB1,NUMIB,ISYMB
c     write(LUPRI,*) '   Norm of J(dj,a;#b): ',xjnorm
c     write(LUPRI,*) '   Norm of I(dj,a;#b): ',xinorm
c     write(LUPRI,*) '   Norm of K=2*J - I : ',xn2jmi
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xniajb,' (packed)'
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xsiajb,' (squared)'
c     write(LUPRI,*) '   Norm of T1AM      : ',t1norm
c     write(LUPRI,*) '   Norm of T2VO      : ',t2norm
c     write(LUPRI,*)
C
      RETURN
      END
C  /* Deck ccho_g1xm */
      SUBROUTINE CCHO_G1XM(XIAJB,XMINT,IB1,NUMIB,ISYMB)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract integrals:
C
C        M(#b,ikc) = L(ib,kc)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XIAJB(*)
      DIMENSION XMINT(NUMIB,*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      DO ISYMCK = 1,NSYM
C
         ISYMBI = ISYMCK
         ISYMI  = MULD2H(ISYMBI,ISYMB)
C
         IF (NRHF(ISYMI) .GT. 0) THEN
C
            DO ISYMK = 1,NSYM
C
               ISYMC  = MULD2H(ISYMK,ISYMCK)
               ISYMIK = MULD2H(ISYMK,ISYMI)
C
               DO K = 1,NRHF(ISYMK)
                  DO C = 1,NVIR(ISYMC)
C
                     NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
                     DO I = 1,NRHF(ISYMI)
C
                        IK  = IMATIJ(ISYMI,ISYMK) + NRHF(ISYMI)*(K - 1)
     &                      + I
                        IKC = IMAIJA(ISYMIK,ISYMC)
     &                      + NMATIJ(ISYMIK)*(C - 1) + IK
C
                        DO B = 1,NUMIB
C
                           NBI = IT1AM(ISYMB,ISYMI)
     &                         + NVIR(ISYMB)*(I - 1) + IB1 + B - 1
C
                           NBICK = IT2AM(ISYMBI,ISYMCK) + INDEX(NBI,NCK)
C
                           XMINT(B,IKC) = XIAJB(NBICK)
C
                        ENDDO
C
                     ENDDO
C
                  ENDDO
               ENDDO
C
            ENDDO
C
         ENDIF
C
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_g1xt */
      SUBROUTINE CCHO_G1XT(T2VO,T2IJCA,NA,IA1,ISYMA,ISYIJC,NIJC)
C
C     Henrik Koch and Alfredo Sanchez.    Aug 1999
C
C     Extract (same as CCHO_TIJCA, except for *not* extracting 2CME):
C
C     T2IJCA(ijc,a) = T2VO(ac,ij)
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION T2VO(*),T2IJCA(NIJC,NA)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DO A = 1,NA
C
         IA  = IA1 + A - 1
C
         DO ISYMC = 1,NSYM
C
            ISYMAC = MULD2H(ISYMA,ISYMC)
            ISYMIJ = MULD2H(ISYMC,ISYIJC)
C
            NIJ = NMATIJ(ISYMIJ)
C
            DO C = 1,NVIR(ISYMC)
C
               NAC = IMATAB(ISYMA,ISYMC) + NVIR(ISYMA)*(C-1) + IA
C
               DO IJ = 1,NIJ
C
                  IJC   = IMAIJA(ISYMIJ,ISYMC)
     &                  + NMATIJ(ISYMIJ)*(C-1) + IJ
C
                  NACIJ = IT2VO(ISYMAC,ISYMIJ)
     &                  + NMATAB(ISYMAC)*(IJ-1) + NAC
C
                  T2IJCA(IJC,A) = T2VO(NACIJ)
C
               ENDDO
C
            ENDDO
C
         ENDDO
C
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_hterm */
      SUBROUTINE CCHO_HTERM(FOCKD,T1AM,XKINT,T2VO,XIAJB,WORK,LWORK,
     &                      CHOELE,NUMCHO,EH,FBATCH,PRINT,
     &                      NONI,NONK)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Calculate the H term:
C
C        EH = EH + Sum(aijk) [ V(ai,jk) * P(ai,kj) + W(ai,jk) * Q(ai,kj) ]
C
C                + Sum(ijkl) X(lk,ji) * R(lk,ji)
C
C                + Sum(ki)   Y(ki) * S(ki)
C
C     where
C
C        V(ai,jk) =   Sum(dl)  d(dlj) * L(ldia) * L(dlkj)
C
C        W(ai,jk) =   Sum(dl)  d(dlj) * L(ldia) * (dl|kj)
C
C                   + Sum(dl)  d(dlj) * L(laid) * (dj|kl)
C
C        X(lk,ji) =   Sum(cd)  d(dlc) * t(ck,dl) * L(icjd)
C
C        Y(ki)    =   Sum(cdl) d(dlc) * t(ck,dl) * L(icld) =  Sum(j) X(jk,ji)
C
C        P(ai,kj) =   Sum(b)   d(aib) * t(ai,bk) * t1(bj)
C
C        Q(ai,kj) = - Sum(ck)  d(aib) * t(ak,bi) * t1(bj)
C
C        R(lk,ji) = - Sum(b)   d(bji) * (bj|ki) * t1(bl)
C
C        S(ki)    =   Sum(em)  d(emi) * L(emki) * t1(em)
C
C     and s(ai,bj) = 2 * t(ai,bj) - t(aj,bi), d(aij) and d(aib) denote
C     the occupied and virtual parts of the Cholesky decomposition of the
C     orbital energy denominator, and L(iajb) = 2 * (ia|jb) - (ja|ib) etc.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION FOCKD(*),T1AM(*)
      DIMENSION XKINT(*),T2VO(*),XIAJB(*)
      DIMENSION WORK(LWORK)
      DIMENSION CHOELE(*)
      LOGICAL   FBATCH,PRINT
      INTEGER IOFM12(8),IOFM3(8)
      INTEGER IOFK12(8),IOFK3(8)
      INTEGER IOFT12(8),IOFT3(8)
      INTEGER BJ
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.00D0, ZERO = 0.00D0)
      PARAMETER (ONE   =  1.00D0, TWO  = 2.00D0)
C 
      TLAST = SECOND()
C 
C
      IF (PRINT) THEN
         TIMT = SECOND()
         WRITE(LUPRI,'(A,/,A,/)')
     &   'Calculation of the H term:',
     &   '=========================='
      ENDIF

c     evp   = 0.0D0
c     ewq   = 0.0D0
c     exr   = 0.0D0
c     eys   = 0.0D0
c     lenk  = NTRAOC(1)
c     liajb = NT2AM(1)
c     lent2 = NT2SQ(1)
c     lent1 = NT1AM(1)
c     xknorm = dnorm2(lenk,XKINT,1)
c     xniajb = dnorm2(liajb,XIAJB,1)
c     t2norm = dnorm2(lent2,T2VO,1)
c     t1norm = dnorm2(lent1,T1AM,1)
c     if (lwork .ge. lent2) then
c        call cc_t2sq(XIAJB,WORK,1)
c        xsiajb = dnorm2(lent2,WORK,1)
c     else
c        write(LUPRI,*) '   HTERM: insuf. core for squaring L(iajb)'
c        write(LUPRI,*) '   HTERM: LWORK = ',LWORK,'. Need = ',lent2
c        xsiajb = 1.0D9
c     endif
c     if (lwork .ge. lenk) then
c        call dzero(WORK,lenk)
c        call daxpy(lenk,2.0D0,XKINT,1,WORK,1)
c        DO ISYME = 1,NSYM
c           ISYIKM = ISYME
c           DO E = 1,NVIR(ISYME)
c              DO ISYMM = 1,NSYM
c                 ISYMIK = MULD2H(ISYMM,ISYIKM)
c                 DO M = 1,NRHF(ISYMM)
c                    DO ISYMK = 1,NSYM
c                       ISYMI  = MULD2H(ISYMK,ISYMIK)
c                       ISYMMK = MULD2H(ISYMK,ISYMM)
c                       DO K = 1,NRHF(ISYMK)
c                          DO I = 1,NRHF(ISYMI)
c                             KIKME = ISJIKA(ISYIKM,ISYME)
c    &                              + NMAJIK(ISYIKM)*(E - 1)
c    &                              + ISJIK(ISYMIK,ISYMM)
c    &                              + NMATIJ(ISYMIK)*(M - 1)
c    &                              + IMATIJ(ISYMI,ISYMK)
c    &                              + NRHF(ISYMI)*(K - 1) + I
c                             KMKIE = ISJIKA(ISYIKM,ISYME)
c    &                              + NMAJIK(ISYIKM)*(E - 1)
c    &                              + ISJIK(ISYMMK,ISYMI)
c    &                              + NMATIJ(ISYMMK)*(I - 1)
c    &                              + IMATIJ(ISYMM,ISYMK)
c    &                              + NRHF(ISYMM)*(K - 1) + M
c                             WORK(KMKIE) = WORK(KMKIE) - XKINT(KIKME)
c                          ENDDO
c                       ENDDO
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c        xk2cme = dnorm2(lenk,WORK,1)
c     else
c        write(LUPRI,*) '   HTERM: insuf. core for K 2CME'
c        write(LUPRI,*) '   HTERM: LWORK = ',LWORK,'. Need = ',lenk
c        xk2cme = 1.0D9
c     endif
c     write(LUPRI,*)
c     write(LUPRI,*) '   Entering HTERM:'
c     write(LUPRI,*) '   ==============='
c     write(LUPRI,*) '   LWORK = ',LWORK
c     write(LUPRI,*) '   Norm of K(lj,k;d) : ',xknorm
c     write(LUPRI,*) '   Norm of K 2CME    : ',xk2cme
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xniajb,' (packed)'
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xsiajb,' (squared)'
c     write(LUPRI,*) '   Norm of T1AM      : ',t1norm
c     write(LUPRI,*) '   Norm of T2VO      : ',t2norm
c     write(LUPRI,*)
C
C----------------------
C     Split work space.
C----------------------
C
      LWRKI = LWORK/2
C
      IF (LWRKI .LE. 0) THEN
         WRITE(LUPRI,*) 'Insufficient memory in CCHO_HTERM'
         WRITE(LUPRI,*) 'Memory available: ',LWORK
         CALL QUIT(' Insufficient memory in CCHO_HTERM ')
      ENDIF
C
C----------------------------
C     Loop over i symmetries.
C----------------------------
C
      DO ISYMI = 1,NSYM
C
         IF (NRHF(ISYMI) .LE. 0) GOTO 1000
C
C---------------------
C        Set up batch.
C---------------------
C
         MINMEM = 3*NCKATR(ISYMI)
         IF (FBATCH) THEN
            NEFI  = MIN(NONI,NRHF(ISYMI))
            LEFF  = NEFI*MINMEM + 1
            LWRKB = MIN(LWRKI,LEFF)
         ELSE
            LWRKB = LWRKI
         ENDIF
         NUMI = MIN(LWRKB/MINMEM,NRHF(ISYMI))
C
         IF (NUMI .LE. 0) THEN
            WRITE(LUPRI,*) 'Insufficient memory for batch in CCHO_HTERM'
            WRITE(LUPRI,*) 'Minimum memory required   : ',MINMEM,
     &                 ' (for ISYMI = ',ISYMI,')'
            WRITE(LUPRI,*) 'Memory available for batch: ',LWRKB
            WRITE(LUPRI,*) 'Memory available in total : ',LWORK
            CALL QUIT(' Insufficient memory in CCHO_HTERM ')
         ENDIF
C
         NBATI = (NRHF(ISYMI) - 1)/NUMI + 1
C
         IF (PRINT) THEN
            WRITE(LUPRI,'(3X,A,I1,A,/,3X,A)')
     &      'Batch over I, symmetry ',ISYMI,':',
     &      '-------------------------'
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10)')
     &      'Minimum work space required   : ',MINMEM,
     &      'Work space available for batch: ',LWRKB
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10,/)')
     &      'Number of occupied orbitals   : ',NRHF(ISYMI),
     &      'Required number of I-batches  : ',NBATI
         ENDIF
C
         DO IBATI = 1,NBATI
C
            NUMII = NUMI
            IF (IBATI .EQ. NBATI) THEN
               NUMII = NRHF(ISYMI) - NUMI*(NBATI - 1)
            ENDIF
C
            II1 = NUMI*(IBATI - 1) + 1
C
            IF (PRINT) THEN
               WRITE(LUPRI,'(6X,A,I10,A,/,6X,A)')
     &         'I-batch number ',IBATI,':',
     &         '--------------------------'
               WRITE(LUPRI,'(6X,A,I10,1X,I10,/)')
     &         'First and last I: ',II1,II1+NUMII-1
            ENDIF
C
C----------------------
C           Allocation.
C----------------------
C
            KM1   = 1
            KM2   = KM1   + NCKATR(ISYMI)*NUMII
            KM3   = KM2   + NCKATR(ISYMI)*NUMII
            KEND1 = KM3   + NCKATR(ISYMI)*NUMII 
            LWRK1 = LWORK - KEND1 + 1
C
            IF (LWRK1 .LT. 0) THEN
               WRITE(LUPRI,*) 'Batching bug in CCHO_HTERM (1)'
               CALL QUIT(' Batching error in CCHO_HTERM ')
            ENDIF
C
C--------------------------------------------------
C           Extract integrals:
C           M1(dl,a#i) = L(ld,ia), ordering: IOFM12
C           M2(dl,a#i) = L(la,id), ordering: IOFM12
C           M3(cd,j#i) = L(jd,ic), ordering: IOFM3
C--------------------------------------------------
C
            CALL CCHO_HXM123(XIAJB,WORK(KM1),WORK(KM2),WORK(KM3),
     &                       II1,NUMII,ISYMI,IOFM12,IOFM3)

c     lenm   = NCKATR(ISYMI)*NUMII
c     xniajb = dnorm2(nt2am(1),XIAJB,1)
c     xm1nrm = dnorm2(lenm,WORK(KM1),1)
c     xm2nrm = dnorm2(lenm,WORK(KM2),1)
c     xm3nrm = dnorm2(lenm,WORK(KM3),1)
c     write(LUPRI,*) '   HTERM: Norm of L(iajb)   : ',xniajb,' (packed)'
c     write(LUPRI,*) '   HTERM: Norm of M1(dl,a#i): ',xm1nrm
c     write(LUPRI,*) '   HTERM: Norm of M2(dl,a#i): ',xm2nrm
c     write(LUPRI,*) '   HTERM: Norm of M3(cd,j#i): ',xm3nrm
c     if ((nsym.eq.1) .and. (numii.eq.nrhft)) then
c        write(LUPRI,*) '   HTERM: Printing L(ia,jb) [packed]:'
c        call outpak(XIAJB,NT1AMX,1,6)
c        write(LUPRI,*) '   HTERM: Printing M1:'
c        call output(WORK(KM1),1,NT1AMX,1,NT1AMX,NT1AMX,NT1AMX,1,6)
c        write(LUPRI,*) '   HTERM: Printing M2:'
c        call output(WORK(KM2),1,NT1AMX,1,NT1AMX,NT1AMX,NT1AMX,1,6)
c     endif
C
C----------------------------------
C           Loop over k symmetries.
C----------------------------------
C
            DO ISYMK = 1,NSYM
C
               IF (NRHF(ISYMK) .LE. 0) GOTO 999
C
               ISYMKI = MULD2H(ISYMK,ISYMI)
C
               MAXAJ = -1
               MAXLJ = -1
               DO ISYMJ = 1,NSYM
                  ISYMA = MULD2H(ISYMJ,ISYMKI)
                  MAXAJ = MAX(MAXAJ,NVIR(ISYMA)*NRHF(ISYMJ))
                  MAXLJ = MAX(MAXLJ,NRHF(ISYMA)*NRHF(ISYMJ))
               ENDDO
C
C----------------------------------------------
C              Allocation for Cholesky vectors.
C----------------------------------------------
C
               KCHOO = KEND1
               KCHOV = KCHOO + NCKI(ISYMK)
               KEND2 = KCHOV + NCKATR(ISYMK)
               LWRK2 = LWORK - KEND2 + 1
C
               IF (LWRK2 .LE. 0) THEN
                  WRITE(LUPRI,*) 'Insufficient memory in CCHO_HTERM'
                  WRITE(LUPRI,*) 'Need (more than): ',KEND2-1
                  WRITE(LUPRI,*) 'Available       : ',LWORK
                  CALL QUIT(' Insufficient memory in CCHO_HTERM ')
               ENDIF
C
C---------------------------
C              Set up batch.
C---------------------------
C
               MINMEM = NCKATR(ISYMK)
     &                + 2*NCKI(ISYMK)
     &                + 2*NMATAB(ISYMKI)*NUMII
     &                + NT1AM(ISYMKI)*NUMII
     &                + 4*MAXAJ*NUMII
     &                + 2*MAXLJ*NUMII
               IF (ISYMK .EQ. ISYMI) THEN
                  MINMEM = MINMEM
     &                   + NT1AM(ISYMKI)*NUMII
     &                   + 2*NUMII
               ENDIF
               IF (FBATCH) THEN
                  NEFK  = MIN(NONK,NRHF(ISYMK))
                  LEFF  = NEFK*MINMEM + 1
                  LWRK2 = MIN(LWRK2,LEFF)
               ENDIF
               NUMK = MIN(LWRK2/MINMEM,NRHF(ISYMK))
C
               IF (NUMK .LE. 0) THEN
                  WRITE(LUPRI,*)
     &            'Insufficient memory for batch in CCHO_HTERM'
                  WRITE(LUPRI,*)
     &            'Minimum memory required   : ',MINMEM,
     &            ' (for ISYMI,ISYMK = ',ISYMI,ISYMK,')'
                  WRITE(LUPRI,*) 'Memory available for batch: ',LWRK2
                  WRITE(LUPRI,*) 'Memory available in total : ',LWORK
                  CALL QUIT(' Insufficient memory in CCHO_HTERM ')
               ENDIF
C
               NBATK = (NRHF(ISYMK) - 1)/NUMK + 1
C
               IF (PRINT) THEN
                  WRITE(LUPRI,'(9X,A,I1,A,/,9X,A)')
     &            'Batch over K, symmetry ',ISYMK,':',
     &            '-------------------------'
                  WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10)')
     &            'Minimum work space required   : ',MINMEM,
     &            'Work space available for batch: ',LWRK2
                  WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10,/)')
     &            'Number of occupied orbitals   : ',NRHF(ISYMK),
     &            'Required number of K-batches  : ',NBATK
               ENDIF
C
               DO IBATK = 1,NBATK
C
                  NUMIK = NUMK
                  IF (IBATK .EQ. NBATK) THEN
                     NUMIK = NRHF(ISYMK) - NUMK*(NBATK - 1)
                  ENDIF
C
                  IK1 = NUMK*(IBATK - 1) + 1
C
                  IF (PRINT) THEN
                     WRITE(LUPRI,'(12X,A,I10,A,/,12X,A)')
     &               'K-batch number ',IBATK,':',
     &               '--------------------------'
                     WRITE(LUPRI,'(12X,A,I10,1X,I10,/)')
     &               'First and last K: ',IK1,IK1+NUMIK-1
                  ENDIF
C
                  LENKI = NUMIK*NUMII
C
C----------------------------
C                 Allocation.
C----------------------------
C
                  KT3   = KEND2
                  KK1   = KT3   + NCKATR(ISYMK)*NUMIK
                  KK2   = KK1   + NCKI(ISYMK)*NUMIK
                  KK3   = KK2   + NCKI(ISYMK)*NUMIK
                  KT1   = KK3   + NT1AM(ISYMKI)*LENKI
                  KT2   = KT1   + NMATAB(ISYMKI)*LENKI
                  KVMAT = KT2   + NMATAB(ISYMKI)*LENKI
                  KWMAT = KVMAT + MAXAJ*LENKI
                  KPMAT = KWMAT + MAXAJ*LENKI
                  KQMAT = KPMAT + MAXAJ*LENKI
                  KXMAT = KQMAT + MAXAJ*LENKI
                  KRMAT = KXMAT + MAXLJ*LENKI
                  KEND3 = KRMAT + MAXLJ*LENKI
                  IF (ISYMK .EQ. ISYMI) THEN
                     KK4   = KEND3
                     KYMAT = KK4   + NT1AM(ISYMKI)*LENKI
                     KSMAT = KYMAT + LENKI
                     KEND3 = KSMAT + LENKI
                  ENDIF
                  LWRK3 = LWORK - KEND3 + 1
C
                  IF (LWRK3 .LT. 0) THEN
                     WRITE(LUPRI,*) 'Batching bug in CCHO_HTERM (2)'
                     CALL QUIT(' Batching error in CCHO_HTERM ')
                  ENDIF
C
C------------------------------------------------------------------
C                 Extract integrals:
C                 K1(dl,j#k) = (dl|kj) = K(lkj,d), ordering: IOFK12
C                 K2(dl,j#k) = (dj|kl) = K(jkl,d), ordering: IOFK12
C------------------------------------------------------------------
C
                  CALL CCHO_HXK12(XKINT,WORK(KK1),WORK(KK2),
     &                            IK1,NUMIK,ISYMK,IOFK12)

c     lenk   = NCKI(ISYMK)*NUMIK
c     xk1nrm = dnorm2(lenk,WORK(KK1),1)
c     xk2nrm = dnorm2(lenk,WORK(KK2),1)
c     write(LUPRI,*) '   HTERM: Norm of K1(dl,j#k),sI,sK : ',
c    &            xk1nrm,ISYMI,ISYMK
c     write(LUPRI,*) '   HTERM: Norm of K2(dl,j#k),sI,sK : ',
c    &            xk2nrm,ISYMI,ISYMK
c     if ((nsym.eq.1) .and. (numik.eq.nrhft)) then
c        write(LUPRI,*) '   HTERM: Printing K1 after extraction: '
c        call output(WORK(KK1),1,NT1AMX,1,NMATIJ(1),NT1AMX,NMATIJ(1),
c    &               1,6)
c        write(LUPRI,*) '   HTERM: Printing K2 after extraction: '
c        call output(WORK(KK2),1,NT1AMX,1,NMATIJ(1),NT1AMX,NMATIJ(1),
c    &               1,6)
c     endif
C
C------------------------------------------------------------------
C                 Extract integrals:
C                 K3(b,#kj#i) = (bj,ki) = K(jki,b), ordering: IOFK3
C------------------------------------------------------------------
C
                  CALL CCHO_HXK3(XKINT,WORK(KK3),IK1,NUMIK,ISYMK,
     &                           II1,NUMII,ISYMI,IOFK3)

c     lenk3  = NT1AM(ISYMKI)*LENKI
c     xk3nrm = dnorm2(lenk3,WORK(KK3),1)
c     write(LUPRI,*) '   HTERM: Norm of K3(b,#kj#i): ',xk3nrm
C
C-----------------------------------------------------
C                 Extract integrals:
C                 K4(em,#k#i) = 2*(em,ki)  - (ei,km)
C                             = 2*K(mki,e) - K(ikm,e).
C-----------------------------------------------------
C
                  IF (ISYMK .EQ. ISYMI) THEN
                     CALL CCHO_HXK4(XKINT,WORK(KK4),IK1,NUMIK,ISYMK,
     &                              II1,NUMII,ISYMI)

c     lenk4  = NT1AM(ISYMKI)*LENKI
c     xk4nrm = dnorm2(lenk4,WORK(KK4),1)
c     write(LUPRI,*) '   HTERM: Norm of K4(b,#kj#i): ',xk4nrm

                  ENDIF
C
C----------------------------------------------------------------------
C                 Extract amplitudes:
C                 T3(cd,l#k) = T(ck,dl) = T2VO(cd,kl), ordering: IOFT3.
C----------------------------------------------------------------------
C
                  CALL CCHO_HXT3(T2VO,WORK(KT3),IK1,NUMIK,ISYMK,IOFT3)

c     lent3  = NCKATR(ISYMK)*NUMIK
c     t3norm = dnorm2(lent3,WORK(KT3),1)
c     write(LUPRI,*) '   HTERM: Norm of T3(cd,l#k) : ',t3norm
C
C-----------------------------------------------------------------------
C                 Extract amplitudes:
C                 T1(a#i,#kb) = T(ai,bk) = T2VO(ab,ik), ordering: IOFT12
C                 T2(a#i,#kb) = T(ak,bi) = T2VO(ab,ki), ordering: IOFT12
C-----------------------------------------------------------------------
C
                  CALL CCHO_HXT12(T2VO,WORK(KT1),WORK(KT2),
     &                            IK1,NUMIK,ISYMK,II1,NUMII,ISYMI,
     &                            IOFT12)

c     lent1  = NMATAB(ISYMKI)*LENKI
c     lent2  = lent1
c     t1norm = dnorm2(lent1,WORK(KT1),1)
c     t2norm = dnorm2(lent2,WORK(KT2),1)
c     write(LUPRI,*) '   HTERM: Norm of T1(a#i,#kb): ',t1norm
c     write(LUPRI,*) '   HTERM: Norm of T2(a#i,#kb): ',t2norm
C
C-------------------------------------
C                 Start Cholesky loop.
C-------------------------------------
C
c     write(LUPRI,*) 'HTERM: ****NOTICE**** fixing ICHO = 1'
c                 do icho = 1,1
                  DO ICHO = 1,NUMCHO
C
C-------------------------------
C                    Initialize.
C-------------------------------
C
                     IF (ISYMK .EQ. ISYMI) THEN
                        CALL DZERO(WORK(KYMAT),LENKI)
                     ENDIF
C
                     HCOR = ZERO
C
C-----------------------------------------------------------
C                    Get Cholesky vectors:
C                    d(dl,j), update vectors, symmetry=ISYMK
C                    d(cd,l), update vectors, symmetry=ISYMK
C-----------------------------------------------------------
C
                     CALL CCHO_DECHO8(FOCKD,CHOELE,NUMCHO,ICHO,
     &                                WORK(KCHOO),WORK(KCHOV),
     &                                ISYMK,ISYMK)

c     write(LUPRI,*) 'HTERM: ****NOTICE**** setting vecs. = 1'
c     do idum = 1,NCKI(ISYMK)
c        koff = KCHOO + idum - 1
c        WORK(koff) = 1.0D0
c     enddo
c     do idum = 1,NCKATR(ISYMK)
c        koff = KCHOV + idum - 1
c        WORK(koff) = 1.0D0
c     enddo

c     lendlj = NCKI(ISYMK)
c     xdlj   = 1.0D0*lendlj
c     lencdl = NCKATR(ISYMK)
c     ocnorm = dnorm2(lendlj,WORK(KCHOO),1)
c     vinorm = dnorm2(lencdl,WORK(KCHOV),1)
c     xcdl   = 1.0D0*lencdl
c     ocnorm = ddot(lendlj,WORK(KCHOO),1,WORK(KCHOO),1)/xdlj
c     vinorm = ddot(lencdl,WORK(KCHOV),1,WORK(KCHOV),1)/xcdl
c     write(LUPRI,*) '   HTERM: N-Norm of occ. update vector ',ICHO,':',
c    &           ocnorm
c     write(LUPRI,*) '   HTERM: N-Norm of vir. update vector ',ICHO,':',
c    &           vinorm
C
C-------------------------------------------------------
C                    Scale integrals:
C                    K1(dl,j#k)  <- d(dl,j) * K1(dl,j#k)
C                    K2(dl,j#k)  <- d(dl,j) * K2(dl,j#k)
C-------------------------------------------------------
C
                     DO ISYMJ = 1,NSYM
C
                        ISYMDL = MULD2H(ISYMJ,ISYMK)
C
                        LENDLJ = NT1AM(ISYMDL)*NRHF(ISYMJ)
C
                        KOFK1 = KK1 + IOFK12(ISYMJ)
                        KOFK2 = KK2 + IOFK12(ISYMJ)
                        KOFFD = KCHOO + ICKI(ISYMDL,ISYMJ)
C
                        CALL CCHO_SCVEC(WORK(KOFK1),WORK(KOFFD),
     &                                  LENDLJ,NUMIK)
                        CALL CCHO_SCVEC(WORK(KOFK2),WORK(KOFFD),
     &                                  LENDLJ,NUMIK)
C
                     ENDDO

c     call chk_hk12(WORK(KK1),WORK(KK2),IK1,NUMIK,ISYMK,IOFK12)
c     lenk   = NCKI(ISYMK)*NUMIK
c     xk1nrm = dnorm2(lenk,WORK(KK1),1)
c     xk2nrm = dnorm2(lenk,WORK(KK2),1)
c     write(LUPRI,*) '   HTERM: After scal.: Norm of K1,sI,sK: ',
c    &           xk1nrm,ISYMI,ISYMK
c     write(LUPRI,*) '   HTERM: After scal.: Norm of K2,sI,sK: ',
c    &           xk2nrm,ISYMI,ISYMK
c     if ((nsym.eq.1) .and. (numik.eq.nrhft)) then
c        write(LUPRI,*) '   HTERM: Printing K1 after scaling: '
c        call output(WORK(KK1),1,NT1AMX,1,NMATIJ(1),NT1AMX,NMATIJ(1),
c    &               1,6)
c        write(LUPRI,*) '   HTERM: Printing K2 after scaling: '
c        call output(WORK(KK2),1,NT1AMX,1,NMATIJ(1),NT1AMX,NMATIJ(1),
c    &               1,6)
c     endif
C
C--------------------------------------------------------
C                    Scale integrals:
C                    K3(b,#kj#i) <- d(bj,i) * K3(b,#kj#i)
C--------------------------------------------------------
C
                     ISYMBJ = ISYMKI
C
                     DO ISYMJ = 1,NSYM
C
                        ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
                        DO I = 1,NUMII
C
                           II = II1 + I - 1
C
                           DO J = 1,NRHF(ISYMJ)
C
                              DO K = 1,NUMIK
C
                                 KJI = NUMIK*NRHF(ISYMJ)*(I - 1)
     &                               + NUMIK*(J - 1) + K
C
                                 DO B = 1,NVIR(ISYMB)
C
                                    KOFK3 = KK3 + IOFK3(ISYMJ)
     &                                    + NVIR(ISYMB)*(KJI - 1)
     &                                    + B - 1
                                    KOFFD = KCHOO + ICKI(ISYMBJ,ISYMI)
     &                                    + NT1AM(ISYMBJ)*(II - 1)
     &                                    + IT1AM(ISYMB,ISYMJ)
     &                                    + NVIR(ISYMB)*(J - 1) + B - 1
C
                                    WORK(KOFK3) = WORK(KOFFD)
     &                                          * WORK(KOFK3)
C
                                 ENDDO
C
                              ENDDO
C
                           ENDDO
C
                        ENDDO
C
                     ENDDO
C
C--------------------------------------------------------
C                    Scale integrals:
C                    K4(bj,#k#i) <- d(bj,i) * K4(bj,#k#i)
C--------------------------------------------------------
C
                     IF (ISYMK .EQ. ISYMI) THEN
C
                        DO I = 1,NUMII
C
                           II = II1 + I - 1
C
                           DO K = 1,NUMIK
C
                              KI = NUMIK*(I - 1) + K
C
                              DO BJ = 1,NT1AM(1)
C
                                 KOFK4 = KK4 + NT1AM(1)*(KI - 1)
     &                                 + BJ - 1
                                 KOFFD = KCHOO + ICKI(1,ISYMI)
     &                                 + NT1AM(1)*(II - 1) + BJ - 1
C
                                 WORK(KOFK4) = WORK(KOFFD)*WORK(KOFK4)
C
                              ENDDO
C
                           ENDDO
C
                        ENDDO
C
                     ENDIF
C
C--------------------------------------------------------
C                    Scale amplitudes:
C                    T1(a#i,#kb) <- d(ab,i) * T1(a#i,#kb)
C                    T2(a#i,#kb) <- d(ab,i) * T2(a#i,#kb)
C--------------------------------------------------------
C
                     ISYMAB = ISYMKI
C
                     DO ISYMB = 1,NSYM
C
                        ISYMA = MULD2H(ISYMB,ISYMAB)
C
                        LENAI  = NVIR(ISYMA)*NUMII
                        LENAIK = LENAI*NUMIK
C
                        DO B = 1,NVIR(ISYMB)
                           DO K = 1,NUMIK
                              DO I = 1,NUMII
C
                                 II = II1 + I - 1
C
                                 DO A = 1,NVIR(ISYMA)
C
                                    KOFT1 = KT1
     &                                    + IOFT12(ISYMB)
     &                                    + LENAIK*(B - 1)
     &                                    + LENAI*(K - 1)
     &                                    + NVIR(ISYMA)*(I - 1)
     &                                    + A - 1
                                    KOFT2 = KT2
     &                                    + IOFT12(ISYMB)
     &                                    + LENAIK*(B - 1)
     &                                    + LENAI*(K - 1)
     &                                    + NVIR(ISYMA)*(I - 1)
     &                                    + A - 1
                                    KOFFD = KCHOV
     &                                    + ICKASR(ISYMAB,ISYMI)
     &                                    + NMATAB(ISYMAB)*(II - 1)
     &                                    + IMATAB(ISYMA,ISYMB)
     &                                    + NVIR(ISYMA)*(B - 1)
     &                                    + A - 1
C
                                    WORK(KOFT1) = WORK(KOFT1)
     &                                          * WORK(KOFFD)
                                    WORK(KOFT2) = WORK(KOFT2)
     &                                          * WORK(KOFFD)
C
                                 ENDDO
C
                              ENDDO
                           ENDDO
                        ENDDO
C
                     ENDDO

c     lent   = NMATAB(ISYMAB)*NUMII*NUMIK
c     xt1nrm = dnorm2(lent,WORK(KT1),1)
c     xt2nrm = dnorm2(lent,WORK(KT2),1)
c     write(LUPRI,*) '   HTERM: After scal.: Norm of T1: ',xt1nrm
c     write(LUPRI,*) '   HTERM: After scal.: Norm of T2: ',xt2nrm
C
C------------------------------------------------------
C                    Scale amplitudes:
C                    T3(cd,l#k) <- d(cd,l) * T3(cd,l#k)
C------------------------------------------------------
C
                     DO ISYML = 1,NSYM
C
                        ISYMCD = MULD2H(ISYML,ISYMK)
                        LENCDL = NMATAB(ISYMCD)*NRHF(ISYML)
C
                        KOFFT = KT3 + IOFT3(ISYML)
                        KOFFD = KCHOV + ICKASR(ISYMCD,ISYML)
C
                        CALL CCHO_SCVEC(WORK(KOFFT),WORK(KOFFD),LENCDL,
     &                                  NUMIK)
C
                     ENDDO


c     lenm   = NCKATR(ISYMI)*NUMII
c     xm1nrm = dnorm2(lenm,WORK(KM1),1)
c     xm2nrm = dnorm2(lenm,WORK(KM2),1)
c     write(LUPRI,*) '   HTERM: Bef. VW: Norm of M1(dl,a#i): ',xm1nrm
c     write(LUPRI,*) '   HTERM: Bef. VW: Norm of M2(dl,a#i): ',xm2nrm
c     lenk   = NCKI(ISYMK)*NUMIK
c     xk1nrm = dnorm2(lenk,WORK(KK1),1)
c     xk2nrm = dnorm2(lenk,WORK(KK2),1)
c     write(LUPRI,*) '   HTERM: Bef. VW: Norm of K1(dl,j#k) : ',xk1nrm
c     write(LUPRI,*) '   HTERM: Bef. VW: Norm of K2(dl,j#k) : ',xk2nrm
C
                     DO ISYMJ = 1,NSYM
C
                        IF (NRHF(ISYMJ) .LE. 0) GOTO 997
C
C-----------------------------------------
C                       Calculate V and W.
C-----------------------------------------
C

c     if ((nsym.eq.1) .and. (numii.eq.nrhft)
c    &    .and. (numik.eq.nrhft)) then
c        write(LUPRI,*)
c        write(LUPRI,*) '   HTERM: calling H_VW_2:'
c        call h_vw_2(WORK(KK1),WORK(KK2),WORK(KM1),WORK(KM2),
c    &               WORK(KVMAT),WORK(KWMAT),NVIRT,NRHFT)
c     endif

                        ISYMDL = MULD2H(ISYMJ,ISYMK)
                        ISYMA  = MULD2H(ISYMDL,ISYMI)
C
                        NDL   = NT1AM(ISYMDL)
                        NTODL = MAX(NDL,1)
                        NAI   = NVIR(ISYMA)*NUMII
                        NTOAI = MAX(NAI,1)
                        NJK   = NRHF(ISYMJ)*NUMIK
C
                        KOFM1 = KM1 + IOFM12(ISYMA)
                        KOFM2 = KM2 + IOFM12(ISYMA)
                        KOFK1 = KK1 + IOFK12(ISYMJ)
                        KOFK2 = KK2 + IOFK12(ISYMJ)

c     write(LUPRI,*)
c     write(LUPRI,*) '   HTERM: About to calculate V and W:'
c     write(LUPRI,*) '   II1,NUMII,IK1,NUMIK           : ',
c    &               II1,NUMII,IK1,NUMIK
c     write(LUPRI,*) '   ISYMA,ISYMI,ISYMJ,ISYMK,ISYMDL: ',
c    &               ISYMA,ISYMI,ISYMJ,ISYMK,ISYMDL
c     write(LUPRI,*) '   NDL,NTODL,NAI,NTOAI,NJK       : ',
c    &               NDL,NTODL,NAI,NTOAI,NJK
c     write(LUPRI,*) '   KM1,KM2,KK1,KK2               : ',
c    &               KM1,KM2,KK1,KK2
c     write(LUPRI,*) '   KVMAT,KWMAT                   : ',
c    &               KVMAT,KWMAT
c     write(LUPRI,*) '   KOFM1,KOFM2,KOFK1,KOFK2       : ',
c    &               KOFM1,KOFM2,KOFK1,KOFK2
c     write(LUPRI,*) '   IOFM12(ISYMA),IOFK12(ISYMJ)   : ',
c    &               IOFM12(ISYMA),IOFK12(ISYMJ)
c     write(LUPRI,*) '   XMONE,ZERO,ONE,TWO            : ',
c    &               XMONE,ZERO,ONE,TWO
c     write(LUPRI,*)
C
                        CALL DGEMM('T','N',NAI,NJK,NDL,ONE,
     &                             WORK(KOFM1),NTODL,WORK(KOFK1),NTODL,
     &                             ZERO,WORK(KVMAT),NTOAI)

c     write(LUPRI,*) '   Norm of V after Coul.: ',
c    &               dnorm2(nai*njk,WORK(KVMAT),1)
C
                        CALL DCOPY(NAI*NJK,WORK(KVMAT),1,WORK(KWMAT),1)

c     write(LUPRI,*) '   Norm of W after Copy : ',
c    &               dnorm2(nai*njk,WORK(KWMAT),1)
C
                        CALL DGEMM('T','N',NAI,NJK,NDL,XMONE,
     &                             WORK(KOFM1),NTODL,WORK(KOFK2),NTODL,
     &                             TWO,WORK(KVMAT),NTOAI)

c     write(LUPRI,*) '   Norm of V after Exch.: ',
c    &               dnorm2(nai*njk,WORK(KVMAT),1)
C
                        CALL DGEMM('T','N',NAI,NJK,NDL,ONE,
     &                             WORK(KOFM2),NTODL,WORK(KOFK2),NTODL,
     &                             ONE,WORK(KWMAT),NTOAI)

c     write(LUPRI,*) '   Norm of W after Exch.: ',
c    &               dnorm2(nai*njk,WORK(KWMAT),1)

c     write(LUPRI,*)
c     write(LUPRI,*) '   HTERM: V(ai,jk):'
c     write(LUPRI,*) '   Dim. of V: ',NAI,' by ',NJK
c     write(LUPRI,*) '   Norm of V: ',dnorm2(nai*njk,WORK(KVMAT),1)
c     call output(WORK(KVMAT),1,NAI,1,NJK,NAI,NJK,1,6)
c     write(LUPRI,*)
c     write(LUPRI,*) '   HTERM: W(ai,jk):'
c     write(LUPRI,*) '   Dim. of W: ',NAI,' by ',NJK
c     write(LUPRI,*) '   Norm of W: ',dnorm2(nai*njk,WORK(KWMAT),1)
c     call output(WORK(KWMAT),1,NAI,1,NJK,NAI,NJK,1,6)
C
C-----------------------------------------
C                       Calculate P and Q.
C-----------------------------------------
C
                        ISYMB = ISYMJ
C
                        NAIK  = NVIR(ISYMA)*NUMII*NUMIK
                        NTAIK = MAX(NAIK,1)
                        NB    = NVIR(ISYMB)
                        NTOTB = MAX(NB,1)
                        NJ    = NRHF(ISYMJ)
C
                        KOFT1 = KT1 + IOFT12(ISYMB)
                        KOFT2 = KT2 + IOFT12(ISYMB)
                        KOFFT = IT1AM(ISYMB,ISYMJ) + 1
C
                        CALL DGEMM('N','N',NAIK,NJ,NB,ONE,
     &                             WORK(KOFT1),NTAIK,T1AM(KOFFT),NTOTB,
     &                             ZERO,WORK(KPMAT),NTAIK)
C
                        CALL DGEMM('N','N',NAIK,NJ,NB,XMONE,
     &                             WORK(KOFT2),NTAIK,T1AM(KOFFT),NTOTB,
     &                             ZERO,WORK(KQMAT),NTAIK)

c     write(LUPRI,*)
c     write(LUPRI,*) '   HTERM: P(aik,j):'
c     write(LUPRI,*) '   Dim. of P: ',NAIK,' by ',NJ
c     write(LUPRI,*) '   Norm of P: ',dnorm2(naik*nj,WORK(KPMAT),1)
c     call output(WORK(KPMAT),1,NAIK,1,NJ,NAIK,NJ,1,6)
c     write(LUPRI,*)
c     write(LUPRI,*) '   HTERM: Q(aik,j):'
c     write(LUPRI,*) '   Dim. of Q: ',NAIK,' by ',NJ
c     write(LUPRI,*) '   Norm of Q: ',dnorm2(naik*nj,WORK(KQMAT),1)
c     call output(WORK(KQMAT),1,NAIK,1,NJ,NAIK,NJ,1,6)
c     write(LUPRI,*)
C
C-----------------------------------------
C                       Calculate X and R.
C-----------------------------------------
C
                        ISYMCD = MULD2H(ISYMJ,ISYMI)
                        ISYML  = MULD2H(ISYMCD,ISYMK)
                        ISYMB  = ISYML
C
                        NCD   = NMATAB(ISYMCD)
                        NTOCD = MAX(NCD,1)
                        NLK   = NRHF(ISYML)*NUMIK
                        NTOLK = MAX(NLK,1)
                        NJI   = NRHF(ISYMJ)*NUMII
                        NB    = NVIR(ISYMB)
                        NTOTB = MAX(NB,1)
                        NL    = NRHF(ISYML)
                        NTOTL = MAX(NL,1)
                        NKJI  = NUMIK*NJI
C
                        KOFT3 = KT3 + IOFT3(ISYML)
                        KOFM3 = KM3 + IOFM3(ISYMJ)
                        KOFK3 = KK3 + IOFK3(ISYMJ)
                        KOFFT = IT1AM(ISYMB,ISYML) + 1
C
                        CALL DGEMM('T','N',NLK,NJI,NCD,ONE,
     &                             WORK(KOFT3),NTOCD,WORK(KOFM3),NTOCD,
     &                             ZERO,WORK(KXMAT),NTOLK)
C
                        CALL DGEMM('T','N',NL,NKJI,NB,XMONE,
     &                             T1AM(KOFFT),NTOTB,WORK(KOFK3),NTOTB,
     &                             ZERO,WORK(KRMAT),NTOTL)
C
C-----------------------------------
C                       Calculate Y.
C-----------------------------------
C
                        IF (ISYMK .EQ. ISYMI) THEN
C
                           NJ = NRHF(ISYMJ)
C
                           DO I = 1,NUMII
                              DO K = 1,NUMIK
                                 KOFFY = KYMAT + NUMIK*(I - 1) + K - 1
                                 DO J = 1,NJ
                                    KOFFX = KXMAT
     &                                    + NJ*NUMIK*NJ*(I - 1)
     &                                    + NJ*NUMIK*(J - 1)
     &                                    + NJ*(K - 1)
     &                                    + J - 1
                                    WORK(KOFFY) = WORK(KOFFY)
     &                                          + WORK(KOFFX)
                                 ENDDO
                              ENDDO
                           ENDDO
C
                        ENDIF
C
C--------------------------------------------------------
C                       Calculate energy contribution XR.
C--------------------------------------------------------
C
                        NLKJI = NRHF(ISYML)*NUMIK*NRHF(ISYMJ)*NUMII
C
                        HCOR = HCOR
     &                       + DDOT(NLKJI,WORK(KXMAT),1,WORK(KRMAT),1)

c     exr = exr + DDOT(NLKJI,WORK(KXMAT),1,WORK(KRMAT),1)
C
C----------------------------------------------------------------
C                       Calculate energy contributions VP and WQ.
C----------------------------------------------------------------
C
                        NAI = NVIR(ISYMA)*NUMII
C
                        DO K = 1,NUMIK
                           DO J = 1,NRHF(ISYMJ)
C
                              JK = NRHF(ISYMJ)*(K - 1) + J
                              KJ = NUMIK*(J - 1) + K
C
                              KOFFV = KVMAT + NAI*(JK - 1)
                              KOFFP = KPMAT + NAI*(KJ - 1)
C
                              HCOR = HCOR
     &                             + DDOT(NAI,WORK(KOFFV),1,
     &                                        WORK(KOFFP),1)

c     evp = evp + DDOT(NAI,WORK(KOFFV),1,WORK(KOFFP),1)
C
                              KOFFW = KWMAT + NAI*(JK - 1)
                              KOFFQ = KQMAT + NAI*(KJ - 1)
C
                              HCOR = HCOR
     &                             + DDOT(NAI,WORK(KOFFW),1,
     &                                        WORK(KOFFQ),1)

c     ewq = ewq + DDOT(NAI,WORK(KOFFW),1,WORK(KOFFQ),1)
C
                           ENDDO
                        ENDDO
C
  997                   CONTINUE
C
                     ENDDO
C
C---------------------------------------------------------------
C                    Calculate energy contribution from Y and S.
C---------------------------------------------------------------
C
                     IF (ISYMK .EQ. ISYMI) THEN
C
                        NEM   = NT1AM(1)
                        NTOEM = MAX(NEM,1)
C
                        CALL DGEMV('T',NEM,LENKI,ONE,WORK(KK4),NTOEM,
     &                             T1AM,1,ZERO,WORK(KSMAT),1)
C
                        HCOR = HCOR + DDOT(LENKI,WORK(KYMAT),1,
     &                                           WORK(KSMAT),1)

c     eys = eys + DDOT(LENKI,WORK(KYMAT),1,WORK(KSMAT),1)
C
                     ENDIF
C
C--------------------------------------
C                    Check convergence.
C--------------------------------------
C
                     EH = EH + HCOR
C 
              ENERGH(ICHO) = ENERGH(ICHO) + HCOR
C 
C Decommented by Domenico
                    TNOW = SECOND()
                    DELTAT = TNOW - TLAST
                    TLAST = TNOW
                    SCNDSH(ICHO) = SCNDSH(ICHO) 
     &                          + DELTAT
C Decommented by Domenico
C 
C
c                    IF (PRINT) THEN
c                       WRITE(LUPRI,'(15X,A,I3,A,/,15X,A)')
c    &                  'Status after Cholesky vector',ICHO,':',
c    &                  '--------------------------------'
c                       IF (ABS(HCOR) .LT. THRCHO) THEN
c                          WRITE(LUPRI,'(15X,A)') 'H term converged'
c                       ELSE
c                          WRITE(LUPRI,'(15X,A)') 'H term not converged'
c                       ENDIF
c                       TIM = SECOND() - TIMT
c                       WRITE(LUPRI,'(15X,A,F10.2,A,/)')
c    &                  'Accumulated H-time: ',TIM,' seconds'
c                    ENDIF
C
                     IF (DABS(HCOR) .LT. THRCHO) GOTO 998

c     write(LUPRI,*) '   HTERM: Final: ICHO,HCOR,EH: ',ICHO,HCOR,EH
C
                  ENDDO
C
  998             CONTINUE
C
               ENDDO
C
  999          CONTINUE
C
            ENDDO
C
         ENDDO
C
 1000    CONTINUE
C
      ENDDO

c     lenk  = NTRAOC(1)
c     liajb = NT2AM(1)
c     lent2 = NT2SQ(1)
c     lent1 = NT1AM(1)
c     xknorm = dnorm2(lenk,XKINT,1)
c     xniajb = dnorm2(liajb,XIAJB,1)
c     t2norm = dnorm2(lent2,T2VO,1)
c     t1norm = dnorm2(lent1,T1AM,1)
c     if (lwork .ge. lent2) then
c        call cc_t2sq(XIAJB,WORK,1)
c        xsiajb = dnorm2(lent2,WORK,1)
c     else
c        write(LUPRI,*) '   HTERM: insuf. core for squaring L(iajb)'
c        write(LUPRI,*) '   HTERM: LWORK = ',LWORK,'. Need = ',lent2
c        xsiajb = 1.0D9
c     endif
c     if (lwork .ge. lenk) then
c        call dzero(WORK,lenk)
c        call daxpy(lenk,2.0D0,XKINT,1,WORK,1)
c        DO ISYME = 1,NSYM
c           ISYIKM = ISYME
c           DO E = 1,NVIR(ISYME)
c              DO ISYMM = 1,NSYM
c                 ISYMIK = MULD2H(ISYMM,ISYIKM)
c                 DO M = 1,NRHF(ISYMM)
c                    DO ISYMK = 1,NSYM
c                       ISYMI  = MULD2H(ISYMK,ISYMIK)
c                       ISYMMK = MULD2H(ISYMK,ISYMM)
c                       DO K = 1,NRHF(ISYMK)
c                          DO I = 1,NRHF(ISYMI)
c                             KIKME = ISJIKA(ISYIKM,ISYME)
c    &                              + NMAJIK(ISYIKM)*(E - 1)
c    &                              + ISJIK(ISYMIK,ISYMM)
c    &                              + NMATIJ(ISYMIK)*(M - 1)
c    &                              + IMATIJ(ISYMI,ISYMK)
c    &                              + NRHF(ISYMI)*(K - 1) + I
c                             KMKIE = ISJIKA(ISYIKM,ISYME)
c    &                              + NMAJIK(ISYIKM)*(E - 1)
c    &                              + ISJIK(ISYMMK,ISYMI)
c    &                              + NMATIJ(ISYMMK)*(I - 1)
c    &                              + IMATIJ(ISYMM,ISYMK)
c    &                              + NRHF(ISYMM)*(K - 1) + M
c                             WORK(KMKIE) = WORK(KMKIE) - XKINT(KIKME)
c                          ENDDO
c                       ENDDO
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c        xk2cme = dnorm2(lenk,WORK,1)
c     else
c        write(LUPRI,*) '   HTERM: insuf. core for K 2CME'
c        write(LUPRI,*) '   HTERM: LWORK = ',LWORK,'. Need = ',lenk
c        xk2cme = 1.0D9
c     endif
c     write(LUPRI,*)
c     write(LUPRI,*) '   Exiting HTERM:'
c     write(LUPRI,*) '   =============='
c     write(LUPRI,*) '   LWORK = ',LWORK
c     write(LUPRI,*) '   Norm of K(lj,k;d) : ',xknorm
c     write(LUPRI,*) '   Norm of K 2CME    : ',xk2cme
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xniajb,' (packed)'
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xsiajb,' (squared)'
c     write(LUPRI,*) '   Norm of T1AM      : ',t1norm
c     write(LUPRI,*) '   Norm of T2VO      : ',t2norm
c     write(LUPRI,*) '   Contribution EVP  : ',evp
c     write(LUPRI,*) '   Contribution EWQ  : ',ewq
c     write(LUPRI,*) '   Contribution EXR  : ',exr
c     write(LUPRI,*) '   Contribution EYS  : ',eys
c     write(LUPRI,*)
C
      RETURN
      END
C  /* Deck ccho_hxt12 */
      SUBROUTINE CCHO_HXT12(T2VO,T1,T2,IK1,NUMIK,ISYMK,II1,NUMII,ISYMI,
     &                      IOFT12)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract amplitudes:
C
C     T1(a#i,#kb) = T(ai,bk) = T2VO(ab,ik), ordering: IOFT12
C     T2(a#i,#kb) = T(ak,bi) = T2VO(ab,ki), ordering: IOFT12
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2VO(*),T1(*),T2(*)
      INTEGER IOFT12(8)
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYMIK = MULD2H(ISYMK,ISYMI)
      ISYMKI = ISYMIK
      ISYMAB = ISYMIK
C
C-------------------
C     Set up IOFT12.
C-------------------
C
      ICOUNT = 0
C
      DO ISYMB = 1,NSYM
C
         ISYMA = MULD2H(ISYMB,ISYMAB)
C
         IOFT12(ISYMB) = ICOUNT
C
         ICOUNT = ICOUNT + NVIR(ISYMA)*NUMII*NUMIK*NVIR(ISYMB)
C
      ENDDO

c     write(LUPRI,*) '      HXT12: II1,NUMII,ISYMI: ',II1,NUMII,ISYMI
c     write(LUPRI,*) '      HXT12: IK1,NUMIK,ISYMK: ',IK1,NUMIK,ISYMK
c     write(LUPRI,*) '      HXT12: IOFT12 (end)   : ',
c    &                        (IOFT12(JSYM),JSYM=1,NSYM),ICOUNT
C
C-------------
C     Extract.
C-------------
C
      DO K = 1,NUMIK
C
         KK = IK1 + K - 1
C
         DO I = 1,NUMII
C
            II = II1 + I - 1
C
            IK = IMATIJ(ISYMI,ISYMK) + NRHF(ISYMI)*(KK - 1) + II
            KI = IMATIJ(ISYMK,ISYMI) + NRHF(ISYMK)*(II - 1) + KK
C
            DO ISYMB = 1,NSYM
C
               ISYMA = MULD2H(ISYMB,ISYMAB)
C
               IF (NVIR(ISYMA) .GT. 0) THEN
C
                  LENAI  = NVIR(ISYMA)*NUMII
                  LENAIK = LENAI*NUMIK
C
                  DO B = 1,NVIR(ISYMB)
C
                     KAB = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + 1
C
                     KOFF1 = IT2VO(ISYMAB,ISYMIK)
     &                     + NMATAB(ISYMAB)*(IK - 1)
     &                     + KAB
                     KOFF2 = IT2VO(ISYMAB,ISYMKI)
     &                     + NMATAB(ISYMAB)*(KI - 1)
     &                     + KAB
                     KOFFT = IOFT12(ISYMB) + LENAIK*(B - 1)
     &                     + LENAI*(K - 1)
     &                     + NVIR(ISYMA)*(I - 1) + 1
C
                     CALL DCOPY(NVIR(ISYMA),T2VO(KOFF1),1,T1(KOFFT),1)
                     CALL DCOPY(NVIR(ISYMA),T2VO(KOFF2),1,T2(KOFFT),1)
C
                  ENDDO
C
               ENDIF
C
            ENDDO
C
         ENDDO
C
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_hxt3 */
      SUBROUTINE CCHO_HXT3(T2VO,T3,IK1,NUMIK,ISYMK,IOFT3)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract amplitudes:
C
C     T3(cd,l#k) = T(ck,dl) = T2VO(cd,kl), ordering: IOFT3.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2VO(*),T3(*)
      INTEGER IOFT3(8)
#include <ccorb.h>
#include <ccsdsym.h>
C
C------------------
C     Set up IOFT3.
C------------------
C
      ICOUNT = 0
C
      DO ISYML = 1,NSYM
C
         ISYMCD = MULD2H(ISYML,ISYMK)
C
         IOFT3(ISYML) = ICOUNT
C
         ICOUNT = ICOUNT + NMATAB(ISYMCD)*NRHF(ISYML)*NUMIK
C
      ENDDO
C
C-------------
C     Extract.
C-------------
C
      DO ISYMKL = 1,NSYM
C
         ISYML  = MULD2H(ISYMKL,ISYMK)
         ISYMCD = ISYMKL
C
         DO L = 1,NRHF(ISYML)
            DO K = 1,NUMIK
C
               IK = IK1 + K - 1
               KL = IMATIJ(ISYMK,ISYML) + NRHF(ISYMK)*(L - 1) + IK
               LK = NRHF(ISYML)*(K - 1) + L
C
               KOFF1 = IT2VO(ISYMCD,ISYMKL) + NMATAB(ISYMCD)*(KL - 1)
     &               + 1
               KOFF2 = IOFT3(ISYML) + NMATAB(ISYMCD)*(LK - 1) + 1
C
               CALL DCOPY(NMATAB(ISYMCD),T2VO(KOFF1),1,T3(KOFF2),1)
C
            ENDDO
         ENDDO
C
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_hxk4 */
      SUBROUTINE CCHO_HXK4(XKINT,XK4,IK1,NUMIK,ISYMK,II1,NUMII,ISYMI)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract integrals:
C
C     K4(bj,#k#i) = 2*(bj|ki)  - (bi|kj)
C                 = 2*K(jki,b) - K(ikj,b).
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XKINT(*),XK4(*)
      INTEGER BJKI
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (TWO = 2.00D0)
C
      ISYMIK = MULD2H(ISYMK,ISYMI)
      ISYMBJ = ISYMIK
C
      DO ISYMB = 1,NSYM
C
         ISYMJ  = MULD2H(ISYMB,ISYMBJ)
         ISYMJK = MULD2H(ISYMJ,ISYMK)
         ISYJKI = MULD2H(ISYMJK,ISYMI)
         ISYIKJ = ISYJKI
C
         IF (NRHF(ISYMJ) .GT. 0) THEN
C
            DO B = 1,NVIR(ISYMB)
               DO I = 1,NUMII
C
                  II = II1 + I - 1
C
                  DO K = 1,NUMIK
C
                     IK = IK1 + K - 1
                     KI = NUMIK*(I - 1) + K
C
                     DO J = 1,NRHF(ISYMJ)
C
                        JKIB = ISJIKA(ISYJKI,ISYMB)
     &                       + NMAJIK(ISYJKI)*(B - 1)
     &                       + ISJIK(ISYMJK,ISYMI)
     &                       + NMATIJ(ISYMJK)*(II - 1)
     &                       + IMATIJ(ISYMJ,ISYMK)
     &                       + NRHF(ISYMJ)*(IK - 1) + J
                        IKJB = ISJIKA(ISYIKJ,ISYMB)
     &                       + NMAJIK(ISYIKJ)*(B - 1)
     &                       + ISJIK(ISYMIK,ISYMJ)
     &                       + NMATIJ(ISYMIK)*(J - 1)
     &                       + IMATIJ(ISYMI,ISYMK)
     &                       + NRHF(ISYMI)*(IK - 1) + II
C
                        BJKI = NT1AM(ISYMBJ)*(KI - 1)
     &                       + IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
     &                       + B
C
                        XK4(BJKI) = TWO*XKINT(JKIB) - XKINT(IKJB)
C
                     ENDDO
C
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDIF
C
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_hxk3 */
      SUBROUTINE CCHO_HXK3(XKINT,XK3,IK1,NUMIK,ISYMK,II1,NUMII,ISYMI,
     &                     IOFK3)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract integrals:
C
C     K3(b,#kj#i) = (bj,ki) = K(jki,b), ordering: IOFK3.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XKINT(*),XK3(*)
      INTEGER IOFK3(8)
      INTEGER BKJI
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYMBJ = MULD2H(ISYMK,ISYMI)
C
C------------------
C     Set up IOFK3.
C------------------
C
      ICOUNT = 0
C
      DO ISYMJ = 1,NSYM
C
         ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
         IOFK3(ISYMJ) = ICOUNT
C
         ICOUNT = ICOUNT + NVIR(ISYMB)*NUMIK*NRHF(ISYMJ)*NUMII
C
      ENDDO
C
C-------------
C     Extract.
C-------------
C
      DO ISYMB = 1,NSYM
C
         ISYMJ  = MULD2H(ISYMB,ISYMBJ)
         ISYMJK = MULD2H(ISYMJ,ISYMK)
         ISYJKI = MULD2H(ISYMJK,ISYMI)
C
         IF (NRHF(ISYMJ) .GT. 0) THEN
C
            DO B = 1,NVIR(ISYMB)
               DO I = 1,NUMII
C
                  II = II1 + I - 1
C
                  DO K = 1,NUMIK
C
                     IK = IK1 + K - 1
C
                     DO J = 1,NRHF(ISYMJ)
C
                        JKIB = ISJIKA(ISYJKI,ISYMB)
     &                       + NMAJIK(ISYJKI)*(B - 1)
     &                       + ISJIK(ISYMJK,ISYMI)
     &                       + NMATIJ(ISYMJK)*(II - 1)
     &                       + IMATIJ(ISYMJ,ISYMK)
     &                       + NRHF(ISYMJ)*(IK - 1) + J
C
                        KJI  = NUMIK*NRHF(ISYMJ)*(I - 1)
     &                       + NUMIK*(J - 1) + K
                        BKJI = IOFK3(ISYMJ) + NVIR(ISYMB)*(KJI - 1) + B
C
                        XK3(BKJI) = XKINT(JKIB)
C
                     ENDDO
C
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDIF
C
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_hxk12 */
      SUBROUTINE CCHO_HXK12(XKINT,XK1,XK2,IK1,NUMIK,ISYMK,IOFK12)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract integrals:
C
C     K1(dl,j#k) = (dl|kj) = K(lkj,d), ordering: IOFK12
C     K2(dl,j#k) = (dj|kl) = K(jkl,d), ordering: IOFK12
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XKINT(*),XK1(*),XK2(*)
      INTEGER IOFK12(8)
#include <ccorb.h>
#include <ccsdsym.h>
      INTEGER DL,DLJK
C
C-------------------
C     Set up IOFK12.
C-------------------
C
      ICOUNT = 0
C
      DO ISYMJ = 1,NSYM
C
         ISYMDL = MULD2H(ISYMJ,ISYMK)
C
         IOFK12(ISYMJ) = ICOUNT
C
         ICOUNT = ICOUNT + NT1AM(ISYMDL)*NRHF(ISYMJ)*NUMIK
C
      ENDDO

c     write(LUPRI,*) '      HXK12: IK1,NUMIK,ISYMK: ',IK1,NUMIK,ISYMK
c     write(LUPRI,*) '      HXK12: IOFK12 (end): ',
c    &                        (IOFK12(JSYM),JSYM=1,NSYM),ICOUNT
C
C-------------
C     Extract.
C-------------
C
c     write(LUPRI,*) 'NTRAOC(1): ',NTRAOC(1)

      DO ISYMJ = 1,NSYM
C
         IF (NRHF(ISYMJ) .GT. 0) THEN
C
            ISYMJK = MULD2H(ISYMJ,ISYMK)
            ISYMDL = ISYMJK
C
            DO K = 1,NUMIK
C
               IK = IK1 + K - 1
C
               DO J = 1,NRHF(ISYMJ)
C
                  JK  = NRHF(ISYMJ)*(K - 1) + J
C
                  DO ISYML = 1,NSYM
C
                     ISYMD  = MULD2H(ISYML,ISYMDL)
                     ISYMLK = MULD2H(ISYML,ISYMK)
                     ISYLKJ = ISYMD
                     ISYJKL = ISYMD
C
                     DO L = 1,NRHF(ISYML)
C
                        LKJ = ISJIK(ISYMLK,ISYMJ)
     &                      + NMATIJ(ISYMLK)*(J - 1)
     &                      + IMATIJ(ISYML,ISYMK)
     &                      + NRHF(ISYML)*(IK - 1) + L
                        JKL = ISJIK(ISYMJK,ISYML)
     &                      + NMATIJ(ISYMJK)*(L - 1)
     &                      + IMATIJ(ISYMJ,ISYMK)
     &                      + NRHF(ISYMJ)*(IK - 1) + J
C
                        DO D = 1,NVIR(ISYMD)
C
                           DLJK = IOFK12(ISYMJ) + NT1AM(ISYMDL)*(JK - 1)
     &                          + IT1AM(ISYMD,ISYML)
     &                          + NVIR(ISYMD)*(L - 1) + D
                           LKJD = ISJIKA(ISYLKJ,ISYMD)
     &                          + NMAJIK(ISYLKJ)*(D - 1) + LKJ
                           JKLD = ISJIKA(ISYJKL,ISYMD)
     &                          + NMAJIK(ISYJKL)*(D - 1) + JKL

c     if (JKLD.GT.NTRAOC(1)) then
c        write(LUPRI,*) '...JKLD = ',JKLD,' but NTRAOC(1) = ',NTRAOC(1)
c        CALL QUIT(' jkld out of bounds in HXK12 ')
c     endif
C
                           XK1(DLJK) = XKINT(LKJD)
                           XK2(DLJK) = XKINT(JKLD)

c     write(LUPRI,*) 'HXK12: DLJK,LKJD,JKLD,K1,K2: ',
c    &           DLJK,LKJD,JKLD,XK1(DLJK),XK2(DLJK)
C
                        ENDDO
C
                     ENDDO
C
                  ENDDO
C
               ENDDO
C
            ENDDO
C
         ENDIF
C
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     ISYDLJ = ISYMK
c     DO ISYMJ = 1,NSYM
c        ISYMDL = MULD2H(ISYMJ,ISYDLJ)
c        DO K = 1,NUMIK
c           DO J = 1,NRHF(ISYMJ)
c              DO ISYML = 1,NSYM
c                 ISYMD = MULD2H(ISYML,ISYMDL)
c                 ISYMDJ = MULD2H(ISYMD,ISYMJ)
c                 DO L = 1,NRHF(ISYML)
c                    DO D = 1,NVIR(ISYMD)
c                       KDLJK = IOFK12(ISYMJ)
c    &                        + NT1AM(ISYMDL)*NRHF(ISYMJ)*(K - 1)
c    &                        + NT1AM(ISYMDL)*(J - 1)
c    &                        + IT1AM(ISYMD,ISYML)
c    &                        + NVIR(ISYMD)*(L - 1) + D
c                       KDJLK = IOFK12(ISYML)
c    &                        + NT1AM(ISYMDJ)*NRHF(ISYML)*(K - 1)
c    &                        + NT1AM(ISYMDJ)*(L - 1)
c    &                        + IT1AM(ISYMD,ISYMJ)
c    &                        + NVIR(ISYMD)*(J - 1) + D
c                       NTST  = NTST + 1
c                       TEST  = XK1(KDLJK) - XK2(KDJLK)
c                       IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*) '      HXK12: ISYMD,ISYML,ISYMJ,ISYMK: ',
c    &                         ISYMD,ISYML,ISYMJ,ISYMK
c     write(LUPRI,*) '      HXK12: D,L,J,K     : ',D,L,J,IK1+K-1
c     write(LUPRI,*) '      HXK12: DLJK,Diff.  : ',KDLJK,TEST
c                          ICOUNE = ICOUNE + 1
c                       ENDIF
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     ndim = NCKI(ISYMK)*NUMIK
c     write(LUPRI,*) '      HXK12: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested (dim = ',ndim,')'
c     if ((nsym.eq.1) .and. (numik.eq.nrhft)) then
c        NERR1 = 0
c        NERR2 = 0
c        DO K = 1,NRHFT
c           DO J = 1,NRHFT
c              DO L = 1,NRHFT
c                 DO D = 1,NVIRT
c                    KOFK1 = NVIRT*NRHFT*NRHFT*(K - 1)
c    &                     + NVIRT*NRHFT*(J - 1)
c    &                     + NVIRT*(L - 1)
c    &                     + D
c                    LKJD  = NRHFT*NRHFT*NRHFT*(D - 1)
c    &                     + NRHFT*NRHFT*(J - 1)
c    &                     + NRHFT*(K - 1)
c    &                     + L
c                    TEST  = XK1(KOFK1) - XKINT(LKJD)
c                    IF (DABS(TEST) .GT. 1.0D-15) THEN
c                       WRITE(LUPRI,*) '      HXK12: K1 is jodido!'
c                       WRITE(LUPRI,*) '      HXK12: D,L,K,J,TEST: ',
c    &                                           D,L,K,J,TEST
c                       NERR1 = NERR1 + 1
c                    ENDIF
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c        WRITE(LUPRI,*) '     HXK12: Detected ',NERR1,' K1 errors'
c        DO K = 1,NRHFT
c           DO J = 1,NRHFT
c              DO L = 1,NRHFT
c                 DO D = 1,NVIRT
c                    KOFK2 = NVIRT*NRHFT*NRHFT*(K - 1)
c    &                     + NVIRT*NRHFT*(J - 1)
c    &                     + NVIRT*(L - 1)
c    &                     + D
c                    JKLD  = NRHFT*NRHFT*NRHFT*(D - 1)
c    &                     + NRHFT*NRHFT*(L - 1)
c    &                     + NRHFT*(K - 1)
c    &                     + J
c                    TEST  = XK2(KOFK2) - XKINT(JKLD)
c                    IF (DABS(TEST) .GT. 1.0D-15) THEN
c                       WRITE(LUPRI,*) '      HXK12: K1 is jodido!'
c                       WRITE(LUPRI,*) '      HXK12: D,L,K,J,TEST: ',
c    &                                           D,L,K,J,TEST
c                       NERR2 = NERR2 + 1
c                    ENDIF
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c        WRITE(LUPRI,*) '     HXK12: Detected ',NERR2,' K2 errors'
c     endif
C
      RETURN
      END
C  /* Deck ccho_hxm123 */
      SUBROUTINE CCHO_HXM123(XIAJB,XM1,XM2,XM3,II1,NUMII,ISYMI,IOFM12,
     &                       IOFM3)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Extract integrals:
C
C     M1(bj,a#i) = L(jb,ia) = X(bj,ai), ordering: IOFM12
C     M2(bj,a#i) = L(ja,ib) = X(bi,aj), ordering: IOFM12
C     M3(ab,j#i) = L(ia,jb) = X(bj,ai), ordering: IOFM3
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XIAJB(*),XM1(*),XM2(*),XM3(*)
      INTEGER IOFM12(8),IOFM3(8)
      INTEGER AI,BJ,BI,AJ,BJAI,BIAJ,AB
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C-------------------------
C     Set up index arrays.
C-------------------------
C
      ICOUN1 = 0
      ICOUN3 = 0
C
      DO ISYM2 = 1,NSYM
C
         ISYM1 = MULD2H(ISYM2,ISYMI)
C
         IOFM12(ISYM2) = ICOUN1
         IOFM3(ISYM2)  = ICOUN3
C
         ICOUN1 = ICOUN1 + NT1AM(ISYM1)*NVIR(ISYM2)*NUMII
         ICOUN3 = ICOUN3 + NMATAB(ISYM1)*NRHF(ISYM2)*NUMII
C
      ENDDO

c     write(LUPRI,*) '      HXM123: II1,NUMII,ISYMI: ',II1,NUMII,ISYMI
c     write(LUPRI,*) '      HXM123: IOFM12 (end): ',
c    &                          (IOFM12(JSYM),JSYM=1,NSYM),ICOUN1
c     write(LUPRI,*) '      HXM123: IOFM3  (end): ',
c    &                          (IOFM3(JSYM),JSYM=1,NSYM),ICOUN3
c     if (icoun1 .ne. icoun3) then
c        write(LUPRI,*) '      HXM123: ICOUN1 .ne. ICOUN3: ',ICOUN1,ICOUN3
c        CALL QUIT(' hxm123: jodido ')
c     endif
C
C-------------
C     Extract.
C-------------
C
      DO ISYMAI = 1,NSYM
C
         ISYMBJ = ISYMAI
         ISYMA  = MULD2H(ISYMAI,ISYMI)
C
         DO I = 1,NUMII
C
            II = II1 + I - 1
C
            DO A = 1,NVIR(ISYMA)
C
               AI  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(II - 1) + A
               KAI = NVIR(ISYMA)*(I - 1) + A
C
               DO ISYMJ = 1,NSYM
C
                  ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                  ISYMAJ = MULD2H(ISYMJ,ISYMA)
                  ISYMBI = ISYMAJ
                  ISYMAB = MULD2H(ISYMB,ISYMA)
C
                  DO J = 1,NRHF(ISYMJ)
C
                     KJI = NRHF(ISYMJ)*(I - 1) + J
                     AJ  = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J - 1) + A
C
                     DO B = 1,NVIR(ISYMB)
C
                        BJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
     &                     + B
                        BI = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(II - 1)
     &                     + B
                        AB = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1)
     &                     + A
C
                        BJAI = IT2AM(ISYMBJ,ISYMAI) + INDEX(BJ,AI)
                        BIAJ = IT2AM(ISYMBI,ISYMAJ) + INDEX(BI,AJ)
C
                        KM12 = IOFM12(ISYMA) + NT1AM(ISYMBJ)*(KAI - 1)
     &                       + BJ
                        KM3  = IOFM3(ISYMJ)  + NMATAB(ISYMAB)*(KJI - 1)
     &                       + AB
C
                        XM1(KM12) = XIAJB(BJAI)
                        XM3(KM3)  = XIAJB(BJAI)
                        XM2(KM12) = XIAJB(BIAJ) 
C
                     ENDDO
C
                  ENDDO
C
               ENDDO
C
            ENDDO
C
         ENDDO
C
      ENDDO

C DEBUG SECTION:
C***************
c     ICOUNE = 0
c     NTST   = 0
c     ISYBJA = ISYMI
c     DO ISYMA = 1,NSYM
c        ISYMBJ = MULD2H(ISYMA,ISYBJA)
c        DO I = 1,NUMII
c           DO A = 1,NVIR(ISYMA)
c              DO ISYMJ = 1,NSYM
c                 ISYMB  = MULD2H(ISYMJ,ISYMBJ)
c                 ISYMAJ = MULD2H(ISYMJ,ISYMA)
c                 DO J = 1,NRHF(ISYMJ)
c                    DO B = 1,NVIR(ISYMB)
c                       KBJAI = IOFM12(ISYMA)
c    &                        + NT1AM(ISYMBJ)*NVIR(ISYMA)*(I - 1)
c    &                        + NT1AM(ISYMBJ)*(A - 1)
c    &                        + IT1AM(ISYMB,ISYMJ)
c    &                        + NVIR(ISYMB)*(J - 1)
c    &                        + B
c                       KAJBI = IOFM12(ISYMB)
c    &                        + NT1AM(ISYMAJ)*NVIR(ISYMB)*(I - 1)
c    &                        + NT1AM(ISYMAJ)*(B - 1)
c    &                        + IT1AM(ISYMA,ISYMJ)
c    &                        + NVIR(ISYMA)*(J - 1)
c    &                        + A
c                       NTST  = NTST + 1
c                       TEST  = XM1(KBJAI) - XM2(KAJBI)
c                       IF (DABS(TEST) .GT. 1.0D-15) THEN
c     write(LUPRI,*) '      HXM12: ISYMB,ISYMJ,ISYMA,ISYMI: ',
c    &                         ISYMB,ISYMJ,ISYMA,ISYMI
c     write(LUPRI,*) '      HXM12: B,J,A,I: ',B,J,A,II1+I-1
c     write(LUPRI,*) '      HXM12: Diff.  : ',TEST
c                          ICOUNE = ICOUNE + 1
c                       ENDIF
c                    ENDDO
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
c     ndim = NCKATR(ISYMI)*NUMII
c     write(LUPRI,*) '      HXM12: ',ICOUNE,' interchange errors out of ',
c    &           NTST,' tested (dim = ',ndim,')'
c     if ((nsym.eq.1) .and. (numii.eq.nrhft)) then
c        NERR1 = 0
c        NERR2 = 0
c        DO NAI = 1,NT1AMX
c           DO NBJ = 1,NT1AMX
c              KOFM1 = NT1AMX*(NAI - 1) + NBJ
c              TEST  = XM1(KOFM1) - XIAJB(INDEX(NAI,NBJ))
c              IF (DABS(TEST) .GT. 1.0D-15) THEN
c                 write(LUPRI,*) '      HXM12: XM1 is jodido!'
c                 write(LUPRI,*) '      HXM12: NAI,NBJ,TEST: ',NAI,NBJ,TEST
c                 NERR1 = NERR1 + 1
c              ENDIF
c           ENDDO
c        ENDDO
c        write(LUPRI,*) '      HXM12: Errors detected for M1: ',NERR1
c        DO I = 1,NRHFT
c           DO A = 1,NVIRT
c              NAI = NVIRT*(I - 1) + A
c              DO J = 1,NRHFT
c                 NAJ = NVIRT*(J - 1) + A
c                 DO B = 1,NVIRT
c                    NBJ = NVIRT*(J - 1) + B
c                    NBI = NVIRT*(I - 1) + B
c                    KOFM2 = NT1AMX*(NAI - 1) + NBJ
c                    TEST  = XM2(KOFM2) - XIAJB(INDEX(NAJ,NBI))
c                    IF (DABS(TEST) .GT. 1.0D-15) THEN
c                       write(LUPRI,*) '      HXM12: XM2 is jodido!'
c                       write(LUPRI,*)'      HXM12: NAI,NBJ,NAJ,NBI,TEST: ',
c    &                                    NAI,NBJ,NAJ,NBI,TEST
c                       NERR2 = NERR2 + 1
c                    ENDIF
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c        write(LUPRI,*) '      HXM12: Errors detected for M2: ',NERR2
c     endif
C
      RETURN
      END
C  /* Deck ccho_h1term */
      SUBROUTINE CCHO_H1TERM(FOCKD,T1AM,XKINT,T2VO,XIAJB,WORK,LWORK,
     &                       CHOELE,NUMCHO,EH1,FBATCH,PRINT,
     &                       NONA,NONI)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Calculate the H1 term:
C
C        EH1 = EH1 + Sum(ai) Z(ai) * U(ai)
C
C     where
C
C        Z(ai)    =   Sum(bj)  d(bja) * s(ai,bj) * t1(bj)
C
C        U(ai)    =   Sum(ckl) d(ckl) * L(kcla) * (ck|il)
C
C     and s(ai,bj) = 2 * t(ai,bj) - t(aj,bi), d(aij) and d(aib) denote
C     the occupied and virtual parts of the Cholesky decomposition of the
C     orbital energy denominator, and L(iajb) = 2 * (ia|jb) - (ja|ib) etc.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION FOCKD(*),T1AM(*)
      DIMENSION XKINT(*),T2VO(*),XIAJB(*)
      DIMENSION WORK(LWORK)
      DIMENSION CHOELE(*)
      LOGICAL   FBATCH,PRINT
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.00D0, ZERO = 0.00D0)
      PARAMETER (ONE   =  1.00D0, TWO  = 2.00D0)
C 
      TLAST = SECOND()
C
C
      IF (PRINT) THEN
         TIMT = SECOND()
         WRITE(LUPRI,'(A,/,A,/)')
     &   'Calculation of the H1 term:',
     &   '==========================='
      ENDIF

c     lenk  = NTRAOC(1)
c     liajb = NT2AM(1)
c     lent2 = NT2SQ(1)
c     lent1 = NT1AM(1)
c     xknorm = dnorm2(lenk,XKINT,1)
c     xniajb = dnorm2(liajb,XIAJB,1)
c     t2norm = dnorm2(lent2,T2VO,1)
c     t1norm = dnorm2(lent1,T1AM,1)
c     if (lwork .ge. lent2) then
c        call cc_t2sq(XIAJB,WORK,1)
c        xsiajb = dnorm2(lent2,WORK,1)
c        call dzero(WORK,lent2)
c        call daxpy(lent2,2.0D0,T2VO,1,WORK,1)
c        DO ISYMJI = 1,NSYM
c           ISYMAB = ISYMJI
c           DO ISYMI = 1,NSYM
c              ISYMJ = MULD2H(ISYMI,ISYMJI)
c              DO I = 1,NRHF(ISYMI)
c                 DO J = 1,NRHF(ISYMJ)
c                    JI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I - 1) + J
c                    IJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
c                    KABIJ = IT2VO(ISYMAB,ISYMJI)
c    &                     + NMATAB(ISYMAB)*(IJ - 1) + 1
c                    KABJI = IT2VO(ISYMAB,ISYMJI)
c    &                     + NMATAB(ISYMAB)*(JI - 1) + 1
c                    CALL DAXPY(NMATAB(ISYMAB),-1.0D0,T2VO(KABJI),1,
c    &                                                WORK(KABIJ),1)
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c        xt2cme = dnorm2(lent2,WORK,1)
c     else
c        write(LUPRI,*) '   H1TERM: insuf. core for squaring L(iajb)'
c        write(LUPRI,*) '   H1TERM: insuf. core for 2CME amplitudes'
c        write(LUPRI,*) '   H1TERM: LWORK = ',LWORK
c        xsiajb = 1.0D9
c        xt2cme = 1.0D9
c     endif
c     write(LUPRI,*)
c     write(LUPRI,*) '   Entering H1TERM:'
c     write(LUPRI,*) '   ================'
c     write(LUPRI,*) '   LWORK = ',LWORK
c     write(LUPRI,*) '   Norm of K(kl,i;c) : ',xknorm
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xniajb,' (packed)'
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xsiajb,' (squared)'
c     write(LUPRI,*) '   Norm of T1AM      : ',t1norm
c     write(LUPRI,*) '   Norm of T2VO      : ',t2norm
c     write(LUPRI,*) '   Norm of T2CME     : ',xt2cme
c     write(LUPRI,*)
C
C--------------------------------------------------
C     Symmetry assignment: T1AM is total symmetric!
C--------------------------------------------------
C
      ISYMBJ = 1
C
C----------------------------
C     Loop over a symmetries.
C----------------------------
C
      DO ISYMA = 1,NSYM
C
         ISYMI = ISYMA
C
         IF ((NVIR(ISYMA).LE.0) .OR. (NRHF(ISYMI).LE.0)) GOTO 1000
C
C-------------------
C        Allocation.
C-------------------
C
         KCHOO = 1
         KEND1 = KCHOO + NMAIJA(ISYMI)
         LWRK1 = LWORK - KEND1 + 1
C
         IF (LWRK1 .LE. 0) THEN
            WRITE(LUPRI,*) 'Insufficient memory in CCHO_H1TERM'
            WRITE(LUPRI,*) 'Need (more than): ',KEND1-1
            WRITE(LUPRI,*) 'Available       : ',LWORK
            CALL QUIT(' Insufficient memory in CCHO_H1TERM ')
         ENDIF
C
         LWRKA = LWRK1/2
C
         IF (LWRKA .LE. 0) THEN
            WRITE(LUPRI,*) 'Insufficient memory in CCHO_H1TERM'
            WRITE(LUPRI,*) 'Memory available for batching: ',LWRKA
            WRITE(LUPRI,*) 'Total memory available       : ',LWORK
            CALL QUIT(' Insufficient memory in CCHO_H1TERM ')
         ENDIF
C
C----------------------------
C        Set up batch over a.
C----------------------------
C
         MINMEM = NMAIJA(ISYMA) + NT1AM(ISYMBJ)
         IF (FBATCH) THEN
            NEFA  = MIN(NONA,NVIR(ISYMA))
            LEFF  = NEFA*MINMEM + 1
            LWRKA = MIN(LWRKA,LEFF)
         ENDIF
         NUMA = MIN(LWRKA/MINMEM,NVIR(ISYMA))
C
         IF (NUMA .LE. 0) THEN
            WRITE(LUPRI,*)
     &      'Insufficient memory for batch in CCHO_H1TERM'
            WRITE(LUPRI,*) 'NUMA = ',NUMA,' (ISYMA = ',ISYMA,')'
            WRITE(LUPRI,*) 'Memory available for batch: ',LWRKA
            WRITE(LUPRI,*) 'Minimum memory required   : ',MINMEM
            CALL QUIT(' Insufficient memory in CCHO_H1TERM ')
         ENDIF
C
         NBATA = (NVIR(ISYMA) - 1)/NUMA + 1
C
         IF (PRINT) THEN
            WRITE(LUPRI,'(3X,A,I1,A,/,3X,A)')
     &      'Batch over A, symmetry ',ISYMA,':',
     &      '-------------------------'
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10)')
     &      'Minimum work space required   : ',MINMEM,
     &      'Work space available for batch: ',LWRKA
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10,/)')
     &      'Number of virtual orbitals    : ',NVIR(ISYMA),
     &      'Required number of A-batches  : ',NBATA
         ENDIF
C
         DO IBATA = 1,NBATA
C
            NUMIA = NUMA
            IF (IBATA .EQ. NBATA) THEN
               NUMIA = NVIR(ISYMA) - NUMA*(NBATA - 1)
            ENDIF
C
            IA1 = NUMA*(IBATA - 1) + 1
C
            IF (PRINT) THEN
               WRITE(LUPRI,'(6X,A,I10,A,/,6X,A)')
     &         'A-batch number ',IBATA,':',
     &         '--------------------------'
               WRITE(LUPRI,'(6X,A,I10,1X,I10,/)')
     &         'First and last A: ',IA1,IA1+NUMIA-1
            ENDIF
C
C----------------------
C           Allocation.
C----------------------
C
            KCHOV = KEND1
            KMINT = KCHOV + NT1AM(ISYMBJ)*NUMIA
            KEND2 = KMINT + NMAIJA(ISYMA)*NUMIA
            LWRK2 = LWORK - KEND2 + 1
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Batching bug in CCHO_H1TERM (1)'
               CALL QUIT(' Batching error in CCHO_H1TERM ')
            ENDIF
C
C-----------------------------------------------
C           Extract integrals:
C           M(kl,c,#a) = L(kcla) = XIAJB(ck,al).
C-----------------------------------------------
C
            CALL CCHO_H1XM(XIAJB,WORK(KMINT),IA1,NUMIA,ISYMA)

c     lenm   = NMAIJA(ISYMA)*NUMIA
c     xmnorm = dnorm2(lenm,WORK(KMINT),1)
c     write(LUPRI,*) '   H1TERM: ISYMA,IA1,NUMIA,IBATA,NBATA: ',
c    &                       ISYMA,IA1,NUMIA,IBATA,NBATA
c     write(LUPRI,*) '   H1TERM: Norm of M(klc;#a): ',xmnorm
C
C-----------------------------------------------------
C           Set up batch loop over i (symmetry=ISYMA).
C-----------------------------------------------------
C
            MINMEM = NMAIJA(ISYMI) + NT1AM(ISYMBJ)*NUMIA + 2*NUMIA
            IF (FBATCH) THEN
               NEFI  = MIN(NONI,NRHF(ISYMI))
               LEFF  = NEFI*MINMEM + 1
               LWRK2 = MIN(LWRK2,LEFF)
            ENDIF
            NUMI = MIN(LWRK2/MINMEM,NRHF(ISYMI))
C
            IF (NUMI .LE. 0) THEN
               WRITE(LUPRI,*)
     &         'Insufficient memory for batch in CCHO_H1TERM'
               WRITE(LUPRI,*) 'NUMI = ',NUMI,' (ISYMI = ',ISYMI,
     &                    ', ISYMA = ',ISYMA,')'
               WRITE(LUPRI,*) 'Minimum memory required   : ',MINMEM
               WRITE(LUPRI,*) 'Memory available for batch: ',LWRK2
               WRITE(LUPRI,*) 'Total memory available    : ',LWORK
               CALL QUIT(' Insufficient memory in CCHO_H1TERM ')
            ENDIF
C
            NBATI = (NRHF(ISYMI) - 1)/NUMI + 1
C
            IF (PRINT) THEN
               WRITE(LUPRI,'(9X,A,I1,A,/,9X,A)')
     &         'Batch over I, symmetry ',ISYMI,':',
     &         '-------------------------'
               WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10)')
     &         'Minimum work space required   : ',MINMEM,
     &         'Work space available for batch: ',LWRK2
               WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10,/)')
     &         'Number of occupied orbitals   : ',NRHF(ISYMI),
     &         'Required number of I-batches  : ',NBATI
            ENDIF
C
            DO IBATI = 1,NBATI
C
               NUMII = NUMI
               IF (IBATI .EQ. NBATI) THEN
                  NUMII = NRHF(ISYMI) - NUMI*(NBATI - 1)
               ENDIF
C
               II1 = NUMI*(IBATI - 1) + 1
C
               IF (PRINT) THEN
                  WRITE(LUPRI,'(12X,A,I10,A,/,12X,A)')
     &            'I-batch number ',IBATI,':',
     &            '--------------------------'
                  WRITE(LUPRI,'(12X,A,I10,1X,I10,/)')
     &            'First and last I: ',II1,II1+NUMII-1
               ENDIF
C
C-------------------------
C              Allocation.
C-------------------------
C
               KKINT = KEND2
               KT2AM = KKINT + NMAIJA(ISYMI)*NUMII
               KZMAT = KT2AM + NT1AM(ISYMBJ)*NUMIA*NUMII
               KUMAT = KZMAT + NUMIA*NUMII
               KEND3 = KUMAT + NUMIA*NUMII
               LWRK3 = LWORK - KEND3 + 1
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Batching bug in CCHO_H1TERM (2)'
                  CALL QUIT(' Batching error in CCHO_H1TERM ')
               ENDIF
C
C---------------------------------------------------------
C              Extract integrals:
C              K(kl,c,#i) = (ck|il) = (kc|li) = K(kl,i,c).
C---------------------------------------------------------
C
               CALL CCHO_H1XK(XKINT,WORK(KKINT),II1,NUMII,ISYMI)

c     lenkk  = NMAIJA(ISYMI)*NUMII
c     xkknrm = dnorm2(lenkk,WORK(KKINT),1)
c     write(LUPRI,*) '   H1TERM: ISYMI,II1,NUMII,IBATI,NBATI: ',
c    &                       ISYMI,II1,NUMII,IBATI,NBATI
c     write(LUPRI,*) '   H1TERM: Norm of K(klc,#i): ',xkknrm
C
C------------------------------------------------
C              Extract amplitudes:
C              T(bj,#a#i) = 2*T(bj,ai) - T(bi,aj)
C------------------------------------------------
C
               CALL CCHO_H1XT(T2VO,WORK(KT2AM),IA1,NUMIA,ISYMA,
     &                        II1,NUMII,ISYMI)

c     lentt  = NT1AM(ISYMBJ)*NUMIA*NUMII
c     xttnrm = dnorm2(lentt,WORK(KT2AM),1)
c     write(LUPRI,*) '   H1TERM: Norm of t(bj,#a#i): ',xttnrm
C
C----------------------------------
C              Start Cholesky loop.
C----------------------------------
C
               DO ICHO = 1,NUMCHO
C
C--------------------------------------------------------
C                 Get Cholesky vectors:
C                 d(kl,c), fixed sym. klc = ISYMI, update
C                 d(bj,#a), fixed sym. bj = 1, update.
C--------------------------------------------------------
C
                  CALL CCHO_DECHO6(FOCKD,CHOELE,NUMCHO,ICHO,
     &                             WORK(KCHOO),WORK(KCHOV),ISYMI,
     &                             ISYMBJ,IA1,NUMIA,ISYMA)
C
C----------------------------------------------------
C                 Scale:
C                 T(bj,#a#i) <- d(bj,#a) * T(bj,#a#i)
C                 K(klc,#i)  <- d(klc) * K(klc,#i)
C----------------------------------------------------
C
                  CALL CCHO_SCVEC(WORK(KT2AM),WORK(KCHOV),
     &                            NT1AM(ISYMBJ)*NUMIA,NUMII)
                  CALL CCHO_SCVEC(WORK(KKINT),WORK(KCHOO),NMAIJA(ISYMI),
     &                            NUMII)
C
C-----------------------------------
C                 Calculate Z and U.
C-----------------------------------
C
                  NAI   = NUMIA*NUMII
                  NBJ   = NT1AM(ISYMBJ)
                  NTOBJ = MAX(NBJ,1)
C
                  CALL DGEMV('T',NBJ,NAI,ONE,WORK(KT2AM),NTOBJ,
     &                       T1AM,1,ZERO,WORK(KZMAT),1)
C
                  NKLC  = NMAIJA(ISYMI)
                  NTKLC = MAX(NKLC,1)
C
                  CALL DGEMM('T','N',NUMIA,NUMII,NKLC,
     &                       ONE,WORK(KMINT),NTKLC,WORK(KKINT),NTKLC,
     &                       ZERO,WORK(KUMAT),NUMIA)
C
C-----------------------------------------------
C                 Calculate energy contribution.
C-----------------------------------------------
C
                  H1COR = DDOT(NAI,WORK(KZMAT),1,WORK(KUMAT),1)
C
C-----------------------------------
C                 Check convergence.
C-----------------------------------
C
                  EH1 = EH1 + H1COR
C 
              ENERGH(ICHO) = ENERGH(ICHO) + H1COR
C 
C Decommented by Domenico
                 TNOW = SECOND()
                 DELTAT = TNOW - TLAST
                 TLAST = TNOW
                 SCNDSH(ICHO) = SCNDSH(ICHO) 
     &                         + DELTAT
C Decommented by Domenico
C
C
c                 IF (PRINT) THEN
c                    WRITE(LUPRI,'(15X,A,I3,A,/,15X,A)')
c    &               'Status after Cholesky vector',ICHO,':',
c    &               '--------------------------------'
c                    IF (ABS(H1COR) .LT. THRCHO) THEN
c                       WRITE(LUPRI,'(15X,A)') 'H1 term converged'
c                    ELSE
c                       WRITE(LUPRI,'(15X,A)') 'H1 term not converged'
c                    ENDIF
c                    TIM = SECOND() - TIMT
c                    WRITE(LUPRI,'(15X,A,F10.2,A,/)')
c    &               'Accumulated H1-time: ',TIM,' seconds'
c                 ENDIF
C
                  IF (DABS(H1COR) .LT. THRCHO) GOTO 999
C
               ENDDO
C
  999          CONTINUE
C
            ENDDO
C
         ENDDO
C
 1000    CONTINUE
C
      ENDDO

c     lenk  = NTRAOC(1)
c     liajb = NT2AM(1)
c     lent2 = NT2SQ(1)
c     lent1 = NT1AM(1)
c     xknorm = dnorm2(lenk,XKINT,1)
c     xniajb = dnorm2(liajb,XIAJB,1)
c     t2norm = dnorm2(lent2,T2VO,1)
c     t1norm = dnorm2(lent1,T1AM,1)
c     if (lwork .ge. lent2) then
c        call cc_t2sq(XIAJB,WORK,1)
c        xsiajb = dnorm2(lent2,WORK,1)
c        call dzero(WORK,lent2)
c        call daxpy(lent2,2.0D0,T2VO,1,WORK,1)
c        DO ISYMJI = 1,NSYM
c           ISYMAB = ISYMJI
c           DO ISYMI = 1,NSYM
c              ISYMJ = MULD2H(ISYMI,ISYMJI)
c              DO I = 1,NRHF(ISYMI)
c                 DO J = 1,NRHF(ISYMJ)
c                    JI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I - 1) + J
c                    IJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
c                    KABIJ = IT2VO(ISYMAB,ISYMJI)
c    &                     + NMATAB(ISYMAB)*(IJ - 1) + 1
c                    KABJI = IT2VO(ISYMAB,ISYMJI)
c    &                     + NMATAB(ISYMAB)*(JI - 1) + 1
c                    CALL DAXPY(NMATAB(ISYMAB),-1.0D0,T2VO(KABJI),1,
c    &                                                WORK(KABIJ),1)
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c        xt2cme = dnorm2(lent2,WORK,1)
c     else
c        write(LUPRI,*) '   H1TERM: insuf. core for squaring L(iajb)'
c        write(LUPRI,*) '   H1TERM: insuf. core for 2CME amplitudes'
c        write(LUPRI,*) '   H1TERM: LWORK = ',LWORK
c        xsiajb = 1.0D9
c        xt2cme = 1.0D9
c     endif
c     write(LUPRI,*)
c     write(LUPRI,*) '   Exiting H1TERM:'
c     write(LUPRI,*) '   ==============='
c     write(LUPRI,*) '   LWORK = ',LWORK
c     write(LUPRI,*) '   Norm of K(kl,i;c) : ',xknorm
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xniajb,' (packed)'
c     write(LUPRI,*) '   Norm of L(ia,jb)  : ',xsiajb,' (squared)'
c     write(LUPRI,*) '   Norm of T1AM      : ',t1norm
c     write(LUPRI,*) '   Norm of T2VO      : ',t2norm
c     write(LUPRI,*) '   Norm of T2CME     : ',xt2cme
c     write(LUPRI,*)
C
      RETURN
      END
C  /* Deck ccho_h1xt */
      SUBROUTINE CCHO_H1XT(T2VO,TAM,IA1,NUMIA,ISYMA,II1,NUMII,ISYMI)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     TAM(bj,#a#i) = 2*T(bj,ai) - T(bi,aj)
C
C                  = 2*T2VO(ba,ji) - T2VO(ba,ij)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2VO(*),TAM(*)
      INTEGER BAIJ,BAJI,AI,BJAI
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (XMONE = -1.00D0, TWO = 2.00D0)
C
      ISYMBJ = MULD2H(ISYMA,ISYMI)
      LENTOT = NT1AM(ISYMBJ)*NUMIA*NUMII
      CALL DZERO(TAM,LENTOT)
C
      DO ISYMIJ = 1,NSYM
C
         ISYMBA = ISYMIJ
         ISYMJ  = MULD2H(ISYMIJ,ISYMI)
         ISYMB  = MULD2H(ISYMBA,ISYMA)
C
         IF (NVIR(ISYMB) .LE. 0) GOTO 100
C
         DO J = 1,NRHF(ISYMJ)
            DO I = 1,NUMII
C
               II = II1 + I - 1
C
               JI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(II - 1) + J
               IJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + II
C
               DO A = 1,NUMIA
C
                  IA = IA1 + A - 1
C
                  AI = NUMIA*(I - 1) + A
C
                  BAJI = IT2VO(ISYMBA,ISYMIJ) + NMATAB(ISYMBA)*(JI - 1)
     &                 + IMATAB(ISYMB,ISYMA) + NVIR(ISYMB)*(IA - 1) + 1
                  BAIJ = IT2VO(ISYMBA,ISYMIJ) + NMATAB(ISYMBA)*(IJ - 1)
     &                 + IMATAB(ISYMB,ISYMA) + NVIR(ISYMB)*(IA - 1) + 1
                  BJAI = NT1AM(ISYMBJ)*(AI - 1)
     &                 + IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + 1
C
                  CALL DAXPY(NVIR(ISYMB),TWO,T2VO(BAJI),1,TAM(BJAI),1)
                  CALL DAXPY(NVIR(ISYMB),XMONE,T2VO(BAIJ),1,TAM(BJAI),1)
C
               ENDDO
C
            ENDDO
         ENDDO
C
  100    CONTINUE
C
      ENDDO 
C
      RETURN
      END
C  /* Deck ccho_h1xm */
      SUBROUTINE CCHO_H1XM(XIAJB,XM,IA1,NUMIA,ISYMA)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     M(kl,c,#a) = XIAJB(ck,al).
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XIAJB(*),XM(*)
      INTEGER AL,CK,CKAL
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      ISYKLC = ISYMA
      NKLC   = NMAIJA(ISYKLC)
C
      DO ISYMAL = 1,NSYM
C
         ISYML  = MULD2H(ISYMAL,ISYMA)
         ISYMCK = ISYMAL
C
         DO L = 1,NRHF(ISYML)
            DO A = 1,NUMIA
C
               IA = IA1 + A - 1
C
               AL = IT1AM(ISYMA,ISYML) + NVIR(ISYMA)*(L - 1) + IA
C
               DO ISYMK = 1,NSYM
C
                  ISYMC  = MULD2H(ISYMK,ISYMCK)
                  ISYMKL = MULD2H(ISYMK,ISYML)
C
                  DO K = 1,NRHF(ISYMK)
C
                     KL = IMATIJ(ISYMK,ISYML) + NRHF(ISYMK)*(L - 1) + K
C
                     DO C = 1,NVIR(ISYMC)
C
                        CK  = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1)
     &                      + C
                        KLC = IMAIJA(ISYMKL,ISYMC)
     &                      + NMATIJ(ISYMKL)*(C - 1) + KL
C
                        KLCA = NKLC*(A - 1) + KLC
                        CKAL = IT2AM(ISYMCK,ISYMAL) + INDEX(CK,AL)
C
                        XM(KLCA) = XIAJB(CKAL)
C
                     ENDDO
C
                  ENDDO
C
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_decho5 */
      SUBROUTINE CCHO_DECHO5(FOCKD,CHOELE,NUMCHO,ICHO,CHOL,OCCHO,VICHO,
     &                       ISYMMN,IA1,NUMA,ISYMA,J1,NUMJ,ISYMJ)
C
C     JLC, BFR, TBP, HK, AS, October 2002.
C
C     Construct vector for Cholesky decomposition:
C
C     For the occupied part,  CHOL(mn,b)  actual vector
C     For the occupied part, OCCHO(ck,#j) update vector
C     For the virtual  part, VICHO(ck,#a) actual vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),CHOL(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C------------------------
C     Occupied part CHOL.
C------------------------
C
C ORIGINAL:
c     IND = 0
c     DO ISYMB = 1,NSYM
c        DO B = 1,NVIR(ISYMB)
c           DO ISYMN = 1,NSYM
c              ISYMM = MULD2H(ISYMN,ISYMMN)
c              DO N = 1,NRHF(ISYMN)
c                 DO M = 1,NRHF(ISYMM)
C
c                    KOFFB = IVIR(ISYMB) + B
c                    KOFFN = IRHF(ISYMN) + N
c                    KOFFM = IRHF(ISYMM) + M
C
c                    OME = FOCKD(KOFFB)-FOCKD(KOFFN)-FOCKD(KOFFM)
C
c                    IND = IND + 1
c                    CHOL(IND) = 
c    &                   SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
c                 END DO
c              END DO
c           END DO
c        END DO
c     ENDDO
c
C NEW:
      IND = 0
      DO ISYMB = 1,NSYM
         DO B = 1,NVIR(ISYMB)
            DO ISYMN = 1,NSYM
               ISYMM = MULD2H(ISYMN,ISYMMN)
               DO N = 1,NRHF(ISYMN)
                  DO M = 1,NRHF(ISYMM)
C
                     KOFFB = IVIR(ISYMB) + B
                     KOFFN = IRHF(ISYMN) + N
                     KOFFM = IRHF(ISYMM) + M
C
                     IND = IND + 1
                     CHOL(IND) = FOCKD(KOFFB)-FOCKD(KOFFN)-FOCKD(KOFFM)
C
                  END DO
               END DO
            END DO
         END DO
      ENDDO
C
C     Construct the vector.
C     ---------------------
C
      DO P = 1,IND
C
         OMEGA = CHOL(P)
C
         CHOL(P) = SQRT(TWO*CHOELE(ICHO))/(OMEGA+CHOELE(ICHO))
C
         DO JCHO = 1,ICHO-1
            CHOL(P) = CHOL(P) 
     &               *(OMEGA-CHOELE(JCHO))/(OMEGA+CHOELE(JCHO))
         END DO
C
      END DO
C
C-------------------------
C     Occupied part OCCHO. 
C-------------------------
C
      IF (ICHO .EQ. 1) THEN
C
         IND = 0
         DO J = J1,J1+NUMJ-1 
            DO ISYMCK = 1,NSYM
               DO ISYMK = 1,NSYM
                  ISYMC = MULD2H(ISYMK,ISYMCK)
                  DO K = 1,NRHF(ISYMK)
                     DO C = 1,NVIR(ISYMC)
C
                        KOFFJ = IRHF(ISYMJ) + J
                        KOFFK = IRHF(ISYMK) + K
                        KOFFC = IVIR(ISYMC) + C
C
                        OME = FOCKD(KOFFC)-FOCKD(KOFFK)-FOCKD(KOFFJ)
C
                        IND = IND + 1
                        OCCHO(IND) = 
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
      ELSE
C
         IND = 0
         DO J = J1,J1+NUMJ-1 
            DO ISYMCK = 1,NSYM
               DO ISYMK = 1,NSYM
                  ISYMC = MULD2H(ISYMK,ISYMCK)
                  DO K = 1,NRHF(ISYMK)
                     DO C = 1,NVIR(ISYMC)
C
                        KOFFJ = IRHF(ISYMJ) + J
                        KOFFK = IRHF(ISYMK) + K
                        KOFFC = IVIR(ISYMC) + C
C
                        OME = FOCKD(KOFFC)-FOCKD(KOFFK)-FOCKD(KOFFJ)
C
                        IND = IND + 1
                        OCCHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(IND,FACTOR,OCCHO,1)
C
      ENDIF
C
C------------------------
C     Virtual part VICHO.
C------------------------
C
C ORIGINAL:
c     IND = 0
c     DO A = IA1,IA1+NUMA-1 
c        DO ISYMCK = 1,NSYM
c           DO ISYMK = 1,NSYM
c              ISYMC = MULD2H(ISYMK,ISYMCK)
c              DO K = 1,NRHF(ISYMK)
c                 DO C = 1,NVIR(ISYMC)
C
c                    KOFFA = IVIR(ISYMA) + A
c                    KOFFK = IRHF(ISYMK) + K
c                    KOFFC = IVIR(ISYMC) + C
C
c                    OME = FOCKD(KOFFA)-FOCKD(KOFFK)+FOCKD(KOFFC)
C
c                    IND = IND + 1
c                    VICHO(IND) = 
c    &                   SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
c                 END DO
c              END DO
c           END DO
c        ENDDO
c     END DO
c
C NEW:
      IND = 0
      DO A = 1,NUMA
         KOFFA = IVIR(ISYMA) + IA1 + A - 1
         DO ISYMCK = 1,NSYM
            DO ISYMK = 1,NSYM
               ISYMC = MULD2H(ISYMK,ISYMCK)
               DO K = 1,NRHF(ISYMK)
                  KOFFK = IRHF(ISYMK) + K
                  DO C = 1,NVIR(ISYMC)
                     KOFFC = IVIR(ISYMC) + C
C
                     IND = IND + 1
                     VICHO(IND) = FOCKD(KOFFA) + FOCKD(KOFFC)
     &                          - FOCKD(KOFFK)
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C     Construct the vector.
C     ---------------------
C
      DO P = 1,IND
C
         OMEGA = VICHO(P)
C
C
         VICHO(P) = SQRT(TWO*CHOELE(ICHO))/(OMEGA+CHOELE(ICHO))
C
         DO JCHO = 1,ICHO-1
            VICHO(P) = VICHO(P) 
     &               *(OMEGA-CHOELE(JCHO))/(OMEGA+CHOELE(JCHO))
         END DO
C
      END DO
C
      RETURN
      END
C  /* Deck ccho_decho7 */
      SUBROUTINE CCHO_DECHO7(FOCKD,CHOELE,NUMCHO,ICHO,CHOL,OCCHO,VICHO,
     &                       ISYDJA,IB1,NUMB,ISYMB,I1,NUMI,ISYMI)
C
C     JLC, BFR, TBP, HK, AS, October 2002.
C
C     Construct vector for Cholesky decomposition:
C
C     For the virtual  part,  CHOL(dj,a)  actual vector
C     For the occupied part, OCCHO(ck,#i) update vector
C     For the virtual  part, VICHO(ck,#b) actual vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),CHOL(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'


         CALL QUIT('Too many vectors required')
      END IF
C
C-----------------------
C     Virtual part CHOL.
C-----------------------
C
      IND = 0
      DO ISYMA = 1,NSYM
         ISYMDJ = MULD2H(ISYMA,ISYDJA)
         DO A = 1,NVIR(ISYMA)
            DO ISYMJ = 1,NSYM
               ISYMD = MULD2H(ISYMJ,ISYMDJ)
               DO J = 1,NRHF(ISYMJ)
                  DO D = 1,NVIR(ISYMD)
C
                     KOFFA = IVIR(ISYMA) + A
                     KOFFJ = IRHF(ISYMJ) + J
                     KOFFD = IVIR(ISYMD) + D
C
                     IND = IND + 1
                     CHOL(IND) = FOCKD(KOFFA)-FOCKD(KOFFJ)+FOCKD(KOFFD)
C
                  END DO
               END DO
            END DO
         END DO
      ENDDO
C
C     Construct the vector.
C     ---------------------
C
      DO P = 1,IND
C
         OMEGA = CHOL(P)
C
         CHOL(P) = SQRT(TWO*CHOELE(ICHO))/(OMEGA+CHOELE(ICHO))
C
         DO JCHO = 1,ICHO-1
            CHOL(P) = CHOL(P) 
     &               *(OMEGA-CHOELE(JCHO))/(OMEGA+CHOELE(JCHO))
         END DO
C
      END DO
C
C-------------------------
C     Occupied part OCCHO. 
C-------------------------
C
      IF (ICHO .EQ. 1) THEN
C
         IND = 0
         DO I = I1,I1+NUMI-1 
            DO ISYMCK = 1,NSYM
               DO ISYMK = 1,NSYM
                  ISYMC = MULD2H(ISYMK,ISYMCK)
                  DO K = 1,NRHF(ISYMK)
                     DO C = 1,NVIR(ISYMC)
C
                        KOFFI = IRHF(ISYMI) + I
                        KOFFK = IRHF(ISYMK) + K
                        KOFFC = IVIR(ISYMC) + C
C
                        OME = FOCKD(KOFFC)-FOCKD(KOFFK)-FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) = 
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
      ELSE
C
         IND = 0
         DO I = I1,I1+NUMI-1 
            DO ISYMCK = 1,NSYM
               DO ISYMK = 1,NSYM
                  ISYMC = MULD2H(ISYMK,ISYMCK)
                  DO K = 1,NRHF(ISYMK)
                     DO C = 1,NVIR(ISYMC)
C
                        KOFFI = IRHF(ISYMI) + I
                        KOFFK = IRHF(ISYMK) + K
                        KOFFC = IVIR(ISYMC) + C
C
                        OME = FOCKD(KOFFC)-FOCKD(KOFFK)-FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(IND,FACTOR,OCCHO,1)
C
      ENDIF
C
C------------------------
C     Virtual part VICHO.
C------------------------
C
      IND = 0
      DO B = IB1,IB1+NUMB-1 
         DO ISYMCK = 1,NSYM
            DO ISYMK = 1,NSYM
               ISYMC = MULD2H(ISYMK,ISYMCK)
               DO K = 1,NRHF(ISYMK)
                  DO C = 1,NVIR(ISYMC)
C
                     KOFFB = IVIR(ISYMB) + B
                     KOFFK = IRHF(ISYMK) + K
                     KOFFC = IVIR(ISYMC) + C
C
                     IND = IND + 1
                     VICHO(IND) = FOCKD(KOFFB)-FOCKD(KOFFK)+FOCKD(KOFFC)
C
                  END DO
               END DO
            END DO
         ENDDO
      END DO
C
C     Construct the vector.
C     ---------------------
C
      DO P = 1,IND
C
         OMEGA = VICHO(P)
C
C
         VICHO(P) = SQRT(TWO*CHOELE(ICHO))/(OMEGA+CHOELE(ICHO))
C
         DO JCHO = 1,ICHO-1
            VICHO(P) = VICHO(P) 
     &               *(OMEGA-CHOELE(JCHO))/(OMEGA+CHOELE(JCHO))
         END DO
C
      END DO
C
      RETURN
      END
C  /* Deck ccho_decho6 */
      SUBROUTINE CCHO_DECHO6(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                       ISYIKC,ISYMDJ,IB1,NUMIB,ISYMB)
C
C     JLC, BFR, TBP, HK, AS, October 2002.
C
C     Construct vector for Cholesky decomposition:
C
C     For the occupied part, OCCHO(ik,c) update vector,
C                                        all ikc of sym. ISYIKC
C     For the virtual  part, VICHO(dj,#b) actual vector,
C                                         fixed dj sym. ISYMDJ
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION FOCKD(*), CHOELE(*), OCCHO(*), VICHO(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C-------------------------
C     Occupied part OCCHO.
C-------------------------
C
      IF (ICHO .EQ. 1) THEN
C
         IND = 0
         DO ISYMC = 1,NSYM
            ISYMIK = MULD2H(ISYMC,ISYIKC)
            DO C = 1,NVIR(ISYMC)
               KOFFC = IVIR(ISYMC) + C
               DO ISYMK = 1,NSYM
                  ISYMI = MULD2H(ISYMK,ISYMIK)
                  DO K = 1,NRHF(ISYMK)
                     KOFFK = IRHF(ISYMK) + K
                     DO I = 1,NRHF(ISYMI)
                        KOFFI = IRHF(ISYMI) + I
C
                        OME = FOCKD(KOFFC) - FOCKD(KOFFK) - FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) = SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
      ELSE
C
         IND = 0
         DO ISYMC = 1,NSYM
            ISYMIK = MULD2H(ISYMC,ISYIKC)
            DO C = 1,NVIR(ISYMC)
               KOFFC = IVIR(ISYMC) + C
               DO ISYMK = 1,NSYM
                  ISYMI = MULD2H(ISYMK,ISYMIK)
                  DO K = 1,NRHF(ISYMK)
                     KOFFK = IRHF(ISYMK) + K
                     DO I = 1,NRHF(ISYMI)
                        KOFFI = IRHF(ISYMI) + I
C
                        OME = FOCKD(KOFFC) - FOCKD(KOFFK) - FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/(OME+CHOELE(ICHO))
C
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(IND,FACTOR,OCCHO,1)
C
      ENDIF
C
C------------------------
C     Virtual part VICHO.
C------------------------
C
      IF (ICHO .EQ. 1) THEN
C
         IND = 0
         DO B = IB1,IB1+NUMIB-1
            KOFFB = IVIR(ISYMB) + B
            DO ISYMJ = 1,NSYM
               ISYMD = MULD2H(ISYMJ,ISYMDJ)
               DO J = 1,NRHF(ISYMJ)
                  KOFFJ = IRHF(ISYMJ) + J
                  DO D = 1,NVIR(ISYMD)
                     KOFFD = IVIR(ISYMD) + D
C
                     OME = FOCKD(KOFFB) - FOCKD(KOFFJ) + FOCKD(KOFFD)
C
                     IND = IND + 1
                     VICHO(IND) = SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
      ELSE
C
         IND = 0
         DO B = IB1,IB1+NUMIB-1
            KOFFB = IVIR(ISYMB) + B
            DO ISYMJ = 1,NSYM
               ISYMD = MULD2H(ISYMJ,ISYMDJ)
               DO J = 1,NRHF(ISYMJ)
                  KOFFJ = IRHF(ISYMJ) + J
                  DO D = 1,NVIR(ISYMD)
                     KOFFD = IVIR(ISYMD) + D
C
                     OME = FOCKD(KOFFB) - FOCKD(KOFFJ) + FOCKD(KOFFD)
C
                     IND = IND + 1
                     VICHO(IND) =
     &                    (OME-CHOELE(ICHO-1))/(OME+CHOELE(ICHO))
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(IND,FACTOR,VICHO,1)
C
      ENDIF
C
      RETURN
      END
C  /* Deck ccho_h1xk */
      SUBROUTINE CCHO_H1XK(XKINT,XKSUB,II1,NUMII,ISYMI)
C
C     JLC, BFR, TBP, HK, and AS. October, 2002.
C
C     Extract occupied integrals for H1 term:
C
C        XKSUB(klc,#i) = (ck|#il) = (kc|l#i)
C                      = XKINT(kl,#i,c)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XKINT(*), XKSUB(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
c     lenk   = NTRAOC(1)
c     xknorm = dnorm2(lenk,XKINT,1)
c     write(LUPRI,*) '      H1XK: Entry; norm of K(kl,i;c): ',xknorm
C
      ISYKLC = ISYMI
C
      DO ISYMC = 1,NSYM
C
         ISYMKL = MULD2H(ISYMC,ISYKLC)
         ISYKLI = MULD2H(ISYMKL,ISYMI)
C
         DO C = 1,NVIR(ISYMC)
            DO I = 1,NUMII
C
               II = II1 + I - 1
C
               KLIC = ISJIKA(ISYKLI,ISYMC)
     &              + NMAJIK(ISYKLI)*(C - 1)
     &              + ISJIK(ISYMKL,ISYMI)
     &              + NMATIJ(ISYMKL)*(II - 1) + 1
               KLCI = NMAIJA(ISYKLC)*(I - 1)
     &              + IMAIJA(ISYMKL,ISYMC)
     &              + NMATIJ(ISYMKL)*(C - 1) + 1
C
               CALL DCOPY(NMATIJ(ISYMKL),XKINT(KLIC),1,XKSUB(KLCI),1)
C
            ENDDO
         ENDDO
C
      ENDDO

c     lenk   = NTRAOC(1)
c     lens   = NMAIJA(ISYMI)*NUMII
c     xknorm = dnorm2(lenk,XKINT,1)
c     xsnorm = dnorm2(lens,XKSUB,1)
c     write(LUPRI,*) '      H1XK: Exit ; norm of K(kl,i;c): ',xknorm
c     write(LUPRI,*) '      H1XK: Exit ; norm of K(klc;#i): ',xsnorm
C
      RETURN
      END
C  /* Deck ccho_decho8 */
      SUBROUTINE CCHO_DECHO8(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                       ISYDLJ,ISYCDL)
C
C     JLC, BFR, TBP, HK, AS, June 2002.
C
C     Construct vector for Cholesky decomposition:
C
C     For the occupied part, OCCHO(dlj),  update vector
C     For the virtual  part, VICHO(cdl),  update vector
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION FOCKD(*), CHOELE(*), OCCHO(*), VICHO(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (TWO = 2.00D0)
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C--------------------------
C     Occupied part, OCCHO.
C--------------------------
C
      CALL CCHO_DECHO3(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                 ISYDLJ,ISYCDL,1,0)

C DEBUG:
c     DO ISYMC = 1,NSYM
c        ISYMDL = MULD2H(ISYMC,ISYCDL)
c        DO C = 1,NVIR(ISYMC)
c           DO ISYML = 1,NSYM
c              ISYMD = MULD2H(ISYML,ISYMDL)
c              DO L = 1,NRHF(ISYML)
c                 DO D = 1,NVIR(ISYMD)
c                    KDLC = ICKATR(ISYMDL,ISYMC)
c    &                    + NT1AM(ISYMDL)*(C - 1)
c    &                    + IT1AM(ISYMD,ISYML)
c    &                    + NVIR(ISYMD)*(L - 1) + D
c                    write(LUPRI,*) 'sym. d,l,c; d,l,c; VICHO:',
c    &                          ISYMD,ISYML,ISYMC,';',D,L,C,';',
c    &                          VICHO(KDLC)
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
C
C-------------------------
C     Virtual part, VICHO.
C-------------------------
C
      IF (ICHO .EQ. 1) THEN
C
         IND = 0
         DO ISYML = 1,NSYM
            ISYMCD = MULD2H(ISYML,ISYCDL)
            DO L = 1,NRHF(ISYML)
               KOFFL = IRHF(ISYML) + L
               DO ISYMD = 1,NSYM
                  ISYMC = MULD2H(ISYMD,ISYMCD)
                  DO D = 1,NVIR(ISYMD)
                     KOFFD = IVIR(ISYMD) + D
                     DO C = 1,NVIR(ISYMC)
                        KOFFC = IVIR(ISYMC) + C
C
                           OME = FOCKD(KOFFC) + FOCKD(KOFFD)
     &                         - FOCKD(KOFFL)
C
                           IND = IND + 1
                           VICHO(IND) =
     &                          SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         NDIMVI = IND
C
      ELSE
C
         IND = 0
         DO ISYML = 1,NSYM
            ISYMCD = MULD2H(ISYML,ISYCDL)
            DO L = 1,NRHF(ISYML)
               KOFFL = IRHF(ISYML) + L
               DO ISYMD = 1,NSYM
                  ISYMC = MULD2H(ISYMD,ISYMCD)
                  DO D = 1,NVIR(ISYMD)
                     KOFFD = IVIR(ISYMD) + D
                     DO C = 1,NVIR(ISYMC)
                        KOFFC = IVIR(ISYMC) + C
C
                           OME = FOCKD(KOFFC) + FOCKD(KOFFD)
     &                         - FOCKD(KOFFL)
C
                           IND = IND + 1
                           VICHO(IND) =
     &                          (OME-CHOELE(ICHO-1))/(OME+CHOELE(ICHO))
C
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(IND,FACTOR,VICHO,1)
C
         NDIMVI = IND
C
      ENDIF

C DEBUG:
c     write(LUPRI,*) '   DECHO8: Norm of VICHO(cd,l), sym. ',ISYCDL,': ',
c    &           dnorm2(NDIMVI,VICHO,1)
c     DO ISYML = 1,NSYM
c        ISYMCD = MULD2H(ISYML,ISYCDL)
c        DO L = 1,NRHF(ISYML)
c           DO ISYMD = 1,NSYM
c              ISYMC = MULD2H(ISYMD,ISYMCD)
c              DO D = 1,NVIR(ISYMD)
c                 DO C = 1,NVIR(ISYMC)
c                    KCDL = ICKASR(ISYMCD,ISYML)
c    &                    + NMATAB(ISYMCD)*(L - 1)
c    &                    + IMATAB(ISYMC,ISYMD)
c    &                    + NVIR(ISYMC)*(D - 1) + C
c                    write(LUPRI,*) 'sym. d,l,c; d,l,c; VICHO:',
c    &                          ISYMD,ISYML,ISYMC,';',D,L,C,';',
c    &                          VICHO(KCDL)
c                 ENDDO
c              ENDDO
c           ENDDO
c        ENDDO
c     ENDDO
C
      RETURN
      END
C  /* Deck ccho_fterm1 */
      SUBROUTINE CCHO_FTERM1(XOINT,T2VO,FOCKD,NUMCHO,CHOELE,WORK,LWORK,
     &                       E4F1,FBATCH,PRINT,NONJ,NONA)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     Calculate (part of F-terms):
C
C        E4F1 = Sum(j) Sum(kai) [ V(ka,i;j) * D(ia,k;j)
C                               + W(ka,i;j) * E(ia,k;j) ]
C
C     where
C
C        V(ka,i;j) = - Sum(em) d(em,j) * (em|ij) * s(em,ak)
C                    - Sum(em) d(em,j) * (ej|im) * s(ek,am)
C
C        W(ka,i;j) = 2 Sum(em) d(em,j) * (em|ij) * s(em,ak)
C                    - Sum(em) d(em,j) * (ej|im) * s(em,ak)
C
C        D(ia,k;j) =   Sum(dl) d(dl,a) * (dl|jk) * s(dl,ai)
C                    - Sum(dl) d(dl,a) * (jd|lk) * t(dl,ai)
C
C        E(ia,k;j) = - Sum(dl) d(dl,a) * (lk|jd) * t(di,al)
C
C    and s(ai,bj) = 2 t(ai,bj) - t(aj,bi), and d(em,j) and d(dl,a) denotes
C    the occupied and virtual parts of the Cholesky decomposition of the
C    orbital energy denominator.
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XOINT(*),T2VO(*),FOCKD(*),CHOELE(*)
      DIMENSION WORK(LWORK)
      LOGICAL   FBATCH,PRINT
C     
      DIMENSION IOFFINT(8),IOFFT(8),IOFFD(8)
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      TLAST = SECOND()
C
      IF (PRINT) THEN
         TIMT = SECOND()
         WRITE(LUPRI,'(A,/,A,/)')
     &   'Calculation of the F1 term:',
     &   '==========================='
      ENDIF
C
*     xitotal = zero
*     xjtotal = zero
*     xt21total = zero
*     xt22total = zero
*     xt23total = zero
*     xt24total = zero
*     xt25total = zero
C                 
      ICOUNT = 0
      DO ISYMDL = 1,NSYM
         IOFFD(ISYMDL) = ICOUNT
         ICOUNT = ICOUNT + NT1AM(ISYMDL)
      ENDDO        
      NTOTT1 = ICOUNT
C      
      DO ISYMJ = 1,NSYM
*        print *, 'ISYMJ = ',isymj
*        print *, '----------------------------------'
C     
         IF ( NRHF(ISYMJ) .EQ. 0 ) GOTO 1234
C      
         ISYMEMI = ISYMJ
C         
         LWORKJ = LWORK/2
C
C         
C-------------------------------
C        Batch over the J index.
C-------------------------------
C         
         LENJ = 4*NCKI(ISYMEMI) + NTOTT1
         IF (FBATCH) THEN
            NEFJ   = MIN(NONJ,NRHF(ISYMJ))
            LEFF   = NEFJ*LENJ + 1
            LWORKJ = MIN(LWORKJ,LEFF)
         ENDIF
         NUMJ = MIN(NRHF(ISYMJ),LWORKJ/LENJ)
         IF (NUMJ .EQ. 0) THEN
            WRITE(LUPRI,*) 'NUMJ .EQ. 0 IN CC_CHOPTF1'
            CALL QUIT('NOT ENOUGH SPACE IN CC_CHOPTF1 !!')
         ENDIF
C         
         NBATJ = (NRHF(ISYMJ)-1)/NUMJ + 1
C
         IF (PRINT) THEN
            WRITE(LUPRI,'(3X,A,I1,A,/,3X,A)')
     &      'Batch over J, symmetry ',ISYMJ,':',
     &      '-------------------------'
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10)')
     &      'Minimum work space required   : ',LENJ,
     &      'Work space available for batch: ',LWORKJ
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10,/)')
     &      'Number of occupied orbitals   : ',NRHF(ISYMJ),
     &      'Required number of J-batches  : ',NBATJ
         ENDIF
C        
         IJ2=0
         DO IBATJ = 1,NBATJ
*           print *, 'IBATJ = ',ibatj
*           xcount = xcount + one
            IJ1 = IJ2 + 1
            IJ2 = IJ2 + NUMJ
            IF (IJ2 .GT. NRHF(ISYMJ)) IJ2 = NRHF(ISYMJ)
            NUMIJ = IJ2 - IJ1 +1
C
            IF (PRINT) THEN
               WRITE(LUPRI,'(6X,A,I10,A,/,6X,A)')
     &         'J-batch number ',IBATJ,':',
     &         '--------------------------'
               WRITE(LUPRI,'(6X,A,I10,1X,I10,/)')
     &         'First and last J: ',IJ1,IJ2
            ENDIF
C        
            ICOUNT = 0
            DO ISYMI = 1,NSYM
               ISYMEM = MULD2H(ISYMI,ISYMJ)
               IOFFINT(ISYMI) = ICOUNT
               ICOUNT = ICOUNT + NT1AM(ISYMEM)*NRHF(ISYMI)*NUMIJ
            ENDDO
*           print *,'IOFFINT : ',ioffint
C            
C------------------------------
C           Dynamic allocation.
C------------------------------
C            
            KCHOO   = 1
            KIINT   = KCHOO  + NTOTT1*NUMIJ
            KJINT   = KIINT  + NCKI(ISYMEMI)*NUMIJ
            KITINT  = KJINT  + NCKI(ISYMEMI)*NUMIJ
            KJTINT  = KITINT + NCKI(ISYMEMI)*NUMIJ
            KENDJ   = KJTINT + NCKI(ISYMEMI)*NUMIJ
            LWRKJ   = LWORKJ - KENDJ + 1
            IF (LWRKJ .LT. 0) THEN
               CALL QUIT('Insufficient memory in CC_CHOPTF1')
            ENDIF
C             
C------------------------------------------------------------------
C           Extract integrals: I(em;i,j)=(im|ej), J(em;i,j)=(em|ij)
C           ordered for a j batch with a fixed ISYMJ
C           for each ISYMI all j values in the batch
C------------------------------------------------------------------
C             
*     print *,'xoint norm:',dnrm2(NTRAOC(1),xoint,1)            
            CALL CCHO_RDINTF1(XOINT,WORK(KIINT),WORK(KJINT),IJ1,
     &                          NUMIJ,ISYMJ,IOFFINT)
*     print *,'xoint norm after rdintf1:',dnrm2(NTRAOC(1),xoint,1)
*     l = ncki(isymemi)*numij
*     print *,'size xoint: ',ntraoc(1),'; size I,J : ',l
*     print *,'I(em,i#j) norm:',dnrm2(l,work(kiint),1)
*     xitotal = xitotal + dnrm2(l,work(kiint),1)**2
*     print *,'J(em,i#j) norm:',dnrm2(l,work(kjint),1)
*     xjtotal = xjtotal + dnrm2(l,work(kjint),1)**2
C        
            DO ISYMA = 1,NSYM
*        print *, 'ISYMA = ',isyma
C
               IF (NVIR(ISYMA) .EQ. 0) GOTO 1235
C              
               ISYMEMK = ISYMA
C                 
C---------------------------
C              Batch over A.
C---------------------------
C                 
               MXKAI = -1
               NKAI = 0
               DO ISYMI = 1,NSYM
                  ISYMKA = MULD2H(ISYMEMI,ISYMI)
                  ISYMK = MULD2H(ISYMKA,ISYMA)
                  NKAI = NRHF(ISYMI)*NRHF(ISYMK)*NUMIJ
                  IF (MXKAI .LT. NKAI) MXKAI = NKAI
               ENDDO
C                 
               LENA = 5*NCKI(ISYMEMK) + 4*MXKAI + NTOTT1
               NAVAIL = LWORK - KENDJ + 1
               IF (FBATCH) THEN
                  NEFA   = MIN(NVIR(ISYMA),NONA)
                  LEFF   = NEFA*LENA + 1
                  NAVAIL = MIN(NAVAIL,LEFF)
               ENDIF
               NUMA = MIN(NVIR(ISYMA),NAVAIL/LENA)
               IF (NUMA .EQ. 0) THEN
                  WRITE(LUPRI,*) 'NUMA .EQ. 0 in CC_CHOPTF1'
                  CALL QUIT('Not enough space in CC_CHOPTF1')
               ENDIF
C                 
               NBATA = (NVIR(ISYMA)-1)/NUMA + 1
C
               IF (PRINT) THEN
                  WRITE(LUPRI,'(9X,A,I1,A,/,9X,A)')
     &            'Batch over A, symmetry ',ISYMA,':',
     &            '-------------------------'
                  WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10)')
     &            'Minimum work space required   : ',LENA,
     &            'Work space available for batch: ',NAVAIL
                  WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10,/)')
     &            'Number of virtual orbitals    : ',NVIR(ISYMA),
     &            'Required number of A-batches  : ',NBATA
               ENDIF
C              
               IA2 = 0
               DO IBATA = 1,NBATA
*      print *, 'IBATA = ',ibata
C                 
                  IA1 = IA2 + 1
                  IA2 = IA2 + NUMA
                  IF (IA2 .GT. NVIR(ISYMA)) IA2 = NVIR(ISYMA)
                  NUMIA = IA2 - IA1 + 1
C
                  IF (PRINT) THEN
                     WRITE(LUPRI,'(12X,A,I10,A,/,12X,A)')
     &               'A-batch number ',IBATA,':',
     &               '--------------------------'
                     WRITE(LUPRI,'(12X,A,I10,1X,I10,/)')
     &               'First and last A: ',IA1,IA2
                  ENDIF
C              
                  ICOUNT = 0
                  DO ISYMK = 1,NSYM
                     ISYMEM = MULD2H(ISYMK,ISYMA)
                     IOFFT(ISYMK) = ICOUNT
                     ICOUNT = ICOUNT + NT1AM(ISYMEM)*NRHF(ISYMK)*NUMIA
                  ENDDO
C                    
C------------------------------------
C                 Dynamic allocation.
C------------------------------------
C                    
                  KCHOV =  KENDJ 
                  KT21  =  KCHOV + NTOTT1*NUMIA
                  KT22  =  KT21  + NCKI(ISYMA)*NUMIA
                  KT23  =  KT22  + NCKI(ISYMA)*NUMIA
                  KT24  =  KT23  + NCKI(ISYMA)*NUMIA
                  KT25  =  KT24  + NCKI(ISYMA)*NUMIA
                  KVMAT =  KT25  + NCKI(ISYMA)*NUMIA
                  KWMAT =  KVMAT + MXKAI*NUMIA
                  KDMAT =  KWMAT + MXKAI*NUMIA
                  KEMAT =  KDMAT + MXKAI*NUMIA
                  KENDA =  KEMAT + MXKAI*NUMIA
                  LWRKA =  LWORK - KENDA + 1
                  IF (LWRKA .LT. 0) THEN
               CALL QUIT('Insufficient memory in CC_CHOPTF1 in BATCH A')
                  ENDIF
C                     
C--------------------------------------------------- 
C                 Extract amplitudes:
C                         T1(em;k,a) = s(em,ak)
C                         T2(em;k,a) = s(ek,am)                     
C                         T3(dl;i,a) = t(dl,ai)                     
C                         T4(dl;i,a) = t(di,al)                     
C                         T5(dl;i,a) = T1(dl,i;a)
C--------------------------------------------------- 
C                     
                  CALL CCHO_FT2(T2VO,WORK(KT21),WORK(KT22),
     &                          WORK(KT23),WORK(KT24),
     &                          IA1,NUMIA,ISYMA,IOFFT)
C                     
                  CALL DCOPY(NCKI(ISYMEMK)*NUMIA,
     &                       WORK(KT21),1,WORK(KT25),1)
*     print *,'T2VO norm:',dnrm2(nt2sq(1),t2vo,1)                  
*     l = NCKI(ISYMA)*NUMIA                  
*     print *,'T21 norm:',dnrm2(l,work(kt21),1)
*     xt21total = xt21total + dnrm2(l,work(kt21),1)**2
*     print *,'T22 norm:',dnrm2(l,work(kt22),1)
*     xt22total = xt22total + dnrm2(l,work(kt22),1)**2
*     print *,'T23 norm:',dnrm2(l,work(kt23),1)
*     xt23total = xt23total + dnrm2(l,work(kt23),1)**2
*     print *,'T24 norm:',dnrm2(l,work(kt24),1)
*     xt24total = xt24total + dnrm2(l,work(kt24),1)**2
*     print *,'T25 norm:',dnrm2(l,work(kt25),1)
*     xt25total = xt25total + dnrm2(l,work(kt25),1)**2
C
                  DO ICHO = 1,NUMCHO
*     print *,'Cholesky vector No. ',icho                  
C                     
                     FCOR = ZERO
*     fcor11 = zero
*     fcor12 = zero
C                     
C----------------------------------------------------
C                    Calculate Cholesky information.
C                    d(em,#j) actual vector
C                    d(dl,#a) update vector
C                    (TODO see if CCHO_DECHO returns the occ. actual 
C                    vector and virtual update vector!!!)
C----------------------------------------------------
C                     
                     CALL CCHO_DECHOF1(FOCKD,CHOELE,NUMCHO,ICHO,
     &                               WORK(KCHOO),WORK(KCHOV),
     &                               IA1,NUMIA,ISYMA,IJ1,NUMIJ,ISYMJ)
*     print *,'d(em#j) norm: ',dnrm2(ntott1*numij,work(kchoo),1)
*     print *,'d(dl#a) norm: ',dnrm2(ntott1*numia,work(kchov),1)
C
C                        
C-------------------------------------------------------------
C                    Scale T2 amplitudes with cholesky update.
C-------------------------------------------------------------
C                        
*     print *,'T2VO norm:',dnrm2(nt2sq(1),t2vo,1)                  
*     l = NCKI(ISYMA)*NUMIA                  
*     print *,'T21 norm:',dnrm2(l,work(kt21),1)                
*     print *,'T22 norm:',dnrm2(l,work(kt22),1)                
*     print *,'T23 norm:',dnrm2(l,work(kt23),1)                
*     print *,'T24 norm:',dnrm2(l,work(kt24),1)                
*     print *,'T25 norm:',dnrm2(l,work(kt25),1)                
                     CALL CCHO_SCFT2V(WORK(KT21),WORK(KT23),
     &                            WORK(KT24),WORK(KCHOV),
     &                            NUMIA,ISYMA,NTOTT1,IOFFT,IOFFD)
*     l = NCKI(ISYMA)*NUMIA                  
*     print *,'T21 norm after scaling:',dnrm2(l,work(kt21),1) 
*     print *,'T22 norm              :',dnrm2(l,work(kt22),1)
*     print *,'T23 norm after scaling:',dnrm2(l,work(kt23),1)
*     print *,'T24 norm after scaling:',dnrm2(l,work(kt24),1)
*     print *,'T25 norm              :',dnrm2(l,work(kt25),1)
C                        
C----------------------------------------------------------------
C                    Scale Integrals with cholesky actual vector.
C----------------------------------------------------------------
C                        
                     CALL CCHO_FSCL(WORK(KITINT),WORK(KJTINT),
     &                              WORK(KIINT),WORK(KJINT),
     &                              WORK(KCHOO),NUMIJ,ISYMJ,
     &                              NTOTT1,IOFFINT,IOFFD)
*     l = ncki(isymemi)*numij
*     print *,'I(em,i#j) norm              :',dnrm2(l,work(kiint),1)
*     print *,'I(em,i#j) norm after scaling:',dnrm2(l,work(kitint),1)
*     print *,'J(em,i#j) norm              :',dnrm2(l,work(kjint),1)
*     print *,'J(em,i#j) norm after scaling:',dnrm2(l,work(kjtint),1)
C
                     DO ISYMI = 1,NSYM
*        print *, 'ISYMI = ',isymi
C                    
                        IF ( NRHF(ISYMI) .EQ. 0 ) GOTO 1236
C                       
                        ISYMEM = MULD2H(ISYMI,ISYMJ)
                        ISYMDL = MULD2H(ISYMI,ISYMA)
                        ISYMKA = ISYMEM
                        ISYMIJ = ISYMEM
                        ISYMK  = MULD2H(ISYMKA,ISYMA)
C
C                        
C-----------------------------------------------------------
C                       Calculate the V and W intermediates.
C-----------------------------------------------------------
C                       
                        KOFFJT   = KJTINT + IOFFINT(ISYMI)
                        KOFFT25  = KT25   + IOFFT(ISYMK)
C
                        NEM   = NT1AM(ISYMEM)
                        NTOEM = MAX(NEM,1)
                        NKA   = NRHF(ISYMK)*NUMIA
                        NTOKA = MAX(NKA,1)
                        NIJ   = NRHF(ISYMI)*NUMIJ
C
                        CALL DGEMM('T','N',NKA,NIJ,NEM,ONE,
     &                             WORK(KOFFT25),NTOEM,WORK(KOFFJT),
     &                             NTOEM,ZERO,WORK(KVMAT),NTOKA)
C
                        NKAIJ = NKA*NIJ
                        CALL DCOPY(NKAIJ,WORK(KVMAT),1,WORK(KWMAT),1)
C
                        KOFFIT   = KITINT + IOFFINT(ISYMI)
                        KOFFT22  = KT22   + IOFFT(ISYMK)
C
                        CALL DGEMM('T','N',NKA,NIJ,NEM,XMONE,
     &                             WORK(KOFFT22),NTOEM,
     &                             WORK(KOFFIT), NTOEM,
     &                             XMONE,WORK(KVMAT),NTOKA)
*     l = nrhf(isymk)*numia*nrhf(isymi)*numij                        
*     print *,'V norm:',dnrm2(l,work(kvmat),1)
*     xvtotal = xvtotal + dnrm2(l,work(kvmat),1)**2
C
                        CALL DGEMM('T','N',NKA,NIJ,NEM,XMONE,
     &                             WORK(KOFFT25),NTOEM,
     &                             WORK(KOFFIT), NTOEM,
     &                             TWO,WORK(KWMAT),NTOKA)
*     l = nrhf(isymk)*numia*nrhf(isymi)*numij                        
*     print *,'W norm:',dnrm2(l,work(kwmat),1)
*     xwtotal = xwtotal + dnrm2(l,work(kwmat),1)**2
C                        
C------------------------------------------------------------
C                        Calculate the D and E intermediates.
C------------------------------------------------------------
C                        
                         KOFFJ   = KJINT  + IOFFINT(ISYMK)
                         KOFFT21 = KT21   + IOFFT(ISYMI)
C
                         NDL   = NT1AM(ISYMDL)
                         NTODL = MAX(NDL,1)
                         NIA   = NRHF(ISYMI)*NUMIA
                         NTOIA = MAX(NIA,1)
                         NKJ   = NRHF(ISYMK)*NUMIJ
C
                         CALL DGEMM('T','N',NIA,NKJ,NDL,ONE,
     &                              WORK(KOFFT21),NTODL,
     &                              WORK(KOFFJ),  NTODL,
     &                              ZERO,WORK(KDMAT),NTOIA)
C
                         KOFFT23 = KT23   + IOFFT(ISYMI)
                         KOFFI   = KIINT  + IOFFINT(ISYMK)
C
                         CALL DGEMM('T','N',NIA,NKJ,NDL,XMONE,
     &                              WORK(KOFFT23),NTODL,
     &                              WORK(KOFFI),  NTODL,
     &                              ONE,WORK(KDMAT),NTOIA)
*     l = nrhf(isymk)*numia*nrhf(isymi)*numij                        
*     print *,'D norm:',dnrm2(l,work(kdmat),1)
*     xdtotal = xdtotal + dnrm2(l,work(kdmat),1)**2
C
                         KOFFT24 = KT24   + IOFFT(ISYMI)
C
                         CALL DGEMM('T','N',NIA,NKJ,NDL,XMONE,
     &                              WORK(KOFFT24),NTODL,
     &                              WORK(KOFFI),  NTODL,
     &                              ZERO,WORK(KEMAT),NTOIA) 
*     l = nrhf(isymk)*numia*nrhf(isymi)*numij                        
*     print *,'E norm:',dnrm2(l,work(kemat),1)
*     xetotal = xetotal + dnrm2(l,work(kemat),1)**2
C                         
C------------------------------------------------------
C                        Calculate energy contribution.
C------------------------------------------------------
C                         
                         DO J = 1,NUMIJ
                            DO K = 1,NRHF(ISYMK)
                               DO A = 1,NUMIA
                                  DO I = 1,NRHF(ISYMI)
C
                                  NKAIJ = NKA*NRHF(ISYMI)*(J-1)
     &                                  + NKA*(I-1) + NRHF(ISYMK)*(A-1)
     &                                  + K
                                  NIAKJ = NIA*NRHF(ISYMK)*(J-1)
     &                                  + NIA*(K-1) + NRHF(ISYMI)*(A-1)
     &                                  + I
C
                                  KOFFV = KVMAT + NKAIJ - 1
                                  KOFFW = KWMAT + NKAIJ - 1
                                  KOFFD = KDMAT + NIAKJ - 1
                                  KOFFE = KEMAT + NIAKJ - 1
C
                                  FCOR = FCOR + WORK(KOFFV)*WORK(KOFFD)
     &                                 + WORK(KOFFW)*WORK(KOFFE)
*     fcor11 = fcor11 + WORK(KOFFV)*WORK(KOFFD)
*     print *,'fcor11 ',fcor11
*     fcor12 = fcor12 + WORK(KOFFW)*WORK(KOFFE)
*     print *,'fcor12 ',fcor12
C
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO

 1236                    CONTINUE
C                        
                     ENDDO
C                    End ISYMI loop
                     !FCOR = TWO*FCOR
                     E4F1 = E4F1 + FCOR
C 
              ENERGF(ICHO) = ENERGF(ICHO) + FCOR
C 
C Decommented by Domenico
                    TNOW = SECOND()
                    DELTAT = TNOW - TLAST
                    TLAST = TNOW
                    SCNDSF(ICHO) = SCNDSF(ICHO) 
     &                          + DELTAT
C Decommented by Domenico
C 
*     e4f11 = e4f11 + fcor11                  
*     e4f12 = e4f12 + fcor12                  
*     print *,'Energy E4F1  : ',two*e4f1                  
*     print *,'Energy E4F11 : ',two*e4f11                 
*     print *,'Energy E4F12 : ',two*e4f12                       
C
c                    IF (PRINT) THEN
c                       WRITE(LUPRI,'(15X,A,I3,A,/,15X,A)')
c    &                  'Status after Cholesky vector',ICHO,':',
c    &                  '--------------------------------'
c                       IF (ABS(FCOR) .LT. THRCHO) THEN
c                          WRITE(LUPRI,'(15X,A)') 'F1 term converged'
c                       ELSE
c                          WRITE(LUPRI,'(15X,A)') 'F1 term not converged'
c                       ENDIF
c                       TIM = SECOND() - TIMT
c                       WRITE(LUPRI,'(15X,A,F10.2,A,/)')
c    &                  'Accumulated F1-time: ',TIM,' seconds'
c                    ENDIF
C
                     IF (DABS(FCOR) .LT. THRCHO) GOTO 1
C
                  ENDDO
C                 End ICHO loop
  1               CONTINUE                 
               ENDDO
C              End IBATA loop
 1235          CONTINUE              
            ENDDO
C           End ISYMA loop
         ENDDO
C        End IBATJ loop
 1234    CONTINUE
      ENDDO
C     End ISYMJ loop
*     print *,'I total norm:',dsqrt(xitotal)
*     print *,'J total norm:',dsqrt(xjtotal)
*     print *,'T21 total norm:',dsqrt(xt21total)/dsqrt(xcount)
*     print *,'T22 total norm:',dsqrt(xt22total)/dsqrt(xcount)
*     print *,'T23 total norm:',dsqrt(xt23total)/dsqrt(xcount)
*     print *,'T24 total norm:',dsqrt(xt24total)/dsqrt(xcount)
*     print *,'T25 total norm:',dsqrt(xt25total)/dsqrt(xcount)
*     print *,'V total sq norm:',xvtotal
*     print *,'W total sq norm:',xwtotal
*     print *,'D total sq norm:',xdtotal
*     print *,'E total sq norm:',xetotal
      RETURN
      END
C     
      SUBROUTINE CCHO_RDINTF1(XOINT,XIINT,XJINT,IJ1,NUMIJ,ISYMJ,IOFFINT)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     Extract integrals : I(em;ij)=(im|ej), J(em;ij)=(em|ij)
C      ordered for a fixed ISYMJ from core XOINT(mi,j;e)=(em|ij)
C      for each ISYMI all j values in the batch
C     
C  
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XOINT(*),XIINT(*),XJINT(*),IOFFINT(*)
C     
*     nsize1 = NCKI(ISYMJ)*NUMIJ
*     do i = 1,nsize1
*        xiint(i) = 999999999.0D0
*        xjint(i) = 999999999.0D0
*     enddo
*     print *,'ISJIKA : ',isjika
*     print *,'NUMIJ ',numij,'; IJ1 ',ij1
      DO ISYME = 1,NSYM
C
         ISYMJIM = ISYME
         ISYMMI  = MULD2H(ISYME,ISYMJ)
*        print *,'ISYME   =',isyme      
*        print *,'ISYMJIM =',isymjim      
*        print *,'ISYMMI  =',isymmi    
C
         DO ISYMI = 1,NSYM
C
            ISYMM = MULD2H(ISYMMI,ISYMI)
            ISYMJI = MULD2H(ISYMJ,ISYMI)
            ISYMEM = MULD2H(ISYME,ISYMM)
*        print *,'ISYMI   =',isymi      
*        print *,'ISYMM =',isymm      
*        print *,'ISYMJI  =',isymji    
*        print *,'ISYMEM  =',isymem    
C
            DO E = 1,NVIR(ISYME)
               DO J = 1,NUMIJ
                  IJ = IJ1 + J - 1
                  DO I = 1,NRHF(ISYMI)
                     DO M = 1,NRHF(ISYMM)
C
                        KJIME = ISJIKA(ISYMJIM,ISYME)
     &                        + NMAIJK(ISYMJIM)*(E - 1)
     &                        + IMAIJK(ISYMJI,ISYMM)
     &                        + NMATIJ(ISYMJI)*(M - 1)
     &                        + IMATIJ(ISYMJ,ISYMI)
     &                        + NRHF(ISYMJ)*(I - 1) 
     &                        + IJ
                        KMIJE = ISJIKA(ISYMJIM,ISYME)
     &                        + NMAIJK(ISYMJIM)*(E - 1)
     &                        + IMAIJK(ISYMMI,ISYMJ)
     &                        + NMATIJ(ISYMMI)*(IJ - 1)
     &                        + IMATIJ(ISYMM,ISYMI)
     &                        + NRHF(ISYMM)*(I - 1)
     &                        + M
*     print *,'ISJIKA(ISYMJIM,ISYME)',ISJIKA(ISYMJIM,ISYME)
*     print *,'NMAIJK(ISYMJIM)',NMAIJK(ISYMJIM)
*     print *,'IMAIJK(ISYMMI,ISYMJ)',IMAIJK(ISYMMI,ISYMJ) 
*     print *,'NMATIJ(ISYMMI)',NMATIJ(ISYMMI) 
*     print *,'IJ - 1', IJ - 1 
*     print *,'NRHF(ISYMM)',NRHF(ISYMM)
C
                        NIJ = NRHF(ISYMI)*(J - 1) + I
                        NEM = IT1AM(ISYME,ISYMM)+NVIR(ISYME)*(M - 1) + E
                        KEMIJ = IOFFINT(ISYMI) 
     &                        + NT1AM(ISYMEM)*(NIJ - 1) + NEM 
*     print *,'IOFFINT(ISYMI)',IOFFINT(ISYMI)
*     print *,'NT1AM(ISYMEM)',NT1AM(ISYMEM)
C
                        XIINT(KEMIJ) = XOINT(KJIME)
                        XJINT(KEMIJ) = XOINT(KMIJE)
*     print *,'KEMIJ =',kemij,';KMIJE = ',kmije,'E,J,I,M',e,j,i,m
C
*     nsize1 = NCKI(ISYMJ)*NUMIJ
*     nsize2 = NTRAOC(1) 
*     if (kemij .gt. nsize1 .or. kjime .gt. nsize2 .or. 
*    &    kmije .gt. nsize2 ) print *, 'OUT OF BOUNDARIES!!!'
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
C           
         ENDDO
      ENDDO
C     
*     print *,'XIINT ------------'
*     write (*,*) (xiint(i),i=1,nsize1)
*     print *,'XJINT ------------'
*     write (*,*) (xjint(i),i=1,nsize1)
*     print *,'XOINT ------------'
*     write (*,*) (xoint(i),i=1,nsize2)
      RETURN
      END
C
      SUBROUTINE CCHO_FT2(T2VO,T21,T22,T23,T24,IA1,NUMIA,ISYMA,IOFFT)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C     
C     Extract amplitudes:
C                         T1(em;ka) = s(em,ak) = 2T(ea,mk) - T(ea,km)
C                         T2(em;ka) = s(ek,am) = 2T(ea,km) - T(ea,mk)
C                         T3(em;ka) = t(em,ak) = T(ea,mk)
C                         T4(em;ka) = t(ek,am) = T(ea,km)
C     
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO = 0.0D0, XMONE = -1.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION T2VO(*),T21(*),T22(*),T23(*),T24(*),IOFFT(*)
C     
      DO ISYMK = 1,NSYM
C
         ISYMAK = MULD2H(ISYMK,ISYMA)
         ISYMEM = ISYMAK  
C
         DO ISYMM = 1,NSYM
C
            ISYMMK = MULD2H(ISYMM,ISYMK)
            ISYMEA = ISYMMK
            ISYME  = MULD2H(ISYMEA,ISYMA)
C
            DO K = 1,NRHF(ISYMK)
               DO M = 1,NRHF(ISYMM)
                  DO A = 1,NUMIA
                     IA = IA1 + A - 1
                     DO E = 1,NVIR(ISYME)
C
                        NMK = IMATIJ(ISYMM,ISYMK) 
     &                      + NRHF(ISYMM)*(K - 1) + M
                        NKM = IMATIJ(ISYMK,ISYMM)
     &                      + NRHF(ISYMK)*(M - 1) + K
                        NEA = IMATAB(ISYME,ISYMA)
     &                      + NVIR(ISYME)*(IA - 1) + E
                        NEAMK = IT2VO(ISYMEA,ISYMMK) 
     &                        + NMATAB(ISYMEA)*(NMK - 1)
     &                        + NEA
                        NEAKM = IT2VO(ISYMEA,ISYMMK) 
     &                        + NMATAB(ISYMEA)*(NKM - 1)
     &                        + NEA
                        NKA = NRHF(ISYMK)*(A - 1) + K
                        NEM = IT1AM(ISYME,ISYMM)
     &                      + NVIR(ISYME)*(M - 1) + E
                        NEMKA = IOFFT(ISYMK) + NT1AM(ISYMEM)*(NKA - 1)
     &                        + NEM
C
                        T23(NEMKA) = T2VO(NEAMK)
                        T24(NEMKA) = T2VO(NEAKM)
C
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
C           
         ENDDO
      ENDDO
C
      NTOT = NCKI(ISYMA)*NUMIA
C
      CALL DCOPY(NTOT,T24,1,T21,1)
      CALL DSCAL(NTOT,XMONE,T21,1)
      CALL DAXPY(NTOT,TWO,T23,1,T21,1)
C
      CALL DCOPY(NTOT,T23,1,T22,1)
      CALL DSCAL(NTOT,XMONE,T22,1)
      CALL DAXPY(NTOT,TWO,T24,1,T22,1)
C
      RETURN
      END
C
C  /* Deck ccho_scft2v */    
      SUBROUTINE CCHO_SCFT2V(T21,T23,T24,VCCHO,NUMIA,ISYMA,NTOTT1,
     &                       IOFFT,IOFFD)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     Scale T21 T23 T24 amplitudes with cholesky update
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T21(*),T23(*),T24(*),VCCHO(*),IOFFT(8),IOFFD(8)
#include <ccorb.h>
#include <ccsdsym.h>
      DO ISYMI = 1,NSYM
C
         ISYMDL = MULD2H(ISYMI,ISYMA)
C
         DO A = 1,NUMIA
            DO I = 1,NRHF(ISYMI)
               DO IDL = 1,NT1AM(ISYMDL)
C
                  NDLA  = NTOTT1*(A -  1) + IOFFD(ISYMDL) + IDL
                  NIA = NRHF(ISYMI)*(A - 1) + I
                  NDLIA = IOFFT(ISYMI) + NT1AM(ISYMDL)*(NIA - 1) + IDL 
C
                  T21(NDLIA) = VCCHO(NDLA)*T21(NDLIA)
                  T23(NDLIA) = VCCHO(NDLA)*T23(NDLIA)
                  T24(NDLIA) = VCCHO(NDLA)*T24(NDLIA)
C
               ENDDO
            ENDDO
         ENDDO
C        
      ENDDO
C     
      RETURN
      END 
C
C  /* Deck ccho_fscl */    
      SUBROUTINE CCHO_FSCL(XITINT,XJTINT,XIINT,XJINT,CCHOO,NUMIJ,ISYMJ,
     &                     NTOTT1,IOFFINT,IOFFD)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     Scale I and J integrals with cholesky actual vector
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION XITINT(*),XJTINT(*),XIINT(*),XJINT(*),
     &          CCHOO(*),IOFFINT(8),IOFFD(8)
#include <ccorb.h>
#include <ccsdsym.h>
      DO ISYMI = 1,NSYM
C
         ISYMEM = MULD2H(ISYMI,ISYMJ)
C
         DO J = 1,NUMIJ
            DO I = 1,NRHF(ISYMI)
               DO IEM = 1,NT1AM(ISYMEM)
C
                  NEMJ  = NTOTT1*(J - 1) + IOFFD(ISYMEM) + IEM
                  NIJ   = NRHF(ISYMI)*(J - 1) + I
                  NEMIJ = IOFFINT(ISYMI) + NT1AM(ISYMEM)*(NIJ - 1) + IEM
C
                  XITINT(NEMIJ) = CCHOO(NEMJ)*XIINT(NEMIJ)
                  XJTINT(NEMIJ) = CCHOO(NEMJ)*XJINT(NEMIJ)
C
               ENDDO
            ENDDO
         ENDDO
C
      ENDDO
C
      RETURN
      END
C  /* Deck ccho_dechof1 */
      SUBROUTINE CCHO_DECHOF1(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                      IA1,NUMIA,ISYMA,IJ1,NUMIJ,ISYMJ)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez Oct 2002.
C
C     Construct vector for Cholesky decomposition:
C
C     For the occupied part, OCCHO(em#j), the actual vector
C     For the virtual  part, VICHO(dl#a), the update vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C-------------------
C     Occupied part.
C-------------------
C
C     Get the omega's.
C     ----------------
C
      IND = 0
      DO J = 1,NUMIJ
         KOFFJ = IRHF(ISYMJ) + IJ1 + J - 1
         DO ISYMEM = 1,NSYM
            DO ISYMM = 1,NSYM
               ISYME = MULD2H(ISYMEM,ISYMM)
               DO M = 1,NRHF(ISYMM)
                  KOFFM = IRHF(ISYMM) + M
                  DO E = 1,NVIR(ISYME)
                     KOFFE = IVIR(ISYME) + E
C
                     IND = IND + 1
                     OCCHO(IND) = FOCKD(KOFFE) - FOCKD(KOFFM) 
     &                          - FOCKD(KOFFJ)
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      NDIMVI = IND
C
C
C     Construct the vector.
C     ---------------------
C
      DO P = 1,NDIMVI
C
         OMEGA = OCCHO(P)
C
         OCCHO(P) = SQRT(TWO*CHOELE(ICHO))/(OMEGA+CHOELE(ICHO))
C
         DO JCHO = 1,ICHO-1
            OCCHO(P) = OCCHO(P) 
     &               *(OMEGA-CHOELE(JCHO))/(OMEGA+CHOELE(JCHO))
         END DO
C
      END DO
C
C
C------------------
C     Virtual part.
C------------------
C
      IF (ICHO .EQ. 1) THEN
C
C        First vector.
C        --------------
C
         IND = 0
         DO A = 1,NUMIA
            KOFFA = IVIR(ISYMA) + IA1 + A - 1
            DO ISYMDL = 1,NSYM
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  DO L = 1,NRHF(ISYML)
                     KOFFL = IRHF(ISYML) + L
                     DO D = 1,NVIR(ISYMD)
                        KOFFD = IVIR(ISYMD) + D
                        OME = FOCKD(KOFFD)-FOCKD(KOFFL)+FOCKD(KOFFA)
C
                        IND = IND + 1
                        VICHO(IND) = 
     &                       SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     ENDDO
                  END DO
               END DO
            END DO
         END DO
C
      ELSE
C
C        Updating vector.
C        ----------------
C
         IND = 0
         DO A = 1,NUMIA
            KOFFA = IVIR(ISYMA) + IA1 + A - 1
            DO ISYMDL = 1,NSYM
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  DO L = 1,NRHF(ISYML)
                     KOFFL = IRHF(ISYML) + L
                     DO D = 1,NVIR(ISYMD)
                        KOFFD = IVIR(ISYMD) + D
                        OME = FOCKD(KOFFD)-FOCKD(KOFFL)+FOCKD(KOFFA)
C
                        IND = IND + 1
                        VICHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     ENDDO
                  END DO
               END DO
            END DO
         END DO
C
         NDIMOC = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMOC,FACTOR,VICHO,1)
C
      END IF
C
      RETURN
      END      
C  /* Deck ccho_fterm2_old */
      SUBROUTINE CCHO_FTERM2_OLD(XOINT,T2VO,FOCKD,NUMCHO,CHOELE,
     &                           WORK,LWORK,E4F2)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     Calculate (part of F-terms):
C
C        E4F2 = 2*2 Sum(a) Sum(ijk) [ P(ja,ki) * Q(jk,ia)
C                               - R(ja,ki) * S(jk,ia) ]
C
C     where
C        
C        P(ja,ki) =     Sum(bc) d(bjc) * M(bc,j;a) * t(bc,ki)
C        R(ja,ki) =     Sum(bc) d(bjc) * N(bc,j;a) * t(bc,ki) 
C        Q(jk,ia) =   - Sum(dl) d(dli) * T1(dl,i;a) * I(dl,kj)
C                 + 2 * Sum(dl) d(dli) * T2(dl,i;a) * I(dl,kj)
C                 -     Sum(dl) d(dli) * T2(dl,i;a) * J(dl,kj)
C        S(jk,ia) =     Sum(dl) d(dli) * T1(dl,i;a) * J(dl,kj)
C
C     and M(bc,j;a)  = 2(bj|ca) - (cj|ba),
C         N(bc,j;a)  = 2(cj|ba) - (bj|ca),
C         T1(dl,i;a) = t(ad,li) 
C         T2(dl,i;a) = t(ad,il) 
C         I(dl,kj)   = (jk|ld) 
C         J(dl,kj)   = (lk|jd)
C     and d(dli) and d(bjc) denotes the occupied and virtual parts of
C     the Cholesky decomposition of the orbital energy denominator.
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XOINT(*),T2VO(*),FOCKD(*),CHOELE(*)
      DIMENSION WORK(LWORK)
      DIMENSION IOFF2(8)
      CHARACTER*7 FILNAM
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>      
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR  =  4.0D0)
C 
      TLAST = SECOND()
C 
*     dimension xchoototal(50),xchovtotal(50)
*     dimension xmttotal(50),xnttotal(50) 
*     dimension xt21ttotal(50),xt22ttotal(50)
*     call dzero(50,xchoototal)
*     call dzero(50,xchovtotal)
*     call dzero(50,xmttotal)
*     call dzero(50,xnttotal)
*     call dzero(50,xt21ttotal)
*     call dzero(50,xt22ttotal)
C
      ISYMT2 = 1
C     
      MXJKI = -1
      DO ISYMI = 1,NSYM
         DO ISYMK = 1,NSYM
            DO ISYMJ = 1,NSYM
               NJKI   = NRHF(ISYMJ)*NRHF(ISYMK)*NRHF(ISYMI)
               IF (NJKI .GT. MXJKI) THEN
                  MXJKI = NJKI
               ENDIF
            ENDDO
         ENDDO
      ENDDO
*     print *,'MXJKI = ',mxjki      
C      
C--------------------------------------------------------
C     IOFF2(ISYMDL) --> offset in OP(dl,kj) (ISYMOP = 1).
C--------------------------------------------------------
C      
      ICOUNT = 0
      DO ISYM = 1,NSYM
         IOFF2(ISYM) = ICOUNT
         ICOUNT = ICOUNT + NT1AM(ISYM)*NMATIJ(ISYM)
      ENDDO
*     print *,IOFF2      
C         
C---------------------------
C        Dynamic allocation.
C---------------------------
C         
         KIINT = 1
         KJINT = KIINT + NCKIJ(1)
         KPMAT = KJINT + NCKIJ(1)
         KQMAT = KPMAT + MXJKI
         KRMAT = KQMAT + MXJKI
         KSMAT = KRMAT + MXJKI
         KEND1 = KSMAT + MXJKI
         LWRK1 = LWORK - KEND1 +1
*     print *,'NCKIJ(1)',NCKIJ(1)        
*     print *,'KIINT = ',kiint         
*     print *,'KJINT = ',kjint         
*     print *,'KPMAT = ',kpmat
*     print *,'KQMAT = ',kqmat         
*     print *,'KRMAT = ',krmat         
*     print *,'KSMAT = ',ksmat         
*     print *,'KEND1 = ',kend1         
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient memory in CC_CHOPTF2')
         ENDIF
C         
C----------------------------------------
C        Extract I(dl,kj) and J(dl,kj).
C----------------------------------------
C
*     print *,'size xoint = ',NTRAOC(1)         
*     print *,'xoint norm:',dnrm2(NTRAOC(1),xoint,1)          
         CALL CCHO_RDINIJF2(XOINT,WORK(KIINT),WORK(KJINT),IOFF2)
*     print *,'xoint norm after rd:',dnrm2(NTRAOC(1),xoint,1)
*     l = nckij(1)
*     print *,'size I and J integrals = ',l
*     print *,'I norm:',dnrm2(l,work(kiint),1)
*     print *,'J norm:',dnrm2(l,work(kjint),1)
C
         DO ISYMA = 1,NSYM
*     print *,'ISYMA = ',isyma         
*     print *,'NBCJ(',isyma,') = ',nckatr(isyma)         
*     print *,'NVIR(',ISYMA,') = ',NVIR(ISYMA)
C
*     if ( nvir(isyma) .eq. 0 ) then
*        do icho = 1,numcho
*        ISYMBCJ = ISYMA
*        ISYMDLI = ISYMA
*        CALL CCHO_DECHOF2(FOCKD,CHOELE,NUMCHO,ICHO,
*    &                     WORK(KCHOO),WORK(KCHOV),
*    &                     ISYMDLI,ISYMBCJ)
*        print *,'d(dl,i) norm: ',dnrm2(ncki(isymdli),work(kchoo),1)
*        print '(4F20.15)',(work(ii),ii=kchoo,kchoo+ncki(isymdli)-1) 
*        print *,'d(bc,j) norm: ',dnrm2(nckatr(isymbcj),work(kchov),1)
*        print '(4F20.15)',(work(ii),ii=kchov,kchov+nckatr(isymbcj)-1) 
*        xchoototal(icho) = xchoototal(icho) 
*    &           + dnrm2(ncki(isymdli),work(kchoo),1)**2
*        xchovtotal(icho) = xchovtotal(icho)
*    &           + dnrm2(nckatr(isymbcj),work(kchov),1)**2
*        enddo
*     endif
            IF ( NVIR(ISYMA) .EQ. 0 ) GOTO 1234
C
            ISYMBCJ = ISYMA
            ISYMDLI = ISYMA
C            
C-----------------------------------------------
C           Allocate space for Cholesky vectors.
C-----------------------------------------------
C            
            LENCHOV = NCKATR(ISYMBCJ)
            LENCHOO = NCKI(ISYMDLI)
C
            KCHOV = KEND1
            KCHOO = KCHOV + LENCHOV
            KEND2 = KCHOO + LENCHOO
            LWRK2 = LWORK - KEND2 + 1
*     print *,'KCHOV = ',kchov,'  KCHOO = ',kchoo,' KEND2 = ',kend2
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Not enough space in CC_CHOPTF2')
            ENDIF
C            
C------------------------
C           Batch over A.
C------------------------
C            
            LENA = 2*NCKATR(ISYMBCJ) + 2*NCKI(ISYMDLI)
            NUMA = MIN(NVIR(ISYMA),LWRK2/LENA)
            IF (NUMA .EQ. 0) THEN
               WRITE(LUPRI,*) 'NUMA .EQ. 0 in CC_CHOPTF2'
               CALL QUIT('Not enough space in CC_CHOPTF2')
            ENDIF
C           
            NBATA = (NVIR(ISYMA)-1)/NUMA + 1
C                 
            IA2 = 0
            DO IBATA = 1,NBATA
C                 
               IA1 = IA2 + 1
               IA2 = IA2 + NUMA
               IF (IA2 .GT. NVIR(ISYMA)) IA2 = NVIR(ISYMA)
               NUMIA = IA2 - IA1 + 1
*     print *,'NUMIA =',numia
*     print *,'NBATA =',nbata
C                    
C---------------------------------
C              Dynamic allocation.
C---------------------------------
C              
               KMINT = KEND2
               KNINT = KMINT  + NCKATR(ISYMBCJ)*NUMIA
               KT21  = KNINT  + NCKATR(ISYMBCJ)*NUMIA
               KT22  = KT21   + NCKI(ISYMDLI)*NUMIA
               KENDA = KT22   + NCKI(ISYMDLI)*NUMIA
               LWRKA = LWORK  - KENDA + 1
*     print *,'KMINT = ',kmint               
*     print *,'KNINT = ',knint               
*     print *,'KT21  = ',kt21
*     print *,'KT22  = ',kt22
*     print *,'KENDA = ',kenda
               IF (LWRKA .LT. 0) THEN
               CALL QUIT('Insufficient memory in CC_CHOPTF2 in BATCH A')
               ENDIF
C               
C---------------------------------------------------------
C              Read the integrals M(bc,j;a) and N(bc,j;a).
C---------------------------------------------------------
C              
               FILNAM = 'CHO_VI1'
               LUNIT  = 61 
               CALL CCHO_RDINMNF2(LUNIT,FILNAM,WORK(KMINT),WORK(KNINT),
     &                        IA1,NUMIA,ISYMA)
*     l = NCKASR(ISYMBCJ)*NUMIA               
*     print *,'M norm: ',dnrm2(l,work(kmint),1)
*     print *,'N norm: ',dnrm2(l,work(knint),1)
*     xmtotal = xmtotal + dnrm2(l,work(kmint),1)**2
*     xntotal = xntotal + dnrm2(l,work(knint),1)**2
C               
C------------------------------------------------
C              Extract T1(dl,i;a) and T2(dl,i;a).
C------------------------------------------------
C               
               CALL CCHO_T2F2(T2VO,WORK(KT21),WORK(KT22),
     &                        IA1,NUMIA,ISYMA)    
*     print *,'T2VO norm:',dnrm2(nt2sq(1),t2vo,1)                  
*     l = NCKI(ISYMA)*NUMIA                  
*     print *,'T21 norm:',dnrm2(l,work(kt21),1)
*     print *,'T22 norm:',dnrm2(l,work(kt22),1)
*     print *,'T21: '
*     print *,(work(i),i=kt21,kt21+l)
*     print *,'T22: '
*     print *,(work(i),i=kt22,kt22+l)
*     xt21total = xt21total + dnrm2(l,work(kt21),1)**2
*     xt22total = xt22total + dnrm2(l,work(kt22),1)**2
C              
               DO ICHO = 1,NUMCHO
*     print *,'Cholesky vector No. ',icho               
C
                  FCOR = ZERO
*     fcor21 = zero                  
*     fcor22 = zero                  
C                  
C------------------------------------------------
C                 Calculate Cholesky information.
C                 d(dl,i) update vector
C                 d(bc,j) update vector
C                 (TODO  check if it is right)
C------------------------------------------------
C 
*     print *,'Cholesky dimensions:',ncki(isymdli),nckatr(isymbcj)
*     print *,'ISYMA = ',isyma,'; ISYMDLI = ',isymdli,
*    &        '; ISYMBCJ = ',isymbcj
*     call ccho_decho3(fockd,choele,numcho,icho,work(kchoo),work(kchov),
*    &                 isymdli,isymbcj,1,1)
*     print *,'d(dl,i) dch3 norm: ',dnrm2(ncki(isymdli),work(kchoo),1)
*     print *,(work(ii),ii=kchoo,kchoo+ncki(isymdli)-1)
*     print *,'d(bc,j) dch3 norm: ',dnrm2(nckatr(isymbcj),work(kchov),1)
*     call ccho_dechof1(fockd,choele,numcho,icho,work(kchoo),
*    &                  work(kchov),1,nvirt,1,1,nrhft,1)
*     print *,'d(dl,i) dchf1 norm: ',dnrm2(ncki(1),work(kchoo),1)
*     print *,'d(bc,j) dchf1 norm: ',dnrm2(nckatr(1),work(kchov),1)
                  CALL CCHO_DECHOF2(FOCKD,CHOELE,NUMCHO,ICHO,
     &                              WORK(KCHOO),WORK(KCHOV),
     &                              ISYMDLI,ISYMBCJ)
*     print *,'d(dl,i) norm: ',dnrm2(ncki(isymdli),work(kchoo),1)
*     print '(4F20.15)',(work(ii),ii=kchoo,kchoo+ncki(isymdli)-1) 
*     print *,'d(bc,j) norm: ',dnrm2(nckatr(isymbcj),work(kchov),1)
*     print '(4F20.15)',(work(ii),ii=kchov,kchov+nckatr(isymbcj)-1) 
*     xchoototal(icho) = xchoototal(icho) 
*    &           + dnrm2(ncki(isymdli),work(kchoo),1)**2
*     xchovtotal(icho) = xchovtotal(icho)
*    &           + dnrm2(nckatr(isymbcj),work(kchov),1)**2
*     print *,'SETTTING CHOLESKY VECTORS TO 1 !!!!'
*     do ii = 1,ncki(isymdli) 
*        kof = kchoo + ii - 1
*        work(kof) = one
*     enddo
*     do ii = 1,nckatr(isymbcj) 
*        kof = kchov + ii - 1
*        work(kof) = one
*     enddo
C                  
C--------------------------------------------------------------
C                 Scale M and N integrals with cholesky update.
C--------------------------------------------------------------
C                  
                  CALL CCHO_FSCMN(WORK(KMINT),WORK(KNINT),WORK(KCHOV),
     &                            NUMIA,ISYMA)
*     l = NCKASR(ISYMBCJ)*NUMIA               
*     print *,'M norm after scaling: ',dnrm2(l,work(kmint),1)
*     print *,'N norm after scaling: ',dnrm2(l,work(knint),1)
*     xmttotal(icho) = xmttotal(icho) 
*    &           + dnrm2(l,work(kmint),1)**2
*     xnttotal(icho) = xnttotal(icho)
*    &           + dnrm2(l,work(knint),1)**2
C                  
C----------------------------------------------------------
C                 Scale T2 amplitudes with cholesky update.
C----------------------------------------------------------
C                  
                  CALL CCHO_F2SCT2O(WORK(KT21),WORK(KT22),WORK(KCHOO),
     &                              NUMIA,ISYMA) 
*     print *,'T2VO norm after f2sct2o:',dnrm2(nt2sq(1),t2vo,1)
*     l = NCKI(ISYMA)*NUMIA                  
*     print *,'T21 norm after scaling:',dnrm2(l,work(kt21),1)
*     print *,'T22 norm after scaling:',dnrm2(l,work(kt22),1)
*     xt21ttotal(icho) = xt21ttotal(icho) 
*    &           + dnrm2(l,work(kt21),1)**2
*     xt22ttotal(icho) = xt22ttotal(icho)
*    &           + dnrm2(l,work(kt22),1)**2
C                 
                  DO ISYMJ = 1,NSYM
*     print *,'ISYMJ = ',isymj         
C                 
                     IF ( NRHF(ISYMJ) .EQ. 0 ) GOTO 1235
C
                     ISYMBC  = MULD2H(ISYMJ,ISYMA)
                     ISYMKI  = ISYMBC
                     ISYMDLK = ISYMJ
C
                     DO ISYMK = 1,NSYM
*     print *,'ISYMK = ',isymk         
C                 
                        IF ( NRHF(ISYMK) .EQ. 0 ) GOTO 1236
C
                        ISYMI   = MULD2H(ISYMKI,ISYMK)
                        ISYMDL  = MULD2H(ISYMI,ISYMA)
                        ISYMKJ  = ISYMDL
C
                        DO A = 1,NUMIA
*     print *,'A = ',a
C                        
C----------------------------------------------------------
C                          Calculate P and R intermediates.
C----------------------------------------------------------
C                        
                           KOFFT = IT2VO(ISYMBC,ISYMKI) 
     &                           + NMATAB(ISYMBC)*IMATIJ(ISYMK,ISYMI)
     &                           + 1                           
                           KOFFM = KMINT + NCKASR(ISYMBCJ)*(A - 1)
     &                           + ICKASR(ISYMBC,ISYMJ)
                           KOFFN = KNINT + NCKASR(ISYMBCJ)*(A - 1)
     &                           + ICKASR(ISYMBC,ISYMJ)
*     print *,'KOFFT = ',kofft,';KOFFM = ',koffm,';KOFFN = ',koffn
C
                           NBC   = NMATAB(ISYMBC)
                           NTOBC = MAX(1,NBC)
                           NJ    = NRHF(ISYMJ)
                           NTOJ  = MAX(1,NJ)
                           NKI   = NRHF(ISYMK)*NRHF(ISYMI)
*     print *,'NBC = ',nbc,';NTOBC = ',ntobc,';NJ = ',nj
*     print *,'NTOJ = ',ntoj,';NKI = ',nki
C
*     write (LUPRI,*) 'Calculating P'                           
                           CALL DGEMM('T','N',NJ,NKI,NBC,ONE,
     &                                WORK(KOFFM),NTOBC,
     &                                T2VO(KOFFT),NTOBC,
     &                                ZERO,WORK(KPMAT),NTOJ)
*     l = nrhf(isymj)*nrhf(isymk)*nrhf(isymi)
*     print *,'P dimension ',l
*     print *,'P norm: ',dnrm2(l,work(kpmat),1)             
C
*     write (LUPRI,*) 'Calculating R'                           
                           CALL DGEMM('T','N',NJ,NKI,NBC,ONE,
     &                                WORK(KOFFN),NTOBC,
     &                                T2VO(KOFFT),NTOBC,
     &                                ZERO,WORK(KRMAT),NTOJ)
*     l = nrhf(isymj)*nrhf(isymk)*nrhf(isymi)                           
*     print *,'R norm: ',dnrm2(l,work(krmat),1)             
*     l = NCKASR(ISYMBCJ)*NUMIA               
*     print *,'M norm after R: ',dnrm2(l,work(kmint),1)
*     print *,'N norm after R: ',dnrm2(l,work(knint),1)
C                           
C----------------------------------------------------------
C                          Calculate Q and S intermediates.
C----------------------------------------------------------
C                           
                           KOFFI   = KIINT + IOFF2(ISYMDL)
     &                             + NT1AM(ISYMDL)*IMATIJ(ISYMK,ISYMJ)
                           KOFFJ   = KJINT + IOFF2(ISYMDL)
     &                             + NT1AM(ISYMDL)*IMATIJ(ISYMK,ISYMJ)
                           KOFFT21 = KT21 + NCKI(ISYMDLI)*(A - 1) 
     &                             + ICKI(ISYMDL,ISYMI)
                           KOFFT22 = KT22 + NCKI(ISYMDLI)*(A - 1)
     &                             + ICKI(ISYMDL,ISYMI)
C
                           NDL   = NT1AM(ISYMDL)
                           NTODL = MAX(1,NDL)
                           NKJ   = NRHF(ISYMK)*NRHF(ISYMJ)
                           NI    = NRHF(ISYMI)
                           NTOI  = MAX(1,NI)
*     print *,'NDL = ',ndl,';NTODL = ',ntodl,';NI = ',nj
*     print *,'NTOI = ',ntoj,';NKJ = ',nkj
C
*     print *,'Calculating Q'                           
                           CALL DGEMM('T','N',NI,NKJ,NDL,XMONE,
     &                                WORK(KOFFT21),NTODL,
     &                                WORK(KOFFI),NTODL,
     &                                ZERO,WORK(KQMAT),NTOI)
C
*     print *,'Calculating Q'                           
                           CALL DGEMM('T','N',NI,NKJ,NDL,TWO,
     &                                WORK(KOFFT22),NTODL,
     &                                WORK(KOFFI),NTODL,
     &                                ONE,WORK(KQMAT),NTOI)
C
*     print *,'Calculating Q'                           
                           CALL DGEMM('T','N',NI,NKJ,NDL,XMONE,
     &                                WORK(KOFFT22),NTODL,
     &                                WORK(KOFFJ),NTODL,
     &                                ONE,WORK(KQMAT),NTOI)
*     l = nrhf(isymj)*nrhf(isymk)*nrhf(isymi)                           
*     print *,'Q norm: ',dnrm2(l,work(kqmat),1)             
C
*     print *,'Calculating S'                           
                           CALL DGEMM('T','N',NI,NKJ,NDL,ONE,
     &                                WORK(KOFFT21),NTODL,
     &                                WORK(KOFFJ),NTODL,
     &                                ZERO,WORK(KSMAT),NTOI)
*     l = nrhf(isymj)*nrhf(isymk)*nrhf(isymi)                           
*     print *,'S norm: ',dnrm2(l,work(ksmat),1)             
C                           
C--------------------------------------------------------
C                          Calculate energy contribution.
C--------------------------------------------------------
C                           
                           DO I = 1,NRHF(ISYMI)
                              DO K = 1,NRHF(ISYMK)
                                 DO J = 1,NRHF(ISYMJ)
C
                                    NJKI = NRHF(ISYMJ)*NRHF(ISYMK)*(I-1)
     &                                   + NRHF(ISYMJ)*(K - 1) + J
                                    NIKJ = NRHF(ISYMI)*NRHF(ISYMK)*(J-1)
     &                                   + NRHF(ISYMI)*(K - 1) + I
C
                                    KOFFP = KPMAT + NJKI - 1
                                    KOFFR = KRMAT + NJKI - 1
                                    KOFFQ = KQMAT + NIKJ - 1
                                    KOFFS = KSMAT + NIKJ - 1
C
                                    FCOR = FCOR
     &                                   + WORK(KOFFP)*WORK(KOFFQ)
     &                                   - WORK(KOFFR)*WORK(KOFFS)
*     fcor21 = fcor21 + WORK(KOFFP)*WORK(KOFFQ)
*     fcor22 = fcor22 - WORK(KOFFR)*WORK(KOFFS)
C
                                 ENDDO
                              ENDDO
                           ENDDO
C
                        ENDDO
C                       End NUMIA loop
 1236                   CONTINUE
                     ENDDO
C                    End ISYMK loop
 1235                CONTINUE
                  ENDDO
C                 End ISYMJ loop
C
                  FCOR = TWO * FCOR
*     fcor21 = two*fcor21                  
*     fcor22 = two*fcor22                  
*     print *,'Energy correction fcor  : ',fcor                  
*     print *,'Energy correction fcor21: ',fcor21                  
*     print *,'Energy correction fcor22: ',fcor22                 
C
                  E4F2 = E4F2 + FCOR
C 
              ENERGF(ICHO) = ENERGF(ICHO) + FCOR
C 
C Decommented by Domenico
                 TNOW = SECOND()
                 DELTAT = TNOW - TLAST
                 TLAST = TNOW
                 SCNDSF(ICHO) = SCNDSF(ICHO) 
     &                         + DELTAT
C Decommented by Domenico
C 
*     e4f21 = e4f21 + fcor21                  
*     e4f22 = e4f22 + fcor22                  
*     print *,'Energy E4F2  : ',e4f2                  
*     print *,'Energy E4F21 : ',e4f21                 
*     print *,'Energy E4F22 : ',e4f22                 
C
                  IF (DABS(FCOR) .LT. THRCHO) GOTO 2
C
               ENDDO
C              End ICHO loop
    2          CONTINUE               
            ENDDO
C           End IBATA loop
 1234       CONTINUE
         ENDDO
C        End ISYMA loop
*     print *,'d(dli) sq norms:',(xchoototal(i),i=1,numcho)
*     print *,'d(bcj) sq norms:',(xchovtotal(i),i=1,numcho)
*     print *,'M sq norm:',xmtotal
*     print *,'N sq norm:',xntotal
*     print *,'T21 sq norm:',xt21total
*     print *,'T22 sq norm:',xt22total
*     print *,'M scal sq norm:'
*     print 1,(xmttotal(i),i=1,numcho)
*     print *,'N scal sq norm:'
*     print 1,(xnttotal(i),i=1,numcho)
*     print *,'T21 scal sq norm:'
*     print 1,(xt21ttotal(i),i=1,numcho)
*     print *,'T22 scal sq norm:'
*     print 1,(xt22ttotal(i),i=1,numcho)
* 1   format(4D20.13)
C
*     print *,'FINAL Energy E4F2  : ',two*e4f2                  
*     print *,'FINAL Energy E4F21 : ',two*e4f21                 
*     print *,'FINAL Energy E4F22 : ',two*e4f22                 
      RETURN
      END
C
C  /* Deck ccho_rdinijf2_old */
      SUBROUTINE CCHO_RDINIJF2_OLD(XOINT,XIINT,XJINT,IOFF2) 
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     
C     Extract occupied integrals: 
C                    I(dl,kj) = (jk|ld) 
C                    J(dl,kj) = (lk|jd) 
C     from XOINT(lj,k;d) = (jk|ld) 
C     
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DIMENSION XOINT(*),XIINT(*),XJINT(*),IOFF2(*) 
C     
      DO ISYMD = 1,NSYM
         ISYMLJK = ISYMD
         DO ISYMK = 1,NSYM
            ISYMLJ = MULD2H(ISYMLJK,ISYMK)
            DO ISYMJ = 1,NSYM
               ISYML = MULD2H(ISYMLJ,ISYMJ)
               ISYMDLK = ISYMJ
               ISYMDL = MULD2H(ISYMDLK,ISYMK)
C
               DO D = 1,NVIR(ISYMD)
                  DO K = 1,NRHF(ISYMK)
                     DO J = 1,NRHF(ISYMJ)
                        DO L = 1,NRHF(ISYML)
C
                           KOFFLJKD = ISJIKA(ISYMLJK,ISYMD)
     &                           + NMAIJK(ISYMLJK)*(D - 1)
     &                           + IMAIJK(ISYMLJ,ISYMK)
     &                           + NMATIJ(ISYMLJ)*(K - 1)
     &                           + IMATIJ(ISYML,ISYMJ)
     &                           + NRHF(ISYML)*(J - 1)
     &                           + L
C
                           KOFFJLKD = ISJIKA(ISYMLJK,ISYMD)
     &                           + NMAIJK(ISYMLJK)*(D - 1)
     &                           + IMAIJK(ISYMLJ,ISYMK)
     &                           + NMATIJ(ISYMLJ)*(K - 1)
     &                           + IMATIJ(ISYMJ,ISYML)
     &                           + NRHF(ISYMJ)*(L - 1)
     &                           + J
C
                           NKJ = IMATIJ(ISYMK,ISYMJ)
     &                         + NRHF(ISYMK)*(J - 1) + K
                           NDL = IT1AM(ISYMD,ISYML)
     &                         + NVIR(ISYMD)*(L - 1) + D
C
                           KOFFDLKJ = IOFF2(ISYMDL)
     &                              + NT1AM(ISYMDL)*(NKJ - 1) 
     &                              + NDL
C
                           XIINT(KOFFDLKJ) = XOINT(KOFFLJKD)
                           XJINT(KOFFDLKJ) = XOINT(KOFFJLKD)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END
C
C  /* Deck ccho_rdinmnf2_old */
      SUBROUTINE CCHO_RDINMNF2_OLD(LUNIT,FILNAM,XMINT,XNINT,
     &                             IA1,NUMIA,ISYMA)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C 
C     Read M and N integrals: 
C           I(bc,j;a) = (bj|ca) = X(bj,c;a) 
C           J(bc,j;a) = (cj|ab) = X(cj,b;a) 
C     X represents the integrals as in CHO_VI1
C     
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C  
      CHARACTER*7 FILNAM   
      DIMENSION XMINT(*),XNINT(*) 
C     
      PARAMETER (XMTHREE = -3.0D0, TWO = 2.0D0, XMONE = -1.0D0)
      PARAMETER (XMHALF = -0.5D0, ONE = 1.0D0)
C
      ISYMBJC = ISYMA
      ISYMBCJ = ISYMA
C
      IOFF = ICKBD(ISYMBJC,ISYMA) + NCKATR(ISYMBJC)*(IA1 - 1) + 1
      LENGTH = NCKATR(ISYMBJC)*NUMIA
C
      IF (LENGTH .GT. 0) THEN
         CALL GETWA2(LUNIT,FILNAM,XNINT,IOFF,LENGTH)
      ENDIF
C
      DO ISYMC = 1,NSYM
         ISYMBJ = MULD2H(ISYMBJC,ISYMC)
         DO ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMBJ,ISYMJ)
            ISYMBC = MULD2H(ISYMB,ISYMC)
C
            DO A = 1,NUMIA
               DO C = 1,NVIR(ISYMC)
                  DO J = 1,NRHF(ISYMJ)
                     DO B = 1,NVIR(ISYMB)
C
                        NBJCA = NCKATR(ISYMBJC)*(A - 1)
     &                        + ICKATR(ISYMBJ,ISYMC)
     &                        + NT1AM(ISYMBJ)*(C - 1)
     &                        + IT1AM(ISYMB,ISYMJ)
     &                        + NVIR(ISYMB)*(J - 1) + B
                        NBCJA = NCKASR(ISYMBCJ)*(A - 1)
     &                        + ICKASR(ISYMBC,ISYMJ)
     &                        + NMATAB(ISYMBC)*(J - 1)
     &                        + IMATAB(ISYMB,ISYMC)
     &                        + NVIR(ISYMB)*(C - 1) + B
C
                        XMINT(NBCJA) = XNINT(NBJCA)
C
                     ENDDO
                  ENDDO
               ENDDO               
            ENDDO
C
         ENDDO
      ENDDO
C
      DO ISYMC = 1,NSYM
         ISYMBJ = MULD2H(ISYMBJC,ISYMC)
         DO ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMBJ,ISYMJ)
            ISYMBC = MULD2H(ISYMB,ISYMC)
C
            DO A = 1,NUMIA
               DO C = 1,NVIR(ISYMC)
                  DO J = 1,NRHF(ISYMJ)
                     DO B = 1,NVIR(ISYMB)
C
                        NBCJA = NCKASR(ISYMBCJ)*(A - 1)
     &                        + ICKASR(ISYMBC,ISYMJ)
     &                        + NMATAB(ISYMBC)*(J - 1)
     &                        + IMATAB(ISYMB,ISYMC)
     &                        + NVIR(ISYMB)*(C - 1) + B
                        NCBJA = NCKASR(ISYMBCJ)*(A - 1)
     &                        + ICKASR(ISYMBC,ISYMJ)
     &                        + NMATAB(ISYMBC)*(J - 1)
     &                        + IMATAB(ISYMC,ISYMB)
     &                        + NVIR(ISYMC)*(B - 1) + C
C
                        XNINT(NBCJA) = XMINT(NCBJA)
C
                     ENDDO
                  ENDDO
               ENDDO               
            ENDDO
C
         ENDDO
      ENDDO
C
      NTOT = NCKASR(ISYMBCJ)*NUMIA
*     l = ntot
*     print *,ntot
*     print *,'XMINT norm inside rdint: ',dnrm2(l,xmint,1)
*     print *,'XNINT norm inside rdint: ',dnrm2(l,xnint,1)
C
      CALL DSCAL(NTOT,TWO,XMINT,1)
*     print *,'XMINT norm by two: ',dnrm2(l,xmint,1)
      CALL DAXPY(NTOT,XMONE,XNINT,1,XMINT,1)
*     print *,'XMINT norm       : ',dnrm2(l,xmint,1)
C
      CALL DSCAL(NTOT,XMTHREE,XNINT,1)
*     print *,'XNINT norm by three: ',dnrm2(l,xnint,1)
      CALL DAXPY(NTOT,ONE,XMINT,1,XNINT,1)
*     print *,'XNINT norm daxpy: ',dnrm2(l,xnint,1)
      CALL DSCAL(NTOT,XMHALF,XNINT,1)
*     print *,'XNINT norm       : ',dnrm2(l,xnint,1)
C
      RETURN
      END
C
C  /* Deck ccho_t2f2_old */
      SUBROUTINE CCHO_T2F2_OLD(T2VO,T21,T22,IA1,NUMIA,ISYMA)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C
C     Extract amplitudes: 
C                 T21(dl,i;a) = t(ad,li) 
C                 T22(dl,i;a) = t(ad,il) 
C     
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DIMENSION T2VO(*),T21(*),T22(*)
*     do i = 1,NCKI(ISYMA)*NUMIA
*        t21(i) = 888888888888D8
*        t22(i) = 888888888888D8
*     enddo
*     print *,'NDLIA must go from 1 to ',NCKI(ISYMA)*NUMIA
C
      ISYMDLI = ISYMA
C     
      DO ISYMI = 1,NSYM
C
         ISYMDL = MULD2H(ISYMDLI,ISYMI)
C
         DO ISYML = 1,NSYM
C
            ISYMLI = MULD2H(ISYMI,ISYML)
            ISYMAD = ISYMLI
            ISYMD  = MULD2H(ISYMAD,ISYMA)
C
            DO I = 1,NRHF(ISYMI)
               DO L = 1,NRHF(ISYML)
                  DO D = 1,NVIR(ISYMD)
                     DO A = 1,NUMIA
C
                        IA = IA1 + A - 1
C
                        NAD = IMATAB(ISYMA,ISYMD)
     &                      + NVIR(ISYMA)*(D - 1) + IA
                        NLI = IMATIJ(ISYML,ISYMI)
     &                      + NRHF(ISYML)*(I - 1) + L
                        NIL = IMATIJ(ISYMI,ISYML)
     &                      + NRHF(ISYMI)*(L - 1) + I
C
                        NADLI = IT2VO(ISYMAD,ISYMLI)
     &                        + NMATAB(ISYMAD)*(NLI - 1)
     &                        + NAD
                        NADIL = IT2VO(ISYMAD,ISYMLI)
     &                        + NMATAB(ISYMAD)*(NIL - 1)
     &                        + NAD
C                       NDLIA = IT2SP(ISYMDLI,ISYMA)
C    &                        + NCKI(ISYMDLI)*(A - 1)
                        NDLIA = NCKI(ISYMDLI)*(A - 1)
     &                        + ICKI(ISYMDL,ISYMI)
     &                        + NT1AM(ISYMDL)*(I - 1)
     &                        + IT1AM(ISYMD,ISYML)
     &                        + NVIR(ISYMD)*(L - 1) + D
C
                        T21(NDLIA) = T2VO(NADLI)
*      print *,'T21(',NDLIA,') = T2VO(',NADLI,') = ',T2VO(NADLI)
                        T22(NDLIA) = T2VO(NADIL)
*      print *,'T22(',NDLIA,') = T2VO(',NADIL,') = ',T2VO(NADIL)
C
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
*     print *,'T21 inside CCHO_T2F2:' 
*     print *,(t21(i),i=1,NCKI(ISYMA)*NUMIA)
*     print *,'T22 inside CCHO_T2F2:' 
*     print *,(t22(i),i=1,NCKI(ISYMA)*NUMIA)
      END
C    
C  /* Deck ccho_fscmn_olf */
      SUBROUTINE CCHO_FSCMN_OLD(XMINT,XNINT,CHOV,NUMIA,ISYMA)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C      
C     Scale M and N integrals with cholesky update
C           M(bc,j;a) = 2(bj|ca) - (cj|ba) 
C           N(bc,j;a) = 2(cj|ba) - (bj|ca)
C      
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C     
      DIMENSION XMINT(*),XNINT(*),CHOV(*) 
C     
      ISYMBCJ = ISYMA
C
      DO A = 1,NUMIA
         DO NBCJ = 1,NCKASR(ISYMBCJ)
C
            NBCJA = NCKASR(ISYMBCJ)*(A - 1) + NBCJ
C
            XMINT(NBCJA) = CHOV(NBCJ)*XMINT(NBCJA)
            XNINT(NBCJA) = CHOV(NBCJ)*XNINT(NBCJA)
C
         ENDDO   
      ENDDO
      RETURN
      END
C    
C  /* Deck ccho_f2sct2o_old */
      SUBROUTINE CCHO_F2SCT2O_OLD(T21,T22,CHOO,NUMIA,ISYMA) 
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     
C     Scale T21 and T22 with cholesky update 
C           T21(dl,i;a) --> d(dl,i) * t(ad,li) 
C           T22(dl,i;a) --> d(dl,i) * t(ad,il)
C     
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C     
      DIMENSION T21(*),T22(*),CHOO(*) 
C     
      ISYMDLI = ISYMA
C
      DO A = 1,NUMIA
         DO NDLI = 1,NCKI(ISYMDLI)
C
            NDLIA = NCKI(ISYMDLI)*(A - 1) + NDLI
C
            T21(NDLIA) = CHOO(NDLI)*T21(NDLIA)
            T22(NDLIA) = CHOO(NDLI)*T22(NDLIA)
C
         ENDDO
      ENDDO
C
      RETURN
      END
C
C  /* Deck ccho_dechof2_old */
      SUBROUTINE CCHO_DECHOF2_OLD(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                      ISYDLI,ISYBCJ)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez Oct 2002.
C
C     Construct vector for Cholesky decomposition:
C
C     For the occupied part, OCCHO(dli), the update vector
C     For the virtual  part, VICHO(bcj), the update vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
*     do i = 1,ncki(isymdli)
*        occho(i) = 8888888888888D28
*     enddo
*     do i = 1,nckatr(isymbcj)
*        vicho(i) = 8888888888888D28
*     enddo
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C-------------------
C     Occupied part.
C-------------------
C
      IF (ICHO .EQ. 1) THEN
         IND = 0
         DO ISYMI = 1,NSYM
            ISYMDL = MULD2H(ISYMI,ISYDLI)
            DO I = 1,NRHF(ISYMI)
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  DO L = 1,NRHF(ISYML)
                     DO D = 1,NVIR(ISYMD)
C
                        KOFFI = IRHF(ISYMI) + I
                        KOFFL = IRHF(ISYML) + L
                        KOFFD = IVIR(ISYMD) + D
C
                        OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) = 
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
         NDIMOC = IND
*     print *,'Cholesky occ. part dim.',ndimoc         
      ELSE
         IND = 0
         DO ISYMI = 1,NSYM
            ISYMDL = MULD2H(ISYMI,ISYDLI)
            DO I = 1,NRHF(ISYMI)
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  DO L = 1,NRHF(ISYML)
                     DO D = 1,NVIR(ISYMD)
C
                        KOFFI = IRHF(ISYMI) + I
                        KOFFL = IRHF(ISYML) + L
                        KOFFD = IVIR(ISYMD) + D
C
                        OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
         NDIMOC = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMOC,FACTOR,OCCHO,1)
C
*     print *,'Cholesky occ. part dim.',ndimoc 
      ENDIF
C
C
C------------------
C     Virtual part.
C------------------
C
      IF (ICHO .EQ. 1) THEN
         IND = 0
         DO ISYMJ = 1,NSYM
            ISYMBC = MULD2H(ISYMJ,ISYBCJ)
            DO J = 1,NRHF(ISYMJ)
               DO ISYMC = 1,NSYM
                  ISYMB = MULD2H(ISYMC,ISYMBC)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
C
                        KOFFJ = IRHF(ISYMJ) + J
                        KOFFC = IVIR(ISYMC) + C
                        KOFFB = IVIR(ISYMB) + B
C
                        OME = FOCKD(KOFFB)+FOCKD(KOFFC)-FOCKD(KOFFJ)
C
                        IND = IND + 1
                        VICHO(IND) = 
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
         NDIMVI = IND
*     print *,'Cholesky vir. part dim.',ndimvi         
C
      ELSE
C
         IND = 0
         DO ISYMJ = 1,NSYM
            ISYMBC = MULD2H(ISYMJ,ISYBCJ)
            DO J = 1,NRHF(ISYMJ)
               DO ISYMC = 1,NSYM
                  ISYMB = MULD2H(ISYMC,ISYMBC)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
C
                        KOFFJ = IRHF(ISYMJ) + J
                        KOFFC = IVIR(ISYMC) + C
                        KOFFB = IVIR(ISYMB) + B
C
                        OME = FOCKD(KOFFB)+FOCKD(KOFFC)-FOCKD(KOFFJ)
C
                        IND = IND + 1
                        VICHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
         NDIMVI = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMVI,FACTOR,VICHO,1)
C         
*     print *,'Cholesky vir. part dim.',ndimvi         
      ENDIF
C
      RETURN
      END
      
C  /* Deck ccho_fterm2 */
      SUBROUTINE CCHO_FTERM2(XIINT,XJINT,XOINT,T2VO,FOCKD,NUMCHO,CHOELE,
     &                       WORK,LWORK,E4F2,IA1,ISYMA,NUMIA,
     &                       PRINT)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     Calculate (part of F-terms):
C
C        E4F2 = 2*2 Sum(a) Sum(ijk) [ P(ja,ki) * Q(jk,ia)
C                               - R(ja,ki) * S(jk,ia) ]
C
C     where
C        
C        P(ja,ki) =     Sum(bc) d(bjc) * M(bc,j;a) * t(bc,ki)
C        R(ja,ki) =     Sum(bc) d(bjc) * N(bc,j;a) * t(bc,ki) 
C        Q(jk,ia) =   - Sum(dl) d(dli) * T1(dl,i;a) * I(dl,kj)
C                 + 2 * Sum(dl) d(dli) * T2(dl,i;a) * I(dl,kj)
C                 -     Sum(dl) d(dli) * T2(dl,i;a) * J(dl,kj)
C        S(jk,ia) =     Sum(dl) d(dli) * T1(dl,i;a) * J(dl,kj)
C
C     and M(bc,j;a)  = 2(bj|ca) - (cj|ba),
C         N(bc,j;a)  = 2(cj|ba) - (bj|ca),
C         T1(dl,i;a) = t(ad,li) 
C         T2(dl,i;a) = t(ad,il) 
C         I(dl,kj)   = (jk|ld) 
C         J(dl,kj)   = (lk|jd)
C     and d(dli) and d(bjc) denotes the occupied and virtual parts of
C     the Cholesky decomposition of the orbital energy denominator.
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XIINT(*),XJINT(*),XOINT(*),T2VO(*),FOCKD(*),CHOELE(*)
      DIMENSION WORK(LWORK)
      DIMENSION IOFF2(8)
      LOGICAL   PRINT
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>      
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR  =  4.0D0)
*     call escribe(lupri,xiint,xjint,xoint,t2vo,fockd,nrhft,nvirt)
C
      TLAST = SECOND()
C
C
      IF (PRINT) THEN
         TIMT = SECOND()
         WRITE(LUPRI,'(6X,A,/,6X,A,/)')
     &   'Calculation of the F2 term:',
     &   '==========================='
      ENDIF
C
      ISYMT2 = 1
C     
      MXJKI = -1
      DO ISYMI = 1,NSYM
         DO ISYMK = 1,NSYM
            DO ISYMJ = 1,NSYM
               NJKI   = NRHF(ISYMJ)*NRHF(ISYMK)*NRHF(ISYMI)
               IF (NJKI .GT. MXJKI) THEN
                  MXJKI = NJKI
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C      
C--------------------------------------------------------
C     IOFF2(ISYMDL) --> offset in OP(dl,kj) (ISYMOP = 1).
C--------------------------------------------------------
C      
      ICOUNT = 0
      DO ISYM = 1,NSYM
         IOFF2(ISYM) = ICOUNT
         ICOUNT = ICOUNT + NT1AM(ISYM)*NMATIJ(ISYM)
      ENDDO
      ISYMBCJ = ISYMA
      ISYMDLI = ISYMA
C         
C---------------------------
C     Dynamic allocation.
C---------------------------
C         
      KIINT = 1
      KJINT = KIINT + NCKIJ(1)
      KPMAT = KJINT + NCKIJ(1)         
      KQMAT = KPMAT + MXJKI
      KRMAT = KQMAT + MXJKI
      KSMAT = KRMAT + MXJKI
      KEND1 = KSMAT + MXJKI
      LWRK1 = LWORK - KEND1 + 1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in CCHO_FTERM2')
      ENDIF
C         
C----------------------------------------
C        Extract I(dl,kj) and J(dl,kj).
C----------------------------------------
C
         CALL CCHO_RDINIJF2(XOINT,WORK(KIINT),WORK(KJINT),IOFF2)
C
C            
C-----------------------------------------------
C           Allocate space for Cholesky vectors.
C-----------------------------------------------
C            
            LENCHOV = NCKATR(ISYMBCJ)
            LENCHOO = NCKI(ISYMDLI)
C
            KCHOV = KEND1
            KCHOO = KCHOV + LENCHOV
            KEND2 = KCHOO + LENCHOO
            LWRK2 = LWORK - KEND2 + 1
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Not enough space in CC_CHOPTF2')
            ENDIF
C                    
C---------------------------------
C              Dynamic allocation.
C---------------------------------
C              
               KMINT = KEND2
               KNINT = KMINT  + NCKATR(ISYMBCJ)*NUMIA
               KT21  = KNINT  + NCKATR(ISYMBCJ)*NUMIA
               KT22  = KT21   + NCKI(ISYMDLI)*NUMIA
               KENDA = KT22   + NCKI(ISYMDLI)*NUMIA
               LWRKA = LWORK  - KENDA + 1
               IF (LWRKA .LT. 0) THEN
               CALL QUIT('Insufficient memory in CC_CHOPTF2 in BATCH A')
               ENDIF
C               
C---------------------------------------------------------
c              Read the integrals M(bc,j;a) and N(bc,j;a).
C              Extract the integrals M(bc,j;a) and N(bc,j;a).
C---------------------------------------------------------
C              
               CALL CCHO_RDINMNF2(XIINT,XJINT,WORK(KMINT),WORK(KNINT),
     &                            IA1,NUMIA,ISYMA)
C               
C------------------------------------------------
C              Extract T1(dl,i;a) and T2(dl,i;a).
C------------------------------------------------
C               
               CALL CCHO_T2F2(T2VO,WORK(KT21),WORK(KT22),
     &                        IA1,NUMIA,ISYMA)    
C              
               DO ICHO = 1,NUMCHO
C
                  FCOR = ZERO
C                  
C------------------------------------------------
C                 Calculate Cholesky information.
C                 d(dl,i) update vector
C                 d(bc,j) update vector
C                 (TODO  check if it is right)
C------------------------------------------------
C 
                  CALL CCHO_DECHOF2(FOCKD,CHOELE,NUMCHO,ICHO,
     &                              WORK(KCHOO),WORK(KCHOV),
     &                              ISYMDLI,ISYMBCJ)
C                  
C--------------------------------------------------------------
C                 Scale M and N integrals with cholesky update.
C--------------------------------------------------------------
C                  
                  CALL CCHO_FSCMN(WORK(KMINT),WORK(KNINT),WORK(KCHOV),
     &                            NUMIA,ISYMA)
C                  
C----------------------------------------------------------
C                 Scale T2 amplitudes with cholesky update.
C----------------------------------------------------------
C                  
                  CALL CCHO_F2SCT2O(WORK(KT21),WORK(KT22),WORK(KCHOO),
     &                              NUMIA,ISYMA) 
C                 
                  DO ISYMJ = 1,NSYM
C                 
                     IF ( NRHF(ISYMJ) .EQ. 0 ) GOTO 1234
C
                     ISYMBC  = MULD2H(ISYMJ,ISYMA)
                     ISYMKI  = ISYMBC
                     ISYMDLK = ISYMJ
C
                     DO ISYMK = 1,NSYM
C                 
                        IF ( NRHF(ISYMK) .EQ. 0 ) GOTO 1235
C
                        ISYMI   = MULD2H(ISYMKI,ISYMK)
                        ISYMDL  = MULD2H(ISYMI,ISYMA)
                        ISYMKJ  = ISYMDL                       
C
                        DO A = 1,NUMIA
C                        
C----------------------------------------------------------
C                          Calculate P and R intermediates.
C----------------------------------------------------------
C                        
                           KOFFT = IT2VO(ISYMBC,ISYMKI) 
     &                           + NMATAB(ISYMBC)*IMATIJ(ISYMK,ISYMI)
     &                           + 1                           
                           KOFFM = KMINT + NCKASR(ISYMBCJ)*(A - 1)
     &                           + ICKASR(ISYMBC,ISYMJ)
                           KOFFN = KNINT + NCKASR(ISYMBCJ)*(A - 1)
     &                           + ICKASR(ISYMBC,ISYMJ)
C
                           NBC   = NMATAB(ISYMBC)
                           NTOBC = MAX(1,NBC)
                           NJ    = NRHF(ISYMJ)
                           NTOJ  = MAX(1,NJ)
                           NKI   = NRHF(ISYMK)*NRHF(ISYMI)
C
                           CALL DGEMM('T','N',NJ,NKI,NBC,ONE,
     &                                WORK(KOFFM),NTOBC,
     &                                T2VO(KOFFT),NTOBC,
     &                                ZERO,WORK(KPMAT),NTOJ)
C
                           CALL DGEMM('T','N',NJ,NKI,NBC,ONE,
     &                                WORK(KOFFN),NTOBC,
     &                                T2VO(KOFFT),NTOBC,
     &                                ZERO,WORK(KRMAT),NTOJ)
C                           
C----------------------------------------------------------
C                          Calculate Q and S intermediates.
C----------------------------------------------------------
C                           
                           KOFFI   = KIINT + IOFF2(ISYMDL)
     &                             + NT1AM(ISYMDL)*IMATIJ(ISYMK,ISYMJ)
                           KOFFJ   = KJINT + IOFF2(ISYMDL)
     &                             + NT1AM(ISYMDL)*IMATIJ(ISYMK,ISYMJ)
                           KOFFT21 = KT21 + NCKI(ISYMDLI)*(A - 1) 
     &                             + ICKI(ISYMDL,ISYMI)
                           KOFFT22 = KT22 + NCKI(ISYMDLI)*(A - 1)
     &                             + ICKI(ISYMDL,ISYMI)
C
                           NDL   = NT1AM(ISYMDL)
                           NTODL = MAX(1,NDL)
                           NKJ   = NRHF(ISYMK)*NRHF(ISYMJ)
                           NI    = NRHF(ISYMI)
                           NTOI  = MAX(1,NI)
C
                           CALL DGEMM('T','N',NI,NKJ,NDL,XMONE,
     &                                WORK(KOFFT21),NTODL,
     &                                WORK(KOFFI),NTODL,
     &                                ZERO,WORK(KQMAT),NTOI)
C
                           CALL DGEMM('T','N',NI,NKJ,NDL,TWO,
     &                                WORK(KOFFT22),NTODL,
     &                                WORK(KOFFI),NTODL,
     &                                ONE,WORK(KQMAT),NTOI)
C
                           CALL DGEMM('T','N',NI,NKJ,NDL,XMONE,
     &                                WORK(KOFFT22),NTODL,
     &                                WORK(KOFFJ),NTODL,
     &                                ONE,WORK(KQMAT),NTOI)
C
                           CALL DGEMM('T','N',NI,NKJ,NDL,ONE,
     &                                WORK(KOFFT21),NTODL,
     &                                WORK(KOFFJ),NTODL,
     &                                ZERO,WORK(KSMAT),NTOI)
C                           
C--------------------------------------------------------
C                          Calculate energy contribution.
C--------------------------------------------------------
C                           
                           DO I = 1,NRHF(ISYMI)
                              DO K = 1,NRHF(ISYMK)
                                 DO J = 1,NRHF(ISYMJ)
C
                                    NJKI = NRHF(ISYMJ)*NRHF(ISYMK)*(I-1)
     &                                   + NRHF(ISYMJ)*(K - 1) + J
                                    NIKJ = NRHF(ISYMI)*NRHF(ISYMK)*(J-1)
     &                                   + NRHF(ISYMI)*(K - 1) + I
C
                                    KOFFP = KPMAT + NJKI - 1
                                    KOFFR = KRMAT + NJKI - 1
                                    KOFFQ = KQMAT + NIKJ - 1
                                    KOFFS = KSMAT + NIKJ - 1
C
                                    FCOR = FCOR
     &                                   + WORK(KOFFP)*WORK(KOFFQ)
     &                                   - WORK(KOFFR)*WORK(KOFFS)
C
                                 ENDDO
                              ENDDO
                           ENDDO
C
                        ENDDO
C                       End NUMIA loop
 1235                   CONTINUE
                     ENDDO
C                    End ISYMK loop
 1234                CONTINUE
                  ENDDO
C                 End ISYMJ loop
C
                  FCOR = TWO * FCOR
C
                  E4F2 = E4F2 + FCOR
C 
              ENERGF(ICHO) = ENERGF(ICHO) + FCOR
C 
C Decommented by Domenico
                    TNOW = SECOND()
                    DELTAT = TNOW - TLAST
                    TLAST = TNOW
                    SCNDSF(ICHO) = SCNDSF(ICHO) 
     &                          + DELTAT
C Decommented by Domenico
C 
C
c                 IF (PRINT) THEN
c                    WRITE(LUPRI,'(15X,A,I3,A,/,15X,A)')
c    &               'Status after Cholesky vector',ICHO,':',
c    &               '--------------------------------'
c                    IF (ABS(FCOR) .LT. THRCHO) THEN
c                       WRITE(LUPRI,'(15X,A)') 'F2 term converged'
c                    ELSE
c                       WRITE(LUPRI,'(15X,A)') 'F2 term not converged'
c                    ENDIF
c                    TIM = SECOND() - TIMT
c                    WRITE(LUPRI,'(15X,A,F10.2,A,/)')
c    &               'Accumulated F2-time: ',TIM,' seconds'
c                 ENDIF
C
                  IF (DABS(FCOR) .LT. THRCHO) GOTO 2
C
               ENDDO
C              End ICHO loop
    2          CONTINUE               
C
      RETURN
      END
C
C  /* Deck ccho_rdinijf2 */
      SUBROUTINE CCHO_RDINIJF2(XOINT,XIINT,XJINT,IOFF2) 
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     
C     Extract occupied integrals: 
C                    I(dl,kj) = (jk|ld) 
C                    J(dl,kj) = (lk|jd) 
C     from XOINT(lj,k;d) = (jk|ld) 
C     
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DIMENSION XOINT(*),XIINT(*),XJINT(*),IOFF2(*) 
C     
      DO ISYMD = 1,NSYM
         ISYMLJK = ISYMD
         DO ISYMK = 1,NSYM
            ISYMLJ = MULD2H(ISYMLJK,ISYMK)
            DO ISYMJ = 1,NSYM
               ISYML = MULD2H(ISYMLJ,ISYMJ)
               ISYMDLK = ISYMJ
               ISYMDL = MULD2H(ISYMDLK,ISYMK)
C
               DO D = 1,NVIR(ISYMD)
                  DO K = 1,NRHF(ISYMK)
                     DO J = 1,NRHF(ISYMJ)
                        DO L = 1,NRHF(ISYML)
C
                           KOFFLJKD = ISJIKA(ISYMLJK,ISYMD)
     &                           + NMAIJK(ISYMLJK)*(D - 1)
     &                           + IMAIJK(ISYMLJ,ISYMK)
     &                           + NMATIJ(ISYMLJ)*(K - 1)
     &                           + IMATIJ(ISYML,ISYMJ)
     &                           + NRHF(ISYML)*(J - 1)
     &                           + L
C
                           KOFFJLKD = ISJIKA(ISYMLJK,ISYMD)
     &                           + NMAIJK(ISYMLJK)*(D - 1)
     &                           + IMAIJK(ISYMLJ,ISYMK)
     &                           + NMATIJ(ISYMLJ)*(K - 1)
     &                           + IMATIJ(ISYMJ,ISYML)
     &                           + NRHF(ISYMJ)*(L - 1)
     &                           + J
C
                           NKJ = IMATIJ(ISYMK,ISYMJ)
     &                         + NRHF(ISYMK)*(J - 1) + K
                           NDL = IT1AM(ISYMD,ISYML)
     &                         + NVIR(ISYMD)*(L - 1) + D
C
                           KOFFDLKJ = IOFF2(ISYMDL)
     &                              + NT1AM(ISYMDL)*(NKJ - 1) 
     &                              + NDL
C
                           XIINT(KOFFDLKJ) = XOINT(KOFFLJKD)
                           XJINT(KOFFDLKJ) = XOINT(KOFFJLKD)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END
C
C  /* Deck ccho_rdinmnf2 */
      SUBROUTINE CCHO_RDINMNF2(XIINT,XJINT,XMINT,XNINT,IA1,NUMIA,ISYMA)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C 
C     Extract M and N integrals from 
C        I(em,d;#b) = (#bm|ed)      
C        J(em,d;#b) = (em|#bd)
C     M(bc,j;#a) = 2 J(bj,c;#a) - J(cj,b;#a) 
C     N(bj,j;#a) = 2 J(cj,b;#a) - J(bj,c;#a)
C     
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C  
      DIMENSION XIINT(*),XJINT(*),XMINT(*),XNINT(*) 
C     
      PARAMETER ( TWO = 2.0D0 )
C
      ISYMBCJ = ISYMA
C
      DO ISYMC = 1,NSYM
         ISYMBJ  = MULD2H(ISYMBCJ,ISYMC)
         DO ISYMJ = 1,NSYM
            ISYMB  = MULD2H(ISYMBJ,ISYMJ)
            ISYMBC = MULD2H(ISYMB,ISYMC) 
            ISYMCJ  = MULD2H(ISYMC,ISYMJ)
C           
            DO C = 1,NVIR(ISYMC)
               DO A = 1,NUMIA
                  DO J = 1,NRHF(ISYMJ)
                     DO B = 1,NVIR(ISYMB)
C
                        NBJCA = NCKATR(ISYMBCJ)*(A - 1)
     &                        + ICKATR(ISYMBJ,ISYMC)
     &                        + NT1AM(ISYMBJ)*(C - 1)
     &                        + IT1AM(ISYMB,ISYMJ)
     &                        + NVIR(ISYMB)*(J - 1)
     &                        + B
                        NCJBA = NCKATR(ISYMBCJ)*(A - 1)
     &                        + ICKATR(ISYMCJ,ISYMB)
     &                        + NT1AM(ISYMCJ)*(B - 1)
     &                        + IT1AM(ISYMC,ISYMJ)
     &                        + NVIR(ISYMC)*(J - 1)
     &                        + C
                        NBCJA = NCKASR(ISYMBCJ)*(A - 1)
     &                        + ICKASR(ISYMBC,ISYMJ)
     &                        + NMATAB(ISYMBC)*(J - 1)
     &                        + IMATAB(ISYMB,ISYMC)
     &                        + NVIR(ISYMB)*(C - 1) + B
C
                        XMINT(NBCJA) = TWO*XJINT(NBJCA) - XJINT(NCJBA)
                        XNINT(NBCJA) = TWO*XJINT(NCJBA) - XJINT(NBJCA)
C
                     ENDDO
                  ENDDO
               ENDDO
          ENDDO
C           
         ENDDO
      ENDDO
C
      RETURN
      END
C
C  /* Deck ccho_t2f2 */
      SUBROUTINE CCHO_T2F2(T2VO,T21,T22,IA1,NUMIA,ISYMA)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C
C     Extract amplitudes: 
C                 T21(dl,i;a) = t(ad,li) 
C                 T22(dl,i;a) = t(ad,il) 
C     
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DIMENSION T2VO(*),T21(*),T22(*)
C
      ISYMDLI = ISYMA
C     
      DO ISYMI = 1,NSYM
C
         ISYMDL = MULD2H(ISYMDLI,ISYMI)
C
         DO ISYML = 1,NSYM
C
            ISYMLI = MULD2H(ISYMI,ISYML)
            ISYMAD = ISYMLI
            ISYMD  = MULD2H(ISYMAD,ISYMA)
C
            DO I = 1,NRHF(ISYMI)
               DO L = 1,NRHF(ISYML)
                  DO D = 1,NVIR(ISYMD)
                     DO A = 1,NUMIA
C
                        IA = IA1 + A - 1
C
                        NAD = IMATAB(ISYMA,ISYMD)
     &                      + NVIR(ISYMA)*(D - 1) + IA
                        NLI = IMATIJ(ISYML,ISYMI)
     &                      + NRHF(ISYML)*(I - 1) + L
                        NIL = IMATIJ(ISYMI,ISYML)
     &                      + NRHF(ISYMI)*(L - 1) + I
C
                        NADLI = IT2VO(ISYMAD,ISYMLI)
     &                        + NMATAB(ISYMAD)*(NLI - 1)
     &                        + NAD
                        NADIL = IT2VO(ISYMAD,ISYMLI)
     &                        + NMATAB(ISYMAD)*(NIL - 1)
     &                        + NAD
C                       NDLIA = IT2SP(ISYMDLI,ISYMA)
C    &                        + NCKI(ISYMDLI)*(A - 1)
                        NDLIA = NCKI(ISYMDLI)*(A - 1)
     &                        + ICKI(ISYMDL,ISYMI)
     &                        + NT1AM(ISYMDL)*(I - 1)
     &                        + IT1AM(ISYMD,ISYML)
     &                        + NVIR(ISYMD)*(L - 1) + D
C
                        T21(NDLIA) = T2VO(NADLI)
                        T22(NDLIA) = T2VO(NADIL)
C
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
      END
C    
C  /* Deck ccho_fscmn */
      SUBROUTINE CCHO_FSCMN(XMINT,XNINT,CHOV,NUMIA,ISYMA)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C      
C     Scale M and N integrals with cholesky update
C           M(bc,j;a) = 2(bj|ca) - (cj|ba) 
C           N(bc,j;a) = 2(cj|ba) - (bj|ca)
C      
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C     
      DIMENSION XMINT(*),XNINT(*),CHOV(*) 
C     
      ISYMBCJ = ISYMA
C
      DO A = 1,NUMIA
         DO NBCJ = 1,NCKASR(ISYMBCJ)
C
            NBCJA = NCKASR(ISYMBCJ)*(A - 1) + NBCJ
C
            XMINT(NBCJA) = CHOV(NBCJ)*XMINT(NBCJA)
            XNINT(NBCJA) = CHOV(NBCJ)*XNINT(NBCJA)
C
         ENDDO   
      ENDDO
      RETURN
      END
C    
C  /* Deck ccho_f2sct2o */
      SUBROUTINE CCHO_F2SCT2O(T21,T22,CHOO,NUMIA,ISYMA) 
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez June 2002.
C
C     
C     Scale T21 and T22 with cholesky update 
C           T21(dl,i;a) --> d(dl,i) * t(ad,li) 
C           T22(dl,i;a) --> d(dl,i) * t(ad,il)
C     
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C     
      DIMENSION T21(*),T22(*),CHOO(*) 
C     
      ISYMDLI = ISYMA
C
      DO A = 1,NUMIA
         DO NDLI = 1,NCKI(ISYMDLI)
C
            NDLIA = NCKI(ISYMDLI)*(A - 1) + NDLI
C
            T21(NDLIA) = CHOO(NDLI)*T21(NDLIA)
            T22(NDLIA) = CHOO(NDLI)*T22(NDLIA)
C
         ENDDO
      ENDDO
C
      RETURN
      END
C
C  /* Deck ccho_dechof2 */
      SUBROUTINE CCHO_DECHOF2(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                      ISYDLI,ISYBCJ)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez Oct 2002.
C
C     Construct vector for Cholesky decomposition:
C
C     For the occupied part, OCCHO(dli), the update vector
C     For the virtual  part, VICHO(bcj), the update vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C-------------------
C     Occupied part.
C-------------------
C
      IF (ICHO .EQ. 1) THEN
         IND = 0
         DO ISYMI = 1,NSYM
            ISYMDL = MULD2H(ISYMI,ISYDLI)
            DO I = 1,NRHF(ISYMI)
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  DO L = 1,NRHF(ISYML)
                     DO D = 1,NVIR(ISYMD)
C
                        KOFFI = IRHF(ISYMI) + I
                        KOFFL = IRHF(ISYML) + L
                        KOFFD = IVIR(ISYMD) + D
C
                        OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) = 
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
         NDIMOC = IND
      ELSE
         IND = 0
         DO ISYMI = 1,NSYM
            ISYMDL = MULD2H(ISYMI,ISYDLI)
            DO I = 1,NRHF(ISYMI)
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  DO L = 1,NRHF(ISYML)
                     DO D = 1,NVIR(ISYMD)
C
                        KOFFI = IRHF(ISYMI) + I
                        KOFFL = IRHF(ISYML) + L
                        KOFFD = IVIR(ISYMD) + D
C
                        OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
         NDIMOC = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMOC,FACTOR,OCCHO,1)
C
      ENDIF
C
C
C------------------
C     Virtual part.
C------------------
C
      IF (ICHO .EQ. 1) THEN
         IND = 0
         DO ISYMJ = 1,NSYM
            ISYMBC = MULD2H(ISYMJ,ISYBCJ)
            DO J = 1,NRHF(ISYMJ)
               DO ISYMC = 1,NSYM
                  ISYMB = MULD2H(ISYMC,ISYMBC)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
C
                        KOFFJ = IRHF(ISYMJ) + J
                        KOFFC = IVIR(ISYMC) + C
                        KOFFB = IVIR(ISYMB) + B
C
                        OME = FOCKD(KOFFB)+FOCKD(KOFFC)-FOCKD(KOFFJ)
C
                        IND = IND + 1
                        VICHO(IND) = 
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
         NDIMVI = IND
C
      ELSE
C
         IND = 0
         DO ISYMJ = 1,NSYM
            ISYMBC = MULD2H(ISYMJ,ISYBCJ)
            DO J = 1,NRHF(ISYMJ)
               DO ISYMC = 1,NSYM
                  ISYMB = MULD2H(ISYMC,ISYMBC)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
C
                        KOFFJ = IRHF(ISYMJ) + J
                        KOFFC = IVIR(ISYMC) + C
                        KOFFB = IVIR(ISYMB) + B
C
                        OME = FOCKD(KOFFB)+FOCKD(KOFFC)-FOCKD(KOFFJ)
C
                        IND = IND + 1
                        VICHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
         NDIMVI = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMVI,FACTOR,VICHO,1)
C         
      ENDIF
C
      RETURN
      END
C  /* Deck cc_choc1term */
      SUBROUTINE CC_CHOC1TERM(XOINT,T2VO,FOCKD,NUMCHO,CHOELE,WORK,LWORK,
     &                       E4C1,FBATCH,IPRINT,NONI,NONJ)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez April 2003.
C
C        E4C1 = Sum(#i#jkl) [ V1(#i#jkl) * ( A(l#jk#i) + B(l#jk#i) ) 
C                           + V2(#i#jlk) * D(l#jk#i) ]
C     
C     where
C        
C        V1(#i#jkl) =   Sum(ab) d(ab#i)*T21(ab,#i#j) * T2VO(ab,kl)
C        V2(#i#jlk) =   Sum(ab) d(ab#i)*T21(ab,#i#j) * T2VO(ab,lk)
C        A(l#jk#i) = - Sum(cm) d(cml)*I(cm,l#j) * L(cm,k#i) 
C        B(l#jk#i) =   Sum(cm) d(cml)*J(cm,l#j) * C(cm,k#i)
C        D(l#jk#i) =   Sum(cm) d(cml)*J(cm,l#j) * M(cm,k#i)
C     
C        T21(ab,ij)  = s(ab,ij) = 2*t(ai,bj) - t(aj,bi)
C                    = 2*T2VO(ab,ij) - T2VO(ab,ji)
C        T2VO(ab,kl) = t(ak,bl)
C        C(cm,k#i)   = (ik|mc) 
C                    = XOINT(mi,k;c)
C        L(cm,k#i)   = 2*(ik|mc) - (ic|mk)
C                    = 2*XOINT(mi,k;c) - XOINT(im,k;c)
C        M(cm,k#i)   = (mk|ic) = XOINT(im,k;c)
C        I(cm,l#j)   = (cm|lj) = XOINT(ml,j;c)
C        J(cm,l#j)   = (cl|mj) = XOINT(lm,j;c)
C     
C        XOINT(ml,j;c) = (cm|lj)
C     
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XOINT(*),T2VO(*),FOCKD(*),CHOELE(*)
      DIMENSION WORK(LWORK)
      LOGICAL   FBATCH,IPRINT
C
      DIMENSION IOFFI(8),IOFFJ(8) 
*     dimension ilabel(25)      
C     
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
C 
      TLAST = SECOND()
C 
*     icounter=0
*     do jj = 1,5
*        do ii = 1,5
*           icounter=icounter+1
*           ilabel(icounter)=jj+ii*10
*        enddo
*     enddo 
      IF (IPRINT) THEN
         TIMT = SECOND()
         WRITE(LUPRI,'(A,/,A,/)')
     &   'Calculation of the C1 term:',
     &   '==========================='
      ENDIF
*     
*     WRITE(LUPRI,*) 'T2VO(AB,IJ) FOR AB=2'
*     n = nrhf(1)
*     nab = nvir(1)*nvir(1)
*     do j=1,n
*        do i=1,n
*           WRITE(LUPRI,*) 'T2VO(AB',I,',',J,') =',xoint( nab*(n*(j-1)+i-1)+2 )
*        enddo
*     enddo
*     WRITE(LUPRI,*) 'T2VO:'
*     ntot=nrhft**2*nvirt**2
*     print '(8F10.7)',(T2VO(ii),ii=1,ntot)
* 
*     xt2norm = dnrm2(nt2sq(1),t2vo,1) 
*     xointnorm = dnrm2(ntraoc(1),xoint,1) 
*     WRITE(LUPRI,*) 'Initial T2VO norm: ',xt2norm
*     WRITE(LUPRI,*) 'Initial XOINT norm: ',xointnorm
*     WRITE(LUPRI,*) 'NUMCHO = ',numcho
*     print '(8F20.15)',(choele(ii),ii=1,numcho)
*     WRITE(LUPRI,*) 'First elements of FOCKD:'
*     print '(8F20.15)',(fockd(ii),ii=1,20)
* 
C     MXAB is the biggest AB symmetry matrix 
      MXAB = -1
      NAB = 0
      DO ISYMAB = 1,NSYM
         NAB = NMATAB(ISYMAB)
         IF (MXAB .LT. NAB) MXAB = NAB
      ENDDO
*     WRITE(LUPRI,*) 'MXAB = ',mxab      
C     
      DO ISYMI = 1,NSYM
*        WRITE(LUPRI,*)  'ISYMI = ',isymi
*        WRITE(LUPRI,*)  '----------------------------------'
C     
         IF ( NRHF(ISYMI) .EQ. 0 ) GOTO 1234
C      
         ISYMABJ = ISYMI
         ISYMCMK = ISYMI
C         
         LWORKI = LWORK/3
C
C         
C-------------------------------
C        Batch over the I index.
C-------------------------------
C         
         LENI = 3*NCKI(ISYMCMK) + MXAB
         IF (FBATCH) THEN
            NEFI   = MIN(NONI,NRHF(ISYMI))
            LEFF   = NEFI*LENI + 1
            LWORKI = MIN(LWORKI,LEFF)
         ENDIF
         NUMI = MIN(NRHF(ISYMI),LWORKI/LENI)
         IF (NUMI .EQ. 0) THEN
            WRITE(LUPRI,*) 'NUMI .EQ. 0 IN CC_CHOPTC1'
            CALL QUIT('NOT ENOUGH SPACE IN CC_CHOPTC1 !!')
         ENDIF
C         
         NBATI = (NRHF(ISYMI)-1)/NUMI + 1
C
C
         IF (IPRINT) THEN
            WRITE(LUPRI,'(3X,A,I1,A,/,3X,A)')
     &      'Batch over I, symmetry ',ISYMI,':',
     &      '-------------------------'
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10)')
     &      'Minimum work space required   : ',LENI,
     &      'Work space available for batch: ',LWORKI
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10,/)')
     &      'Number of occupied orbitals   : ',NRHF(ISYMI),
     &      'Required number of J-batches  : ',NBATI
         ENDIF
C        
         II2=0
         DO IBATI = 1,NBATI
*           WRITE(LUPRI,*)  'IBATI = ',ibati
*           xcount = xcount + one
            II1 = II2 + 1
            II2 = II2 + NUMI
            IF (II2 .GT. NRHF(ISYMI)) II2 = NRHF(ISYMI)
            NUMII = II2 - II1 + 1
C
            IF (IPRINT) THEN
               WRITE(LUPRI,'(6X,A,I10,A,/,6X,A)')
     &         'I-batch number ',IBATI,':',
     &         '--------------------------'
               WRITE(LUPRI,'(6X,A,I10,1X,I10,/)')
     &         'First and last I: ',II1,II2
            ENDIF
C              
C           IOFFI(8) offset
C           ---------------
            ICOUNT = 0
            DO ISYMK = 1,NSYM
               ISYMCM = MULD2H(ISYMI,ISYMK)
               IOFFI(ISYMK) = ICOUNT
               ICOUNT = ICOUNT + NT1AM(ISYMCM)*NRHF(ISYMK)*NUMII
            ENDDO
*           WRITE(LUPRI,*) 'IOFFI : ',ioffi
C            
C------------------------------
C           Dynamic allocation.
C------------------------------
C
            KCMAT = 1
            KLINT = KCMAT + NCKI(ISYMCMK)*NUMII
            KMINT = KLINT + NCKI(ISYMCMK)*NUMII
            KENDI = KMINT + NCKI(ISYMCMK)*NUMII
            LWRKI = LWORKI - KENDI + 1
            IF (LWRKI .LT. 0) THEN
            CALL QUIT('Allocation bug in CC_CHOC1TERM in the I Cbatch')
            ENDIF
C             
C------------------------------------------------------------------
C           Extract integrals: 
C              C(cm,ki) =  (ik|mc) 
C              L(cm,ki) = 2(ik|mc) - (ic|mk)
C              M(cm,ki) =  (mk|ic)
C           ordered for a i batch with a fixed ISYMI using IOFFI
C------------------------------------------------------------------
C            
            CALL CCHO_RDINTC11(XOINT,WORK(KCMAT),WORK(KLINT),
     &                        WORK(KMINT),II1,NUMII,ISYMI,IOFFI)
*           xointnorm = dnrm2(ntraoc(1),xoint,1)
*           ntot = ncki(isymcmk)*numii
*           xcnorm    = dsqrt(xcnorm**2
*    &                +       dnrm2(ntot,work(kcmat),1)**2) 
*           xlnorm    = dsqrt(xlnorm**2
*    &                +       dnrm2(ntot,work(klint),1)**2) 
*           xmnorm    = dsqrt(xmnorm**2
*    &                +       dnrm2(ntot,work(kmint),1)**2) 
*           WRITE(LUPRI,*) 'After CCHO_RDINTC11:' 
*           WRITE(LUPRI,*) 'XOINT norm:',xointnorm
*           WRITE(LUPRI,*) 'C norm:',xcnorm
*           WRITE(LUPRI,*) 'L norm:',xlnorm
*           WRITE(LUPRI,*) 'M norm:',xmnorm
*           call dcopy(ntraoc(1),work(kcmat),1,work(kendi),1)
*           call daxpy(ntraoc(1),one,work(klint),1,work(kendi),1) 
*           WRITE(LUPRI,*) 'C + L norm / 3:',
*    &              dnrm2(ntraoc(1),work(kendi),1)/3.0d0
C           
*           xiintnrm = zero
*           xjintnrm = zero
*           xiintsclnrm = zero
*           xjintsclnrm = zero
* 
            DO ISYMJ = 1,NSYM
*        WRITE(LUPRI,*)  'ISYMJ = ',isymj
*        WRITE(LUPRI,*)  '----------------------------------'
*
*        if(nrhf(isymj).eq.0) then
*              KCHOV = KENDI
*              KCHOO = KCHOV + NMATAB(ISYMAB)*NUMII
*              KENDO = KCHOO + NCKI(ISYMCML)
*                 icho = 1
*                    CALL CCHO_DECHOC1(FOCKD,CHOELE,NUMCHO,ICHO,
*    &                                 WORK(KCHOO),WORK(KCHOV),
*    &                                 ISYMCML,ISYMABI,
*    &                                 II1,NUMII,ISYMI)
*                 ntoto = ncki(isymcml) 
*                 ntotv = nmatab(isymab)*numii
*                 xchoonrm = dnrm2(ntoto,work(kchoo),1) 
*                 xchovnrm = dnrm2(ntotv,work(kchov),1)
*                 WRITE(LUPRI,*) 'Cholesky vector #',icho
*                 WRITE(LUPRI,*) 'd(ab,#i) update norm:',xchovnrm
*     print '(5F20.15)',(work(ii),ii=kchov,kchov+ntotv-1)
*                 WRITE(LUPRI,*) 'd(cm,l) update norm:',xchoonrm
*           endif
C           
               IF(NRHF(ISYMJ) .EQ. 0) GOTO 1235
C              
               ISYMAB = MULD2H(ISYMI,ISYMJ)
C               
               IF(NMATAB(ISYMAB) .EQ. 0) GOTO 1235
C               
               ISYMCML = ISYMJ
               ISYMABI = MULD2H(ISYMAB,ISYMI)
               ISYMKL = MULD2H(ISYMI,ISYMJ)
C              
C              Dynamic allocation for d(cml) and d(ab#i)
C              -----------------------------------------
               KCHOV = KENDI
               KCHOO = KCHOV + NMATAB(ISYMAB)*NUMII
               KENDO = KCHOO + NCKI(ISYMCML)
               LWRKO = LWORK - KENDO + 1
               IF (LWRKO .LT. 0) THEN
                CALL QUIT('Insufficient memory in CC_CHOC1TERM in C#J1')
               ENDIF
C              
C              MXKL is the biggest block inside the KL matrix with
C              ISYMKL symmetry               
C              ---------------------------------------------------
               MXKL = -1
               DO ISYMK = 1,NSYM
                  ISYML = MULD2H(ISYMKL,ISYMK)
                  NKL = NRHF(ISYMK)*NRHF(ISYML)
                  IF (MXKL .LT. NKL) MXKL = NKL
               ENDDO
C               
C---------------------------
C              Batch over J.
C---------------------------
C              
               LENJ = NMATAB(ISYMAB)*NUMII + 2*NCKI(ISYMCML)
     &              + 4*MXKL*NUMII
C              
               LWORKJ = LWORK - KENDO + 1
               IF (FBATCH) THEN
                  NEFJ = MIN(NRHF(ISYMJ),NONJ)
                  LEFF = NEFJ*LENJ + 1
                  LWORKJ = MIN(LWORKJ,LEFF)
               ENDIF
               NUMJ = MIN(NRHF(ISYMJ),LWORKJ/LENJ)
               IF (NUMJ .EQ. 0) THEN
                  WRITE(LUPRI,*) 'NUMJ .EQ. 0 in CC_CHOPTC1'
                  CALL QUIT('Not enough space in CC_CHOPTC1')
               ENDIF
C                 
               NBATJ = (NRHF(ISYMJ)-1)/NUMJ + 1
C
               IF (IPRINT) THEN
                  WRITE(LUPRI,'(9X,A,I1,A,/,9X,A)')
     &            'Batch over J, symmetry ',ISYMJ,':',
     &            '-------------------------'
                  WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10)')
     &            'Minimum work space required   : ',LENJ,
     &            'Work space available for batch: ',LWORKJ
                  WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10,/)')
     &            'Number of occupied orbitals    : ',NRHF(ISYMJ),
     &            'Required number of J-batches  : ',NBATJ
               ENDIF
C              
               IJ2 = 0
               DO IBATJ = 1,NBATJ
                  IJ1 = IJ2 + 1
                  IJ2 = IJ2 + NUMJ
                  IF(IJ2 .GT. NRHF(ISYMJ)) IJ2 = NRHF(ISYMJ)
                  NUMIJ = IJ2 - IJ1 + 1
C
                  IF (IPRINT) THEN
                     WRITE(LUPRI,'(12X,A,I10,A,/,12X,A)')
     &               'J-batch number ',IBATJ,':',
     &               '--------------------------'
                     WRITE(LUPRI,'(12X,A,I10,1X,I10,/)')
     &               'First and last J: ',IJ1,IJ2
                  ENDIF
C                  
C                 IOFFJ(8) offset
C                 ---------------
                  ICOUNT = 0
                  DO ISYML = 1,NSYM
                     ISYMCM = MULD2H(ISYML,ISYMJ)
                     IOFFJ(ISYML) = ICOUNT
                     ICOUNT = ICOUNT + NT1AM(ISYMCM)*NRHF(ISYML)*NUMIJ
                  ENDDO
*                 WRITE(LUPRI,*) 'IOFFJ : ',ioffj
C                  
C------------------------------------
C                 Dynamic allocation.
C------------------------------------
C                  
                  KT21   = KENDO
                  KIINT  = KT21   + NMATAB(ISYMAB)*NUMII*NUMIJ
                  KJINT  = KIINT  + NCKI(ISYMCML)*NUMIJ
                  KV1MAT = KJINT  + NCKI(ISYMCML)*NUMIJ
                  KV2MAT = KV1MAT + MXKL*NUMII*NUMIJ
                  KAMAT  = KV2MAT + MXKL*NUMII*NUMIJ
                  KBMAT  = KAMAT  + MXKL*NUMII*NUMIJ
                  KENDJ  = KBMAT  + MXKL*NUMII*NUMIJ
                  LWRKJ  = LWORK  - KENDJ + 1
                  IF (LWRKJ .LT. 0) THEN
                 CALL QUIT('Insufficient memory in CC_CHOC1TERM in #J2')
                  ENDIF
C             
C------------------------------------------------------------------
C                 Extract integrals: 
C                    I(cm,l#j) = (cm|lj)
C                    J(cm,l#j) = (cl|jm)
C                 ordered for a j batch with a fixed ISYMJ
C------------------------------------------------------------------
C
*                 WRITE(LUPRI,*) 'Extracting I and J integrals'
                  CALL CCHO_RDINTC12(XOINT,WORK(KIINT),WORK(KJINT),IJ1,
     &                               NUMIJ,ISYMJ,IOFFJ)
*                   
*                 xointnorm = dnrm2(ntraoc(1),xoint,1) 
*                 ntot=ncki(isymcml)*numij
*                 xiintnrm = dsqrt( xiintnrm**2
*    &                     + dnrm2(ntot,work(kiint),1)**2 )
*                 xjintnrm = dsqrt ( xjintnrm**2
*    &                     + dnrm2(ntot,work(kjint),1)**2 )
*                 WRITE(LUPRI,*) 'After CCHO_RDINTC12:'
*                 WRITE(LUPRI,*) 'XOINT norm: ',xointnorm
*                 print '(8F20.15)',(xoint(ii),ii=1,ntraoc(1))
*                 WRITE(LUPRI,*) 'I norm: ',xiintnrm
*                 print '(8F20.15)',(work(ii),
*    &                              ii=kiint,kiint+ntraoc(1)-1 )
*                 WRITE(LUPRI,*) 'J norm: ',xjintnrm
*                 print '(8F20.15)',(work(ii),
*    &                              ii=kjint,kjint+ntraoc(1)-1 )
C            
C------------------------------------------------------------------
C                 Extract amplitudes:                        
C                    T21(ab,#i#j) = s(ai,bj)
C------------------------------------------------------------------
C            
*                 WRITE(LUPRI,*) 'Extracting T21'
                  CALL CCHO_RDAMPC1(T2VO,WORK(KT21),II1,NUMII,ISYMI,
     &                              IJ1,NUMIJ,ISYMJ)
*                 ntot=nmatab(isymab)*numii*numij
*                 xt21nrm=dsqrt( xt21nrm**2
*    &                   +       dnrm2(ntot,work(kt21),1)**2 )
*                 WRITE(LUPRI,*) 'T21 norm: ',xt21nrm
*     WRITE(LUPRI,*) 'T21:'
*     nnab=nmatab(isymab)
*     do jj=1,nrhf(isymj)
*        do ii=1,nrhf(isymi)
*     ioff=nnab*(numii*(jj-1)+ii-1)
*     print '(2I1,4F10.7)',ii,jj,(work(kk),
*    &                   kk=kt21+ioff,kt21+ioff+nnab-1)
*        enddo
*     enddo
*     ntot=nmatab(isymab)*numii*numij
*     print '(25I10)',(ILABEL(II),II=1,25)
*     print '(25G10.3)',(work(ii),ii=kt21,kt21+ntot-1)
C                 
                  DO ICHO = 1,NUMCHO
C                 
                     EC1COR = ZERO
*                    EC11COR = ZERO
*                    EC12COR = ZERO
C                     
C---------------------------------------------------
C                    Calculate Cholesky information.
C                       d(ab,#i) update vector 
C                       d(cml)   update vector
C---------------------------------------------------
C                     
*                 WRITE(LUPRI,*) 'Calculating Cholesky information'
*     WRITE(LUPRI,*) 'NUMCHO = ',numcho
*     print '(8F20.15)',(choele(ii),ii=1,numcho)
                     CALL CCHO_DECHOC1(FOCKD,CHOELE,NUMCHO,ICHO,
     &                                 WORK(KCHOO),WORK(KCHOV),
     &                                 ISYMCML,ISYMABI,
     &                                 II1,NUMII,ISYMI)
*                 ntoto = ncki(isymcml) 
*                 ntotv = nmatab(isymab)*numii
*                 xchoonrm = dnrm2(ntoto,work(kchoo),1) 
*                 xchovnrm = dnrm2(ntotv,work(kchov),1)
*                 WRITE(LUPRI,*) 'Cholesky vector #',icho
*                 WRITE(LUPRI,*) 'd(ab,#i) update norm:',xchovnrm
*     print '(5F20.15)',(work(ii),ii=kchov,kchov+ntotv-1)
*                 WRITE(LUPRI,*) 'd(cm,l) update norm:',xchoonrm
*     print '(5F20.15)',(work(ii),ii=kchoo,kchoo+ntoto-1)
*
*     WRITE(LUPRI,*) 'Resetting CHOV and CHOO to 1'
*     do ii=kchov,kchov+ntotv-1
*        work(ii)=one
*     enddo
*     do ii=kchoo,kchoo+ntoto-1
*        work(ii)=one
*     enddo 
* 
*                 WRITE(LUPRI,*) 'Recalculating using ccho_decho8 without sym'
*     WRITE(LUPRI,*) 'NUMCHO = ',numcho
*     print '(8F20.15)',(choele(ii),ii=1,numcho)
*                 call ccho_decho8(fockd,choele,numcho,icho,
*    &                             work(kchoo),work(kchov),1,1)
*                 ntot = nmatab(1)*numii
*                 xchoonrm = dnrm2(ntot,work(kchoo),1) 
*                 ntot = ncki(1) 
*                 xchovnrm = dnrm2(ntot,work(kchov),1)
*                 WRITE(LUPRI,*) 'Cholesky vector #',icho
*                 WRITE(LUPRI,*) 'd(ab,#i) update norm:',xchoonrm
*     print '(8F20.15)',(work(ii),ii=kchov,kchov+15)
*                 WRITE(LUPRI,*) 'd(cm,l) update norm:',xchovnrm
*     print '(8F20.15)',(work(ii),ii=kchoo,kchoo+15)
*     WRITE(LUPRI,*) 'T21 before scaling:'
*     ntot=nmatab(isymab)*numii*numij
*     print '(25I10)',(ILABEL(II),II=1,25)
*     print '(25G10.3)',(work(ii),ii=kt21,kt21+ntot-1)

C                     
C-----------------------------------------------------------------------
C                    Scale amplitudes and integrals with cholesky update
C-----------------------------------------------------------------------
C                     
*                 WRITE(LUPRI,*) 'Scale amplitudes'
                     CALL CCHO_SCLAMPC1(WORK(KT21),WORK(KCHOV),
     &                                  ISYMAB,NUMII,ISYMI,NUMIJ,ISYMJ)
*                 WRITE(LUPRI,*) 'Scale integrals'
                     CALL CCHO_SCLINTC1(WORK(KIINT),WORK(KJINT),
     &                                  WORK(KCHOO),ISYMCML,NUMIJ,IOFFJ)
*     ntot=nmatab(isymab)*numii*numij
*     WRITE(LUPRI,*) 'T21 after scaling:'
*     nnab=nmatab(isymab)
*     do jj=1,nrhf(isymj)
*        do ii=1,nrhf(isymi)
*     ioff=nnab*(numii*(jj-1)+ii-1)
*     print '(2I1,4F10.7)',ii,jj,(work(kk),
*    &                   kk=kt21+ioff,kt21+ioff+nnab-1)
*        enddo
*     enddo
*     print '(25I10)',(ILABEL(II),II=1,25)
*     print '(25G10.3)',(work(ii),ii=kt21,kt21+ntot-1)
*         WRITE(LUPRI,*) 'NTOT = ',ntot
*     if(icho .eq. 1) then
*         xt21sclnrm = dsqrt( xt21sclnrm**2
*    &               +        dnrm2(ntot,work(kt21),1)**2 )
*         WRITE(LUPRI,*) 'norm d(abi)*T21(ab,ij) = ',xt21sclnrm
*         ntot=ncki(isymcml)*numij
*         WRITE(LUPRI,*) 'NTOT = ',ntot
*         xiintsclnrm = dsqrt( xiintsclnrm**2
*    &                +        dnrm2(ntot,work(kiint),1)**2 )
*         xjintsclnrm = dsqrt( xjintsclnrm**2
*    &                +        dnrm2(ntot,work(kjint),1)**2 )
*         WRITE(LUPRI,*) 'norm d(cml)*I(cmlj) = ',xiintsclnrm
*         WRITE(LUPRI,*) 'norm d(cml)*J(cmlj) = ',xjintsclnrm
*     endif
C                    
                     DO ISYMK = 1,NSYM
*        WRITE(LUPRI,*)  'ISYMK = ',isymk
*        WRITE(LUPRI,*)  '----------------------------------'
C                    
                        IF( NRHF(ISYMK) .EQ. 0 ) GOTO 1236
C                       
                        ISYML = MULD2H(ISYMKL,ISYMK) 
                        ISYMCM = MULD2H(ISYML,ISYMJ)
C                        
                        IF( NT1AM(ISYMCM) .EQ. 0 ) GOTO 1236
C                        
C----------------------------------------------------------------
C                       Calculate the V, A, B and D intermediates
C----------------------------------------------------------------
C                       
*                 WRITE(LUPRI,*) 'Calculate the V, A and B'
                        NAB = NMATAB(ISYMAB)
                        NTOAB = MAX(NAB,1)
                        NIJ = NUMII*NUMIJ
                        NKL = NRHF(ISYMK)*NRHF(ISYML)
                        KOFFT2VO = IT2VO(ISYMAB,ISYMKL)
     &                        + NAB*IMATIJ(ISYMK,ISYML) + 1
C                       
                        CALL DGEMM('T','N',NIJ,NKL,NAB,ONE,
     &                             WORK(KT21),NTOAB,T2VO(KOFFT2VO),
     &                             NTOAB,ZERO,WORK(KV1MAT),NIJ)
C                       
                        NLK = NKL
                        KOFFT2VO = IT2VO(ISYMAB,ISYMKL)
     &                        + NAB*IMATIJ(ISYML,ISYMK) + 1
C                       
                        CALL DGEMM('T','N',NIJ,NLK,NAB,ONE,
     &                             WORK(KT21),NTOAB,T2VO(KOFFT2VO),
     &                             NTOAB,ZERO,WORK(KV2MAT),NIJ)
*     if(icho .eq. 1) then
*                       xvnrm=dsqrt( xvnrm**2
*    &                       +       dnrm2(nij*nkl,work(kvmat),1)**2 )
*                       WRITE(LUPRI,*) 'V norm: ',xvnrm
*     endif
C                       
                        NLJ = NRHF(ISYML)*NUMIJ
                        NTOLJ = MAX(NLJ,1)
                        NKI = NRHF(ISYMK)*NUMII
                        NCM = NT1AM(ISYMCM)
                        NTOCM = MAX(NCM,1)
                        KOFFIINT = KIINT + IOFFJ(ISYML)
                        KOFFLINT = KLINT + IOFFI(ISYMK)
                        KOFFJINT = KJINT + IOFFJ(ISYML)
                        KOFFCMAT = KCMAT + IOFFI(ISYMK)
                        KOFFMINT = KMINT + IOFFI(ISYMK)
C 
                        CALL DGEMM('T','N',NLJ,NKI,NCM,XMONE,
     &                             WORK(KOFFIINT),NTOCM,WORK(KOFFLINT),
     &                             NTOCM,ZERO,WORK(KAMAT),NTOLJ)
*     if(icho .eq. 1) then
*                       xanrm=dsqrt( xanrm**2
*    &                       +       dnrm2(nlj*nki,work(kamat),1)**2 )
*                       WRITE(LUPRI,*) 'A norm: ',xanrm
*     endif
C 
                        CALL DGEMM('T','N',NLJ,NKI,NCM,ONE,
     &                             WORK(KOFFJINT),NTOCM,WORK(KOFFCMAT),
     &                             NTOCM,ZERO,WORK(KBMAT),NTOLJ)
*     if(icho .eq. 1) then
*                       xbnrm=dsqrt( xbnrm**2
*    &                       +       dnrm2(nlj*nki,work(kbmat),1)**2 )
*                       WRITE(LUPRI,*) 'B norm: ',xbnrm
*     endif
C                         
C------------------------------------------------------
C                       Calculate energy contribution.
C------------------------------------------------------
C                       
*                 WRITE(LUPRI,*) 'Calculate energy contribution'
                        NIJ  = NUMII*NUMIJ  
                        NIJK = NIJ*NRHF(ISYMK)
                        NIJL = NIJ*NRHF(ISYML)
                        NLJ  = NRHF(ISYML)*NUMIJ
                        NLJK = NLJ*NRHF(ISYMK) 
C
C                       A <- A + B
C                       ----------
                        NTOT = NLJK*NUMII
*                    WRITE(LUPRI,*) 'NTOT = ',NTOT,' ONE = ',one
                        CALL DAXPY(NTOT,ONE,WORK(KBMAT),1,WORK(KAMAT),1)
*                    WRITE(LUPRI,*) 'A+B norm: ',dnrm2(nij*nkl,work(kamat),1)
C                       
C                       B <- D
C                       ------
                        CALL DGEMM('T','N',NLJ,NKI,NCM,ONE,
     &                             WORK(KOFFJINT),NTOCM,WORK(KOFFMINT),
     &                             NTOCM,ZERO,WORK(KBMAT),NTOLJ)
*     if(icho .eq. 1) then
*                       xdnrm=dsqrt( xdnrm**2
*    &                       +       dnrm2(nlj*nki,work(kbmat),1)**2 )
*                       WRITE(LUPRI,*) 'D norm: ',xdnrm
*     endif
C                       
                        DO I = 1,NUMII
                           DO K = 1,NRHF(ISYMK)
                              DO J = 1,NUMIJ
                                 DO L = 1,NRHF(ISYML) 
C                                   
                                    NIJKL = NIJK*(L - 1)
     &                                    + NIJ*(K - 1)
     &                                    + NUMII*(J - 1) + I
                                    NIJLK = NIJL*(K - 1)
     &                                    + NIJ*(L - 1)
     &                                    + NUMII*(J - 1) + I
                                    NLJKI = NLJK*(I - 1)
     &                                    + NLJ*(K - 1)
     &                                    + NRHF(ISYML)*(J - 1) + L
C                                   
                                    KOFFV1 = KV1MAT + NIJKL - 1
                                    KOFFV2 = KV2MAT + NIJLK - 1
                                    KOFFA  = KAMAT + NLJKI - 1
                                    KOFFB  = KBMAT + NLJKI - 1
C                                   
                                    EC1COR = EC1COR
     &                                     + WORK(KOFFV1)
     &                                      *WORK(KOFFA)
     &                                     + WORK(KOFFV2)
     &                                      *WORK(KOFFB)
C                                    
*                                   EC11COR = EC11COR
*    &                                     + WORK(KOFFV1)
*    &                                      *WORK(KOFFA)
*                                   EC12COR = EC12COR
*    &                                     + WORK(KOFFV2)
*    &                                      *WORK(KOFFB)
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
C                       
 1236                   CONTINUE
                     ENDDO
C                    enddo isymk
C                    
                     E4C1 = E4C1 + EC1COR
C 
              ENERGC(ICHO) = ENERGC(ICHO) + EC1COR
C 
C Decommented by Domenico
                    TNOW = SECOND()
                    DELTAT = TNOW - TLAST
                    TLAST = TNOW
                    SCNDSC(ICHO) = SCNDSC(ICHO) 
     &                          + DELTAT
C Decommented by Domenico
C 
*                    E4C11= E4C11+ EC11COR
*                    E4C12= E4C12+ EC12COR
*     if(icho .eq. 1) then
*        ECorIcho1 = ECorIcho1 + EC1COR
*        WRITE(LUPRI,*) 'Accumulated EC1COR(icho=1):',ECorIcho1
*     endif
*                    WRITE(LUPRI,*) 'E4C1 = ',E4C1,'  EC1COR = ',EC1COR
*                    WRITE(LUPRI,*) 'E4C11= ',E4C11,'  EC11COR = ',EC11COR
*                    WRITE(LUPRI,*) 'E4C12= ',E4C12,'  EC12COR = ',EC12COR
                     IF (DABS(EC1COR) .LT. THRCHO) GOTO 1
C                    
                  ENDDO
C                 enddo icho
  1               CONTINUE                  
               ENDDO
C              enddo ibatj
 1235          CONTINUE
            ENDDO
C           enddo isymj
         ENDDO
C        enddo ibati
 1234    CONTINUE
      ENDDO
C     enddo isymi
C
      RETURN
      END
C  /* Deck ccho_rdintc11 */
      SUBROUTINE CCHO_RDINTC11(XOINT,XCMAT,XLINT,XMINT,
     &                         II1,NUMII,ISYMI,IOFFI)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez April 2003.
C
C     Extract integrals: 
C        C(cm,k#i) =  (ik|mc) =  XOINT(mi,k;c) 
C        L(cm,k#i) = 2(ik|mc) - (mk|ic) = 2XOINT(mi,k;c) - XOINT(ik,m;c)
C        M(cm,k#i) =  (mk|ic) =  XOINT(ik,m;c)
C     ordered for a i batch with a fixed ISYMI      
C     from core XOINT(mi,k;c)=(cm|ik)
C        
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (ONE = 1.0D0,XMONE = -1.0D0,TWO = 2.0D0,THREE = 3.0D0)
      DIMENSION XOINT(*),XCMAT(*),XLINT(*),XMINT(*),IOFFI(*)
C     
      DO ISYMC = 1,NSYM
C     
         ISYMMK  = MULD2H(ISYMI,ISYMC)
         ISYMMIK = ISYMC
C        
         DO ISYMK = 1,NSYM
C        
            ISYMM  = MULD2H(ISYMMK,ISYMK)
            ISYMMI = MULD2H(ISYMK,ISYMC)
            ISYMIK = MULD2H(ISYMM,ISYMC)
            ISYMCM = MULD2H(ISYMC,ISYMM)
C           
               DO C = 1,NVIR(ISYMC)
                  DO K = 1,NRHF(ISYMK)
                     DO I = 1,NUMII
                        II = II1 + I - 1
                        DO M = 1,NRHF(ISYMM)
C                       
                           KMIKC = ISJIKA(ISYMMIK,ISYMC)
     &                           + NMAIJK(ISYMMIK)*(C - 1)
     &                           + IMAIJK(ISYMMI,ISYMK)
     &                           + NMATIJ(ISYMMI)*(K - 1)
     &                           + IMATIJ(ISYMM,ISYMI)
     &                           + NRHF(ISYMM)*(II - 1)
     &                           + M
                           KIKMC = ISJIKA(ISYMMIK,ISYMC)
     &                           + NMAIJK(ISYMMIK)*(C - 1)
     &                           + IMAIJK(ISYMIK,ISYMM)
     &                           + NMATIJ(ISYMIK)*(M - 1)
     &                           + IMATIJ(ISYMI,ISYMK)
     &                           + NRHF(ISYMI)*(K - 1)
     &                           + II
                           KCMKI = IOFFI(ISYMK)
     &                           + NT1AM(ISYMCM)
     &                             *(NRHF(ISYMK)*(I-1) + K - 1) 
     &                           + IT1AM(ISYMC,ISYMM)
     &                           + NVIR(ISYMC)*(M-1) + C
C   
C                          C <- (ik|mc)
C                          ------------
                           XCMAT(KCMKI) = XOINT(KMIKC)
C   
C                          L <- -(mk|ic)
C                          ------------
                           XLINT(KCMKI) = -XOINT(KIKMC)
C                          
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C              
         ENDDO
C        
      ENDDO
C     
      NTOT = NMAIJA(ISYMI)*NUMII
*     WRITE(LUPRI,*) 'C <- (ik|mc) norm: ',dnrm2(ntot,xcmat,1) 
*     WRITE(LUPRI,*) 'L <- -(mk|ic) norm: ',dnrm2(ntot,xlint,1)
C 
C     M <- -L = (mk|ic) 
C     -----------------
      CALL DCOPY(NTOT,XLINT,1,XMINT,1)
      CALL DSCAL(NTOT,XMONE,XMINT,1)
*     WRITE(LUPRI,*) 'M <- (mk|ic) norm: ',dnrm2(ntot,xmint,1)
C
C     L <- 2*C + L = 2*(ik|mc) - (mk|ic) 
C     ----------------------------------
      CALL DAXPY(NTOT,TWO,XCMAT,1,XLINT,1)
*     WRITE(LUPRI,*) 'L <- 2*(ik|mc) - (mk|ic) norm: ',dnrm2(ntot,xlint,1)
*     call dscal(ntot,10.0D10,xcmat,1)      
*     WRITE(LUPRI,*) 'C <- (ik|mc)*10^10 norm: ',dnrm2(ntot,xcmat,1) 
C     
      RETURN
      END
C  /* Deck ccho_rdintc12 */
      SUBROUTINE CCHO_RDINTC12(XOINT,XIINT,XJINT,IJ1,NUMIJ,ISYMJ,IOFFJ)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez April 2003.
C
C     Extract integrals: 
C        I(cm,l#j) = (cm|lj) = XOINT(ml,j;c)
C        J(cm,l#j) = (cl|jm) = XOINT(lj,m;c)
C     ordered for a j batch with a fixed ISYMJ      
C     from core XOINT(ml,j;c)=(cm|lj)
C        
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (ONE = 1.0D0, XMONE = -1.0D0, THREE = 3.0D0)
      DIMENSION XOINT(*),XIINT(*),XJINT(*),IOFFJ(*)
C     
      DO ISYMC = 1,NSYM
C     
         ISYMML = MULD2H(ISYMC,ISYMJ)
         ISYMMLJ = ISYMC
C        
         DO ISYML = 1,NSYM
C        
            ISYMLJ = MULD2H(ISYML,ISYMJ)
            ISYMCM = MULD2H(ISYML,ISYMJ)
            ISYMM  = MULD2H(ISYMCM,ISYMC)
C           
               DO C = 1,NVIR(ISYMC)
                  DO J = 1,NUMIJ
                     IJ = IJ1 + J - 1
                     DO L = 1,NRHF(ISYML)
                        DO M = 1,NRHF(ISYMM)
C                       
                           KMLJC = ISJIKA(ISYMMLJ,ISYMC)
     &                           + NMAIJK(ISYMMLJ)*(C - 1)
     &                           + IMAIJK(ISYMML,ISYMJ)
     &                           + NMATIJ(ISYMML)*(IJ - 1)
     &                           + IMATIJ(ISYMM,ISYML)
     &                           + NRHF(ISYMM)*(L - 1)
     &                           + M
                           KLJMC = ISJIKA(ISYMMLJ,ISYMC)
     &                           + NMAIJK(ISYMMLJ)*(C - 1)
     &                           + IMAIJK(ISYMLJ,ISYMM)
     &                           + NMATIJ(ISYMLJ)*(M - 1)
     &                           + IMATIJ(ISYML,ISYMJ)
     &                           + NRHF(ISYML)*(IJ - 1)
     &                           + L
                           KCMLJ = IOFFJ(ISYML)
     &                           + NT1AM(ISYMCM)
     &                             *(NRHF(ISYML)*(J-1) + L - 1)
     &                           + IT1AM(ISYMC,ISYMM)
     &                           + NVIR(ISYMC)*(M-1) + C
C                          
                           XIINT(KCMLJ) = XOINT(KMLJC)
                           XJINT(KCMLJ) = XOINT(KLJMC)
C                           
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C                    
         ENDDO         
      ENDDO
C     
      RETURN
      END
C  /* Deck ccho_rdampc1 */
      SUBROUTINE CCHO_RDAMPC1(T2VO,T21,II1,NUMII,ISYMI,IJ1,NUMIJ,ISYMJ)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez April 2003.
C
C     Extract amplitudes:
C        T21(ab,#i#j) = s(ai,bj) = 2*T2VO(ab,ij) - T2VO(ab,ji)
C        
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION T2VO(*),T21(*)
C      
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (ONE = 1.0D0, XMONE = -1.0D0, TWO = 2.0D0)
C     
      ISYMIJ = MULD2H(ISYMI,ISYMJ)
      ISYMAB = MULD2H(1,ISYMIJ)
C     
      NAB = NMATAB(ISYMAB)
C     
      DO I = 1,NUMII
         II = II1 + I - 1
         DO J = 1,NUMIJ
            IJ = IJ1 + J - 1
C           
            KIJ = IMATIJ(ISYMI,ISYMJ)
     &          + NRHF(ISYMI)*(IJ - 1)
     &          + II
            KABIJ = IT2VO(ISYMAB,ISYMIJ) + (KIJ - 1)*NAB + 1
C           
            KJI = IMATIJ(ISYMJ,ISYMI)
     &          + NRHF(ISYMJ)*(II - 1)
     &          + IJ
            KABJI = IT2VO(ISYMAB,ISYMIJ) + (KJI - 1)*NAB + 1
C           
            KT21IJ = NUMII*(J - 1) + I
            IOFFT21 = (KT21IJ - 1)*NAB + 1
C           
            CALL DCOPY(NAB,T2VO(KABJI),1,T21(IOFFT21),1)
*      WRITE(LUPRI,*) 'RDAMP: COPY -- I=',I,' J=',J            
*      print '(2I1,4F10.7)',i,j,(t21(kk),kk=iofft21,iofft21+nab-1)
            CALL DSCAL(NAB,XMONE,T21(IOFFT21),1)
*      WRITE(LUPRI,*) 'RDAMP: DSCAL -- I=',I,' J=',J            
*      print '(2I1,4F10.7)',i,j,(t21(kk),kk=iofft21,iofft21+nab-1)
            CALL DAXPY(NAB,TWO,T2VO(KABIJ),1,T21(IOFFT21),1)
*      WRITE(LUPRI,*) 'RDAMP: DAXPY -- I=',I,' J=',J            
*      print '(2I1,4F10.7)',i,j,(t21(kk),kk=iofft21,iofft21+nab-1)
C           
         ENDDO
      ENDDO
      RETURN
      END
C  /* Deck ccho_dechoc1 */
      SUBROUTINE CCHO_DECHOC1(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                       ISYDLI,ISYBCJ,IJ1,NUMIJ,ISYMJ)
C
C     JLC, BFR, TBP, HK, AS, April 2003.
C
C     Construct vector for Cholesky decomposition:
C
C        For the occupied part, OCCHO(dli)  update vector
C        For the virtual  part, VICHO(bc#j)  update vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C-------------------
C     Occupied part.
C-------------------
C
      IF (ICHO .EQ. 1) THEN
         IND = 0
         DO ISYMI = 1,NSYM
            ISYMDL = MULD2H(ISYMI,ISYDLI)
            DO I = 1,NRHF(ISYMI)
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  DO L = 1,NRHF(ISYML)
                     DO D = 1,NVIR(ISYMD)
C
                        KOFFI = IRHF(ISYMI) + I
                        KOFFL = IRHF(ISYML) + L
                        KOFFD = IVIR(ISYMD) + D
C
                        OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) = 
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
         NDIMOC = IND
      ELSE
         IND = 0
         DO ISYMI = 1,NSYM
            ISYMDL = MULD2H(ISYMI,ISYDLI)
            DO I = 1,NRHF(ISYMI)
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  DO L = 1,NRHF(ISYML)
                     DO D = 1,NVIR(ISYMD)
C
                        KOFFI = IRHF(ISYMI) + I
                        KOFFL = IRHF(ISYML) + L
                        KOFFD = IVIR(ISYMD) + D
C
                        OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
         NDIMOC = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMOC,FACTOR,OCCHO,1)
C
      ENDIF
C
C------------------
C     Virtual part.
C------------------
C
      IF (ICHO .EQ. 1) THEN
         IND = 0
C        DO ISYMJ = 1,NSYM
            ISYMBC = MULD2H(ISYMJ,ISYBCJ)
            DO J = 1,NUMIJ
               DO ISYMC = 1,NSYM
                  ISYMB = MULD2H(ISYMC,ISYMBC)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
C
                        KOFFJ = IRHF(ISYMJ) + IJ1 + J - 1
                        KOFFC = IVIR(ISYMC) + C
                        KOFFB = IVIR(ISYMB) + B
C
                        OME = FOCKD(KOFFB)+FOCKD(KOFFC)-FOCKD(KOFFJ)
C
                        IND = IND + 1
                        VICHO(IND) = 
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
C        ENDDO
         NDIMVI = IND
C
      ELSE
C
         IND = 0
C        DO ISYMJ = 1,NSYM
            ISYMBC = MULD2H(ISYMJ,ISYBCJ)
            DO J = 1,NUMIJ
               DO ISYMC = 1,NSYM
                  ISYMB = MULD2H(ISYMC,ISYMBC)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
C
                        KOFFJ = IRHF(ISYMJ) + IJ1 + J - 1
                        KOFFC = IVIR(ISYMC) + C
                        KOFFB = IVIR(ISYMB) + B
C
                        OME = FOCKD(KOFFB)+FOCKD(KOFFC)-FOCKD(KOFFJ)
C
                        IND = IND + 1
                        VICHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
C        ENDDO
C
         NDIMVI = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMVI,FACTOR,VICHO,1)
C         
      ENDIF
C
      RETURN
      END
C  /* Deck ccho_sclampc1 */
      SUBROUTINE CCHO_SCLAMPC1(T21,CHOV,ISYMAB,NUMII,ISYMI,NUMIJ,ISYMJ)
C
C     JLC, BFR, TBP, HK, AS, April 2003.
C
C     Scale s(ab,#i#j) with cholesky update vector: 
C        s(ab,#i#j)d(ab#i)
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION T21(*),CHOV(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      NABI = NMATAB(ISYMAB)*NUMII
C     
      DO J = 1,NUMIJ
         DO IABI = 1,NABI
C        
            KABIJ = NABI*(J-1) + IABI
C           
            T21(KABIJ) = T21(KABIJ)*CHOV(IABI)
C           
         ENDDO
      ENDDO
C     
      RETURN
      END
C  /* Deck ccho_sclintc1 */
      SUBROUTINE CCHO_SCLINTC1(XIINT,XJINT,CHOO,ISYMCML,NUMIJ,IOFFJ)
C
C     JLC, BFR, TBP, HK, AS, April 2003.
C
C     Scale (cm|jl) and (cl|jm) with cholesky update vector: 
C        d(cml)*I(cm,l#j)
C        d(cml)*J(cm,l#j)
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XIINT(*),XJINT(*),CHOO(*),IOFFJ(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      DO ISYML = 1,NSYM
         ISYMCM = MULD2H(ISYMCML,ISYML)
         NCM = NT1AM(ISYMCM)
*        WRITE(LUPRI,*) 'ISYMCM = ',ISYMCM,' NCM = ',NCM
*        WRITE(LUPRI,*) 'IOFFJ(',ISYML,') = ',ioffj(isyml)
         DO J = 1,NUMIJ
*           WRITE(LUPRI,*) 'J = ',j,' of ',numij         
            DO L = 1,NRHF(ISYML)
*              WRITE(LUPRI,*) 'L = ',l
               DO ICM = 1,NCM
*                 WRITE(LUPRI,*) 'ICM = ',icm
*                 WRITE(LUPRI,*) 'IOFFJ(',ISYML,') = ',ioffj(isyml)
*                 call flush_(6)
C           
                  KCMLJ = IOFFJ(ISYML) + NCM*(NRHF(ISYML)*(J-1) + L - 1)
     &                  + ICM
*                 WRITE(LUPRI,*) 'KCMLJ = ',KCMLJ
*                 call flush_(6)
                  KCML  = ICKI(ISYMCM,ISYML) + NCM*(L-1) + ICM
*                 WRITE(LUPRI,*) 'KCML = ',KCML
*                 call flush_(6)
C              
                  XIINT(KCMLJ) = XIINT(KCMLJ)*CHOO(KCML)
                  XJINT(KCMLJ) = XJINT(KCMLJ)*CHOO(KCML) 
C              
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      RETURN
      END
C  /* Deck cc_choc2term */
      SUBROUTINE CC_CHOC2TERM(XOINT,T2VO,FOCKD,NUMCHO,CHOELE,WORK,LWORK,
     &                       E4C2,FBATCH,IPRINT,NONI,NONJ)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez April 2003.
C
C        E4C2 = Sum(#i#j) A(#i#j) * B(#i#j)
C     
C     where
C        
C        A(#i#j) =   Sum(abk) T21(abk;#i)* d(abk)*T22(abk;#j)
C        B(#i#j) = - Sum(clm) I(clm;#i)* d(clm)*L(clm;#j) 
C     
C        T21(abk;i) = s(ak,bi) = 2*t(ak,bi) - t(ai,bk)
C                   = 2*T2VO(ab,ki) - T2VO(ab,ik)      
C        T22(abk;j) = t(ak,bj) = T2VO(ab,kj)
C        I(clm;i)   = (cl|im) = XOINT(li,m;c)
C        L(clm;j)   = L_{mjlc} = 2*(mj|lc) - (mc|lj) 
C                   = 2*XOINT(lm,j;c) - XOINT(ml,j;c)
C     
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XOINT(*),T2VO(*),FOCKD(*),CHOELE(*)
      DIMENSION WORK(LWORK)
      LOGICAL   FBATCH,IPRINT
C
C     DIMENSION IOFFT21(8),IOFFI(8),IOFFT22(8),IOFFL(8) 
C     
#include <ccorb.h>
#include <ccsdsym.h>
#include <cc_cho.h>
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
C 
      TLAST = SECOND()
C 
      IF (IPRINT) THEN
         TIMT = SECOND()
         WRITE(LUPRI,'(A,/,A,/)')
     &   'Calculation of the C2 term:',
     &   '==========================='
      ENDIF
*     WRITE(LUPRI,*) 'T2VO:'
*     ntot=nrhft**2*nvirt**2
*     print '(8F10.7)',(T2VO(ii),ii=1,ntot)
*     nnab=nmatab(1)
*     do jj=1,nrhft
*        do ii=1,nrhft
*     ioff=nnab*(nrhft*(jj-1)+ii-1)
*     print '(2I1,4F10.7)',ii,jj,(t2vo(kk),
*    &                   kk=ioff+1,ioff+nnab)
*        enddo
*     enddo
* 
*     xt2norm = dnrm2(nt2sq(1),t2vo,1) 
*     xointnorm = dnrm2(ntraoc(1),xoint,1) 
*     WRITE(LUPRI,*) 'Initial T2VO norm: ',xt2norm
*     WRITE(LUPRI,*) 'Initial XOINT norm: ',xointnorm
C     
      DO ISYMI = 1,NSYM
*        WRITE(LUPRI,*)  'ISYMI = ',isymi
*        WRITE(LUPRI,*)  '----------------------------------'
C     
         IF ( NRHF(ISYMI) .EQ. 0 ) GOTO 1234
C      
         ISYMABK = ISYMI
         ISYMCLM = ISYMI
         ISYMJ   = ISYMI
C            
C        Dynamic allocation for d(abk) and d(clm)
C        ----------------------------------------
         KCHOV = 1
         KCHOO = KCHOV + NCKASR(ISYMABK)
         KEND1 = KCHOO + NCKI(ISYMCLM)
         LWRK1 = LWORK - KEND1 + 1
         IF (LWRK1 .LT. 0) THEN
             CALL QUIT('Insufficient memory in CC_CHOC2TERM d')
         ENDIF
C                
         LWORKI = LWRK1/2
C         
C-------------------------------
C        Batch over the I index.
C-------------------------------
C         
         LENI = NCKASR(ISYMABK) + NCKI(ISYMCLM)
         IF (FBATCH) THEN
            NEFI   = MIN(NONI,NRHF(ISYMI))
            LEFF   = NEFI*LENI + 1
            LWORKI = MIN(LWORKI,LEFF)
         ENDIF
         NUMI = MIN(NRHF(ISYMI),LWORKI/LENI)
         IF (NUMI .EQ. 0) THEN
            WRITE(LUPRI,*) 'NUMI .EQ. 0 IN CC_CHOPTC1'
            CALL QUIT('NOT ENOUGH SPACE IN CC_CHOPTC1 !!')
         ENDIF
C         
         NBATI = (NRHF(ISYMI)-1)/NUMI + 1
C
C
         IF (IPRINT) THEN
            WRITE(LUPRI,'(3X,A,I1,A,/,3X,A)')
     &      'Batch over I, symmetry ',ISYMI,':',
     &      '-------------------------'
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10)')
     &      'Minimum work space required   : ',LENI,
     &      'Work space available for batch: ',LWORKI
            WRITE(LUPRI,'(3X,A,I10,/,3X,A,I10,/)')
     &      'Number of occupied orbitals   : ',NRHF(ISYMI),
     &      'Required number of J-batches  : ',NBATI
         ENDIF
C        
         II2=0
         DO IBATI = 1,NBATI
*           WRITE(LUPRI,*)  'IBATI = ',ibati
*           xcount = xcount + one
            II1 = II2 + 1
            II2 = II2 + NUMI
            IF (II2 .GT. NRHF(ISYMI)) II2 = NRHF(ISYMI)
            NUMII = II2 - II1 + 1
C
            IF (IPRINT) THEN
               WRITE(LUPRI,'(6X,A,I10,A,/,6X,A)')
     &         'I-batch number ',IBATI,':',
     &         '--------------------------'
               WRITE(LUPRI,'(6X,A,I10,1X,I10,/)')
     &         'First and last I: ',II1,II2
            ENDIF
C           
C           IOFFT21(8) offset
C           ---------------
C           ICOUNT = 0
C           DO ISYMK = 1,NSYM
C              ISYMAB = MULD2H(ISYMI,ISYMK)
C              IOFFT21(ISYMK) = ICOUNT
C              ICOUNT = ICOUNT + NMATAB(ISYMAB)*NRHF(ISYMK)*NUMII
C           ENDDO
*           WRITE(LUPRI,*) 'IOFFT21 : ',iofft21
C              
C           IOFFI(8) offset
C           ---------------
C           ICOUNT = 0
C           DO ISYMM = 1,NSYM
C              ISYMCL = MULD2H(ISYMI,ISYMM)
C              IOFFI(ISYMM) = ICOUNT
C              ICOUNT = ICOUNT + NT1AM(ISYMCL)*NRHF(ISYMM)*NUMII
C           ENDDO
*           WRITE(LUPRI,*) 'IOFFI : ',ioffi
C            
C------------------------------
C           Dynamic allocation.
C------------------------------
C
            KT21  = KEND1
            KIINT = KT21 + NCKASR(ISYMABK)*NUMII
            KENDI = KIINT + NCKI(ISYMCLM)*NUMII
            LWRKI = LWORK - KENDI + 1
            IF (LWRKI .LT. 0) THEN
         CALL QUIT('Insufficient memory in CC_CHOC2TERM in the I batch')
            ENDIF
C            
C-----------------------------------------
C           Extract T21(abk;#i) = s(ak,bi)
C-----------------------------------------
C            
            CALL CCHO_RDT21C2(T2VO,WORK(KT21),II1,NUMII,ISYMI)
*     WRITE(LUPRI,*) 'T21:'
*     nnab=nmatab(1)
*     do jj=1,nrhf(isymj)
*        do ii=1,nrhf(isymi)
*     ioff=nnab*(numii*(jj-1)+ii-1)
*     print '(2I1,4F10.7)',ii,jj,(work(kk),
*    &                   kk=kt21+ioff,kt21+ioff+nnab-1)
*        enddo
*     enddo
*              ntot=nckasr(isymabk)*numii
*              xt21nrm=dnrm2(ntot,work(kt22),1) 
*              WRITE(LUPRI,*) 'T21(abk;#i) norm:',xt22nrm
C            
C--------------------------------------
C           Extract I(clm;#i) = (cl|im)
C--------------------------------------
C            
            CALL CCHO_RDIC2(XOINT,WORK(KIINT),II1,NUMII,ISYMI)
*                 xiintnorm = dnrm2(ntraoc(1),work(kiint),1)
*                 WRITE(LUPRI,*) 'XIINT norm: ',xiintnorm
C               
C---------------------------
C           Batch over J.
C---------------------------
C              
            LENJ = NCKASR(ISYMABK) + NCKI(ISYMCLM) + 2*NUMII
C              
            LWORKJ = LWORK - KENDI + 1
            IF (FBATCH) THEN
               NEFJ = MIN(NRHF(ISYMJ),NONJ)
               LEFF = NEFJ*LENJ + 1
               LWORKJ = MIN(LWORKJ,LEFF)
            ENDIF
            NUMJ = MIN(NRHF(ISYMJ),LWORKJ/LENJ)
            IF (NUMJ .EQ. 0) THEN
               WRITE(LUPRI,*) 'NUMJ .EQ. 0 in CC_CHOPTC2'
               CALL QUIT('Not enough space in CC_CHOPTC2')
            ENDIF
C                 
            NBATJ = (NRHF(ISYMJ)-1)/NUMJ + 1
C
            IF (IPRINT) THEN
               WRITE(LUPRI,'(9X,A,I1,A,/,9X,A)')
     &         'Batch over J, symmetry ',ISYMJ,':',
     &         '-------------------------'
               WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10)')
     &         'Minimum work space required   : ',LENJ,
     &         'Work space available for batch: ',LWORKJ
               WRITE(LUPRI,'(9X,A,I10,/,9X,A,I10,/)')
     &         'Number of occupied orbitals    : ',NRHF(ISYMJ),
     &         'Required number of J-batches  : ',NBATJ
            ENDIF
C              
            IJ2 = 0
            DO IBATJ = 1,NBATJ
               IJ1 = IJ2 + 1
               IJ2 = IJ2 + NUMJ
               IF(IJ2 .GT. NRHF(ISYMJ)) IJ2 = NRHF(ISYMJ)
               NUMIJ = IJ2 - IJ1 + 1
C
               IF (IPRINT) THEN
                  WRITE(LUPRI,'(12X,A,I10,A,/,12X,A)')
     &            'J-batch number ',IBATJ,':',
     &            '--------------------------'
                  WRITE(LUPRI,'(12X,A,I10,1X,I10,/)')
     &            'First and last J: ',IJ1,IJ2
               ENDIF
C           
C              IOFFT22(8) offset
C              ---------------
C              ICOUNT = 0
C              DO ISYMK = 1,NSYM
C                 ISYMAB = MULD2H(ISYMJ,ISYMK)
C                 IOFFT22(ISYMK) = ICOUNT
C                 ICOUNT = ICOUNT + NMATAB(ISYMAB)*NRHF(ISYMK)*NUMIJ
C              ENDDO
*              WRITE(LUPRI,*) 'IOFFT22 : ',iofft22
C                  
C              IOFFL(8) offset
C              ---------------
C              ICOUNT = 0
C              DO ISYMM = 1,NSYM
C                 ISYMCL = MULD2H(ISYMM,ISYMJ)
C                 IOFFL(ISYMM) = ICOUNT
C                 ICOUNT = ICOUNT + NT1AM(ISYMCL)*NRHF(ISYMM)*NUMIJ
C              ENDDO
*              WRITE(LUPRI,*) 'IOFFL : ',ioffl
C                  
C------------------------------------
C              Dynamic allocation.
C------------------------------------
C                  
               KT22  = KENDI
               KLINT = KT22 + NCKASR(ISYMABK)*NUMIJ
               KAMAT = KLINT + NCKI(ISYMCLM)*NUMIJ
               KBMAT = KAMAT + NUMII*NUMIJ
               KENDJ = KBMAT + NUMII*NUMIJ
               LWRKJ = LWORK - KENDJ + 1
               IF (LWRKJ .LT. 0) THEN
                  CALL QUIT('Insufficient memory in CC_CHOC2TERM in #J')
               ENDIF
C               
C--------------------------------------------
C              Extract T22(abk;#j) = t(ak,bj)
C--------------------------------------------
C               
               CALL CCHO_RDT22C2(T2VO,WORK(KT22),IJ1,NUMIJ,ISYMJ)
*              ntot=nckasr(isymabk)*numij
*              xt22nrm=dnrm2(ntot,work(kt22),1) 
*              WRITE(LUPRI,*) 'T22(abk;#j) norm:',xt22nrm
C               
C--------------------------------------------------
C              Extract L(clm;j) = 2*(mj|lc)-(mc|lj)
C--------------------------------------------------
C               
               CALL CCHO_RDLC2(XOINT,WORK(KLINT),IJ1,NUMIJ,ISYMJ)
*              ntot=ncki(isymclm)*numij
*              xlintnorm = dnrm2(ntot,work(klint),1)
*              WRITE(LUPRI,*) 'L(clm;#j) norm: ',xlintnorm
C              
               DO ICHO = 1,NUMCHO
C              
                  EC2COR = ZERO
C                  
C-----------------------------------------------
C                 Calculate Cholesky information
C                    d(abk) update vector
C                    d(clm) update vector
C-----------------------------------------------
C                  
                  CALL CCHO_DECHOC2(FOCKD,CHOELE,NUMCHO,ICHO,
     &                              WORK(KCHOO),WORK(KCHOV),
     &                              ISYMCLM,ISYMABK)
*                 ntoto = ncki(isymclm) 
*                 WRITE(LUPRI,*) 'NCKI(isymclm) = ',ncki(isymclm)
*                 ntotv = nckasr(isymabk)
*                 WRITE(LUPRI,*) 'NCKASR(isymabk) = ',nckasr(isymabk)
*                 xchoonrm = dnrm2(ntoto,work(kchoo),1) 
*                 xchovnrm = dnrm2(ntotv,work(kchov),1)
*                 WRITE(LUPRI,*) 'Cholesky vector #',icho
*                 WRITE(LUPRI,*) 'd(abk) update norm:',xchovnrm
*                 WRITE(LUPRI,*) 'd(clm) update norm:',xchoonrm
C                  
C-----------------------------------------------------
C                 Scale T22 and L with cholesky update
C-----------------------------------------------------
C                  
                  CALL CCHO_SCLT22C2(WORK(KT22),WORK(KCHOV),
     &                               ISYMJ,NUMIJ)
                  CALL CCHO_SCLLC2(WORK(KLINT),WORK(KCHOO),
     &                             ISYMJ,NUMIJ) 
*                 ntot=nckasr(isymabk)*numij
*                 xt22sclnrm=dnrm2(ntot,work(kt22),1) 
*                 WRITE(LUPRI,*) 'd(abk)*T22(abk;#j) norm:',xt22sclnrm
*                 ntot=ncki(isymclm)*numij
*                 xlsclnrm=dnrm2(ntot,work(klint),1) 
*                 WRITE(LUPRI,*) 'd(clm)*L(clm;#j) norm:',xlsclnrm
C                  
C-------------------------------------------------
C                 Calculate A and B intermeditates
C-------------------------------------------------
C                  
                  NABK   = NCKASR(ISYMABK)
                  NTOABK = MAX(1,NABK)
C                 
                  CALL DGEMM('T','N',NUMII,NUMIJ,NABK,ONE,
     &                       WORK(KT21),NTOABK,WORK(KT22),
     &                       NTOABK,ZERO,WORK(KAMAT),NUMII)
C                 
                  NCLM   = NCKI(ISYMCLM)
                  NTOCLM = MAX(1,NCLM)
                  CALL DGEMM('T','N',NUMII,NUMIJ,NCLM,XMONE,
     &                       WORK(KIINT),NTOCLM,WORK(KLINT),
     &                       NTOCLM,ZERO,WORK(KBMAT),NUMII)
C                  
C----------------------------------------------
C                 Calculate energy contribution
C----------------------------------------------
C                  
                  NIJ = NUMII*NUMIJ
                  EC2COR = EC2COR
     &                   + DDOT(NIJ,WORK(KAMAT),1,WORK(KBMAT),1)
C                 
                  E4C2 = E4C2 + EC2COR
C 
              ENERGC(ICHO) = ENERGC(ICHO) + EC2COR
C 
C Decommented by Domenico
                 TNOW = SECOND()
                 DELTAT = TNOW - TLAST
                 TLAST = TNOW
                 SCNDSC(ICHO) = SCNDSC(ICHO) 
     &                         + DELTAT
C Decommented by Domenico 
*                 xanrm=dnrm2(nij,work(kamat),1)
*                 WRITE(LUPRI,*) 'A norm:',xanrm
*                 xbnrm=dnrm2(nij,work(kbmat),1)
*                 WRITE(LUPRI,*) 'B norm:',xbnrm
*                 WRITE(LUPRI,*) 'E4C2 = ',E4C2,'  EC2COR = ',EC2COR
                  IF (DABS(EC2COR) .LT. THRCHO) GOTO 1
C                 
               ENDDO
C              enddo icho
  1            CONTINUE
            ENDDO
C           enddo ibatj
         ENDDO
C        enddo ibati
 1234    CONTINUE
      ENDDO
C     enddo isymi
C     
      RETURN
      END

C  /* Deck ccho_rdt21c2 */
      SUBROUTINE CCHO_RDT21C2(T2VO,T21,II1,NUMII,ISYMI)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez April 2003.
C
C     Extract amplitudes:
C        T21(abk;#i) = s(ak,bi) = 2*T2VO(ab,ki) - T2VO(ab,ik)
C        
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (ONE = 1.0D0, XMONE = -1.0D0, TWO = 2.0D0)
      DIMENSION T2VO(*),T21(*)
C     
      ISYMT = 1
      ISYMABK = MULD2H(ISYMT,ISYMI)
C     
      DO ISYMK = 1,NSYM
C     
         ISYMKI = MULD2H(ISYMK,ISYMI)
         ISYMAB = MULD2H(ISYMT,ISYMKI)
C        
         NAB = NMATAB(ISYMAB) 
C        
         DO I = 1,NUMII
            II = II1 + I - 1
            DO K = 1,NRHF(ISYMK)
C           
               NKI = IMATIJ(ISYMK,ISYMI) + NRHF(ISYMK)*(II - 1) + K
               NABKI = IT2VO(ISYMAB,ISYMKI) + NAB*(NKI - 1) + 1
C              
               NIK = IMATIJ(ISYMI,ISYMK) + NRHF(ISYMI)*(K - 1) + II
               NABIK = IT2VO(ISYMAB,ISYMKI) + NAB*(NIK - 1) + 1
C              
               NT21 = NCKASR(ISYMABK)*(I - 1)
     &              + ICKASR(ISYMAB,ISYMK)
     &              + NAB*(K - 1) + 1
C              
               CALL DCOPY(NAB,T2VO(NABIK),1,T21(NT21),1)
               CALL DSCAL(NAB,XMONE,T21(NT21),1)
               CALL DAXPY(NAB,TWO,T2VO(NABKI),1,T21(NT21),1) 
C              
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C  /* Deck ccho_rdic2 */
      SUBROUTINE CCHO_RDIC2(XOINT,XIINT,II1,NUMII,ISYMI)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez April 2003.
C
C           Extract I(clm;#i) = (cl|im) = XOINT(li,m;c)
C     
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      DIMENSION XOINT(*),XIINT(*)
C     
      ISYMCLM = ISYMI
C     
      DO ISYMC = 1,NSYM
C     
         ISYMLIM = ISYMC
C        
         DO ISYMM = 1,NSYM
C        
            ISYMCL = MULD2H(ISYMCLM,ISYMM)
            ISYML  = MULD2H(ISYMCL,ISYMC)
            ISYMLI = MULD2H(ISYMLIM,ISYMM)
            ISYMI  = MULD2H(ISYMLI,ISYML)
C           
            DO C = 1,NVIR(ISYMC)
               DO M = 1,NRHF(ISYMM)
                  DO I = 1,NUMII
                     II = II1 + I - 1
                        DO L = 1,NRHF(ISYML)
C                       
                           NLIMC = ISJIKA(ISYMLIM,ISYMC)
     &                           + NMAIJK(ISYMLIM)*(C - 1)
     &                           + IMAIJK(ISYMLI,ISYMM)
     &                           + NMATIJ(ISYMLI)*(M - 1)
     &                           + IMATIJ(ISYML,ISYMI)
     &                           + NRHF(ISYML)*(II - 1)
     &                           + L
                           NCLMI = NCKI(ISYMCLM)*(I - 1)
     &                           + ICKI(ISYMCL,ISYMM)
     &                           + NT1AM(ISYMCL)*(M - 1)
     &                           + IT1AM(ISYMC,ISYML)
     &                           + NVIR(ISYMC)*(L - 1)
     &                           + C
C                          
                           XIINT(NCLMI) = XOINT(NLIMC)
C                          
                        ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      RETURN
      END
C  /* Deck ccho_rdt22c2 */
      SUBROUTINE CCHO_RDT22C2(T2VO,T22,IJ1,NUMIJ,ISYMJ)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez April 2003.
C
C     Extract T22(abk;#j) = t(ak,bj) = T2VO(ab,kj)
C        
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (ONE = 1.0D0, XMONE = -1.0D0, TWO = 2.0D0)
      DIMENSION T2VO(*),T22(*)
C     
      ISYMT = 1
      ISYMABK = MULD2H(ISYMT,ISYMJ)
C     
      DO ISYMK = 1,NSYM
C     
         ISYMAB = MULD2H(ISYMABK,ISYMK)
         ISYMKJ = MULD2H(ISYMK,ISYMJ) 
C        
         NAB = NMATAB(ISYMAB) 
C        
         DO J = 1,NUMIJ
            IJ = IJ1 + J - 1
            DO K = 1,NRHF(ISYMK) 
C           
               NKJ = IMATIJ(ISYMK,ISYMJ) + NRHF(ISYMK)*(IJ - 1) + K
               NABKJ = IT2VO(ISYMAB,ISYMKJ) + NAB*(NKJ - 1) + 1
C              
               NT22 = NCKASR(ISYMABK)*(J - 1)
     &              + ICKASR(ISYMAB,ISYMK)
     &              + NAB*(K - 1) + 1
C              
               CALL DCOPY(NAB,T2VO(NABKJ),1,T22(NT22),1) 
C              
            ENDDO
         ENDDO
      ENDDO
C     
      RETURN
      END

C  /* Deck ccho_rdlc2 */
      SUBROUTINE CCHO_RDLC2(XOINT,XLINT,IJ1,NUMIJ,ISYMJ)
C
C     Javier Lopez Cacheiro, Berta Fernandez, Thomas Bondo Pedersen,
C     Henrik Koch, Alfredo Sanchez April 2003.
C
C           Extract L(clm;#j) = L_{mjlc} = 2*(mj|lc) - (mc|lj) 
C                             = 2*XOINT(lm,j;c) - XOINT(ml,j;c)
C     
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XOINT(*),XLINT(*)
C     
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C      
      ISYMCLM = ISYMJ
C     
      DO ISYMC = 1,NSYM
C     
         ISYMLMJ = ISYMC
         ISYMLM  = MULD2H(ISYMLMJ,ISYMJ)
         ISYMMLJ = ISYMLMJ
         ISYMML  = ISYMLM
C        
         DO ISYML = 1,NSYM
C        
            ISYMM = MULD2H(ISYMML,ISYML) 
            ISYMCL = MULD2H(ISYMC,ISYML)
*           WRITE(LUPRI,*) 'ISYMM =',isymm
*           WRITE(LUPRI,*) 'NRHF(ISYMM)=',nrhf(isymm)
C           
            DO C = 1,NVIR(ISYMC)
               DO J = 1,NUMIJ
                  IJ = IJ1 + J - 1
                  DO L = 1,NRHF(ISYML)
                     DO M = 1,NRHF(ISYMM) 
C                    
                        NMLJC = ISJIKA(ISYMMLJ,ISYMC)
     &                        + NMAIJK(ISYMMLJ)*(C - 1)
     &                        + IMAIJK(ISYMML,ISYMJ)
     &                        + NMATIJ(ISYMML)*(IJ - 1)
     &                        + IMATIJ(ISYMM,ISYML)
     &                        + NRHF(ISYMM)*(L - 1)
     &                        + M
                        NLMJC = ISJIKA(ISYMLMJ,ISYMC)
     &                        + NMAIJK(ISYMLMJ)*(C - 1)
     &                        + IMAIJK(ISYMLM,ISYMJ)
     &                        + NMATIJ(ISYMLM)*(IJ - 1)
     &                        + IMATIJ(ISYML,ISYMM)
     &                        + NRHF(ISYML)*(M - 1)
     &                        + L
                        NCLMJ = NCKI(ISYMCLM)*(J - 1)
     &                        + ICKI(ISYMCL,ISYMM)
     &                        + NT1AM(ISYMCL)*(M - 1)
     &                        + IT1AM(ISYMC,ISYML)
     &                        + NVIR(ISYMC)*(L - 1)
     &                        + C
C                       
                        XLINT(NCLMJ) = TWO*XOINT(NLMJC) - XOINT(NMLJC)
*                       WRITE(LUPRI,*) 'C=',C,'L=',L,'M=',M,'J=',J
*                       WRITE(LUPRI,*) 'NCLMJ = ',NCLMJ,' NLMJC = ',NLMJC, 
*    &                          ' NMLJC = ',NMLJC
C                        
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
C        
         ENDDO
      ENDDO
C     
      RETURN
      END
C  /* Deck ccho_dechoc2 */
      SUBROUTINE CCHO_DECHOC2(FOCKD,CHOELE,NUMCHO,ICHO,OCCHO,VICHO,
     &                       ISYDLI,ISYBCJ)
C
C     JLC, BFR, TBP, HK, AS, April 2003.
C
C     Construct vector for Cholesky decomposition:
C
C        For the occupied part, OCCHO(dli)  update vector
C        For the virtual  part, VICHO(bcj)  update vector
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),CHOELE(*),OCCHO(*),VICHO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      IF (ICHO .GT. NUMCHO) THEN
         WRITE(LUPRI,*) 'Only ',NUMCHO,' vectors were calculated',
     &              ' to converge the decomposition to 1.0D-13'
         CALL QUIT('Too many vectors required')
      END IF
C
C-------------------
C     Occupied part.
C-------------------
C
      IF (ICHO .EQ. 1) THEN
         IND = 0
         DO ISYMI = 1,NSYM
            ISYMDL = MULD2H(ISYMI,ISYDLI)
            DO I = 1,NRHF(ISYMI)
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  DO L = 1,NRHF(ISYML)
                     DO D = 1,NVIR(ISYMD)
C
                        KOFFI = IRHF(ISYMI) + I
                        KOFFL = IRHF(ISYML) + L
                        KOFFD = IVIR(ISYMD) + D
C
                        OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) = 
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
         NDIMOC = IND
      ELSE
         IND = 0
         DO ISYMI = 1,NSYM
            ISYMDL = MULD2H(ISYMI,ISYDLI)
            DO I = 1,NRHF(ISYMI)
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  DO L = 1,NRHF(ISYML)
                     DO D = 1,NVIR(ISYMD)
C
                        KOFFI = IRHF(ISYMI) + I
                        KOFFL = IRHF(ISYML) + L
                        KOFFD = IVIR(ISYMD) + D
C
                        OME = FOCKD(KOFFD)-FOCKD(KOFFL)-FOCKD(KOFFI)
C
                        IND = IND + 1
                        OCCHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
         NDIMOC = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMOC,FACTOR,OCCHO,1)
C
      ENDIF
C
C------------------
C     Virtual part.
C------------------
C
      IF (ICHO .EQ. 1) THEN
         IND = 0
         DO ISYMJ = 1,NSYM
            ISYMBC = MULD2H(ISYMJ,ISYBCJ)
            DO J = 1,NRHF(ISYMJ)
               DO ISYMC = 1,NSYM
                  ISYMB = MULD2H(ISYMC,ISYMBC)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
C
                        KOFFJ = IRHF(ISYMJ) + J
                        KOFFC = IVIR(ISYMC) + C
                        KOFFB = IVIR(ISYMB) + B
C
                        OME = FOCKD(KOFFB)+FOCKD(KOFFC)-FOCKD(KOFFJ)
C
                        IND = IND + 1
                        VICHO(IND) = 
     &                      SQRT(TWO*CHOELE(1))/(OME+CHOELE(1))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
         NDIMVI = IND
C
      ELSE
C
         IND = 0
         DO ISYMJ = 1,NSYM
            ISYMBC = MULD2H(ISYMJ,ISYBCJ)
            DO J = 1,NRHF(ISYMJ)
               DO ISYMC = 1,NSYM
                  ISYMB = MULD2H(ISYMC,ISYMBC)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
C
                        KOFFJ = IRHF(ISYMJ) + J
                        KOFFC = IVIR(ISYMC) + C
                        KOFFB = IVIR(ISYMB) + B
C
                        OME = FOCKD(KOFFB)+FOCKD(KOFFC)-FOCKD(KOFFJ)
C
                        IND = IND + 1
                        VICHO(IND) =
     &                       (OME-CHOELE(ICHO-1))/ (OME+CHOELE(ICHO))
C
                     END DO
                  END DO
               END DO
            END DO
         ENDDO
C
         NDIMVI = IND
         FACTOR = SQRT(CHOELE(ICHO)/CHOELE(ICHO-1))
C
         CALL DSCAL(NDIMVI,FACTOR,VICHO,1)
C         
      ENDIF
C
      RETURN
      END
C  /* Deck ccho_sclt22c2 */
      SUBROUTINE CCHO_SCLT22C2(T22,CHOV,ISYMJ,NUMIJ)
C
C     JLC, BFR, TBP, HK, AS, April 2003.
C
C     Scale T22(abk;#j) with cholesky update vector: 
C        d(abk)*T22(abk;#j)
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION T22(*),CHOV(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      ISYMABK = ISYMJ
      NABK = NCKASR(ISYMABK)
C     
      DO J = 1,NUMIJ
         DO IABK = 1,NABK
C        
            IOFF = NABK*(J - 1) + IABK
            T22(IOFF) = CHOV(IABK)*T22(IOFF)
C           
         ENDDO
      ENDDO
      RETURN
      END
C  /* Deck ccho_scllc2 */
      SUBROUTINE CCHO_SCLLC2(XLINT,CHOO,ISYMJ,NUMIJ)
C
C     JLC, BFR, TBP, HK, AS, April 2003.
C
C     Scale L(clm;#j) with cholesky update vector: 
C        d(clm)*L(clm;#j)
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XLINT(*),CHOO(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      ISYMCLM = ISYMJ
      NCLM = NCKI(ISYMCLM)
C     
      DO J = 1,NUMIJ
         DO ICLM = 1,NCLM
C        
            IOFF = NCLM*(J - 1) + ICLM
            XLINT(IOFF) = CHOO(ICLM)*XLINT(IOFF)
C           
         ENDDO
      ENDDO
      RETURN
      END
*     subroutine escribe(LUPRI,xiint,xjint,xoint,t2vo,fockd,nrhft,nvirt)
*        implicit none
*        integer lupri
*        integer nrhft,nvirt,i,j,k,m,a,b,e,d
*        double precision xiint(nvirt,nrhft,nvirt,nvirt)
*        double precision xjint(nvirt,nrhft,nvirt,nvirt)
*        double precision xoint(nrhft,nrhft,nrhft,nvirt)
*        double precision t2vo(nvirt,nrhft,nvirt,nrhft)
*        double precision fockd(*)
*        write(lupri,*) 'XIINT:'
*        do b = 1,nvirt
*           do d = 1,nvirt
*              do m = 1,nrhft
*                 do e = 1,nvirt
*                    write(lupri,'(4I2,A,D11.4)') e,m,d,b,':',
*    &                                              xiint(e,m,d,b)
*                 enddo
*              enddo
*           enddo
*        enddo
*        write(lupri,*) 'XJINT:'
*        do b = 1,nvirt
*           do d = 1,nvirt
*              do m = 1,nrhft
*                 do e = 1,nvirt
*                    write(lupri,'(4I2,A,D11.4)') e,m,d,b,':',
*    &                                              xjint(e,m,d,b)
*                 enddo
*              enddo
*           enddo
*        enddo
*        write(lupri,*) 'XOINT:'
*        do a = 1,nvirt
*           do k = 1,nrhft
*              do j = 1,nrhft
*                 do i = 1,nrhft
*                    write(lupri,'(4I2,A,D11.4)') i,j,k,a,':',
*    &                                              xoint(i,j,k,a)
*                 enddo
*              enddo
*           enddo
*        enddo
*        write(lupri,*) 'T2VO:'
*        do j = 1,nrhft
*           do i = 1,nrhft
*              do b = 1,nvirt
*                 do a = 1,nvirt
*                    write(lupri,'(4I2,A,D11.4)') a,b,i,j,':',
*    &                                              xoint(a,b,i,j)
*                 enddo
*              enddo
*           enddo
*        enddo
*        write(lupri,*) 'FOCKD:'
*        write(lupri,'(4D11.4)') (FOCKD(i),i=1,nrhft+nvirt)
*        return
*        end
