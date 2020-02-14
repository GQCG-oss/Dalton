c/* deck cc_lanczos_lrinp */
      SUBROUTINE CC_lanczos_lrinp(WORD)
C---------------------------------------------------------------------*
C
C    Purpose: Read input for damped linear response using 
C             Lanczos algorithm
C             WORD='*CCLRLA'(nczos) need to decide 
C
C    Sonia Coriani, Aug. 2010
C    first revision March 2012
C=====================================================================*
      IMPLICIT NONE

#include "priunit.h"
#include "ccsections.h"
#include "ccrspprp.h"
#include "cclrlancz.h"
#include "codata.h"
#include "cclrinf.h"

* local parameters:
      CHARACTER SECNAM*(13)
      PARAMETER (SECNAM='CC_LANCZOS_LRINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 15)
      INTEGER IJUMP,I, IFREQ
      INTEGER IDIP(3),IANG(3)

* variables:

      CHARACTER WORD*(7)
      CHARACTER LABHELP*(80), LABELA*(8), LABELB*(8)
      CHARACTER TABLE(NTABLE)*(8)
      LOGICAL LRELAX

#if defined (SYS_CRAY)
      REAL d0, d1
#else
      DOUBLE PRECISION d0, d1
#endif
      PARAMETER (D0=0.0D0, D1=1.0D0)

* data:
!      DATA SETGSTOPA /.FALSE./
!      DATA SETGSTTPA /.FALSE./
!      DATA SETXSTOPA /.FALSE./
      DATA TABLE / '.CHAINL','.OPERAT','.DIPLEN','.DIPVEL','.ANGMOM',
     &             '.PRINT ','.DAMPIN','.FREQ I','.JRESTA','.ANALYZ',
     &             '.EIGFIL','.ALLEIG','.SUMRUL','.REDMEM','.DBGSYM'/

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCLRLA') THEN
         write(lupri,*)'INPUT SESSION FOR LANCZOS LINEAR RESPONSE'
      ELSE
         CALL QUIT('CC_LANCZOS_LRINP 4 wrong section:'//WORD(1:7))
      END IF

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*

      JCHAIN  = 1
      JCHAINOLD   = 0
      JCHAINNEW   = 1
      LCHAINADD   = .FALSE.
      ABSANALYZE  = .FALSE.
      DUMP_EIGFIL = .FALSE.
      DUMP_ALLFIL = .FALSE.
      SUM_RULES   = .FALSE.
      !Warning
      REDMEML     = .true.
      Debug_sym   = .true.
C
C     Common default value of the damping parameter is set to 
C     be 1000 cm-1 = 4.556333D-3 a.u.
C
      LRELAX=.false.
      NDAMP = 1   !number of gammas
      DAMPING(NDAMP) = 1000/XTKAYS
      LABELO = 'XDIPLEN '

      FREQ_RANGE(1)= D0
      FREQ_RANGE(2)= D0
      FREQ_RANGE(3)= D1

      EIG_RANGE(1)= D0
      EIG_RANGE(2)= D0
 
*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*

100   CONTINUE

! get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
           READ (LUCMD,'(A7)') WORD
           CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
c         jump to the appropriate input section:
          IF (IJUMP .LE. NTABLE) THEN
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), IJUMP
            CALL QUIT('Illegal address GOTO in CC_LANCZOS_LRINP')

C           ---------------------------------------------------------
C           .CHAINLENGTH
C            Manually select the chain length J, i.e. nr of q and p
C            Lanczos vectors in Lanczos tridiagonalization procedure
C           ---------------------------------------------------------

1           CONTINUE
            READ (LUCMD,*) JCHAIN
            GO TO 100

C           ------------------------
C           .OPERAT: operator labels 
C           For the time being we limit ourselves to <<OPER,OPER>>_w+ig
C           ------------------------
2           CONTINUE
            READ (LUCMD,'(2A)') LABELA, LABELB
            DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (NLROP.LT.MXLROP) THEN
                    CALL CC_LRINPREQ(LABELA,LABELB,1,1,.TRUE.,LRELAX)
                    !hack for Lanczos drive
                    IF (LABELA(1:4).EQ.LABELB(1:4)) LABELO=LABELA
                    !end hack for lanczos drive
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF OPERATOR DOUBLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLROP
                CALL QUIT('TOO MANY OPER-DOUBLES CC_LANCZOS_LRINP')
                  END IF
                END IF
                READ (LUCMD,'(2A)') LABELA, LABELB
             END DO
             BACKSPACE(LUCMD)

            GO TO 100
 
C           -----------------------------------------------------
C           .DIPLEN: calculate complete dipole transition vectors
C                    in length gauge 
C           Not yet implemented!!!!
C           -----------------------------------------------------
3           CONTINUE
!            IDIP(1) = INDPRP_CC('XDIPLEN ')
!            IDIP(2) = INDPRP_CC('YDIPLEN ')
!            IDIP(3) = INDPRP_CC('ZDIPLEN ')
!            do i=1,3
!               iop_lancz_A(i) = idip(i)
!               iop_lancz_B(i) = idip(i)
!               lrlx_lancz_A(i) = .FALSE.
!               lrlx_lancz_B(i) = .FALSE.
!            end do
!            nlanczop = 3

            GO TO 100

C           -----------------------------------------------------
C           .DIPVEL: calculate complete dipole transition vectors
C                    in velocity gauge 
C           Not yet implemented!!!!
C           -----------------------------------------------------
4           CONTINUE
            GO TO 100

C           ------------------------------------------------------
C           .ANGMOM: calculate complete magnetic dipole transition 
C                    vectors and if possible rotatory strenghts
C           Not yet implemented!!!!
C           ------------------------------------------------------
5           CONTINUE
            GO TO 100

C           ------------
C           .PRINT 
C           ------------
6           CONTINUE
            READ (LUCMD,*) IPRLRLCZ
            GO TO 100
C           -------------------------------------------------------
C           .DAMPING: manually specify the damping factor if needed
C           specify NDAMP first (nr of damping gammas)
C           then the values on a row
C
C           -------------------------------------------------------
7           CONTINUE
            READ (LUCMD,*) NDAMP
            READ (LUCMD,*) (DAMPING(I),I=1,NDAMP)
            GO TO 100
C           -------------------------------------------------------
C           .FREQ INTERVAL: specify the frequency range of the calc
C           1=Fstart,2=Fstop,3=Fstep
C           -------------------------------------------------------
8           CONTINUE
            ABS_RANGE = .TRUE.
            READ(LUCMD,*) (FREQ_RANGE(I), I=1,3)
            GO TO 100
C           -------------------------------------------------------
C           .JRESTART: restart the calculation from the previous
C            Q vectors (for previous Chain length)
C            Specify old value (Jstart) and new value (Jend)
C           -------------------------------------------------------
9           CONTINUE
            LCHAINADD = .TRUE.
            READ(LUCMD,*) JCHAINOLD
            READ(LUCMD,*) JCHAINNEW
            JCHAIN = JCHAINNEW
            GO TO 100
C           -------------------------------------------------------
C           .ANALYZE : analyze the 'orbital' nature of the band
C           -------------------------------------------------------
10          CONTINUE
            ABSANALYZE=.true.
            GO TO 100
C           -------------------------------------------------------
C           .EIGFIL  : save the R and L Lanczos pseudoeigenvectors
C                      on file for restart of CC exci calculation 
C           Use together with ANALYZE option to make sure you compute
C           the pseudoeigenvectors in full-space 
C           Select the freq range of the eigenvectors to be dumped on
C           file
C           -------------------------------------------------------
11          CONTINUE
            ABSANALYZE  = .true.
            DUMP_EIGFIL = .true.
            READ(LUCMD,*) (EIG_RANGE(I), I=1,2)
            GO TO 100

12          CONTINUE
C           -------------------------------------------------------
C           .ALLEIG  : save ALL R and L Lanczos pseudoeigenvectors
C                      on file for restart of CC exci calculation 
C           to be used together with ANALYZE option  (JCHAIN in tot)
C           -------------------------------------------------------
            ABSANALYZE  = .true.
            DUMP_EIGFIL = .false.
            DUMP_ALLFIL = .true.
            GO TO 100

13          CONTINUE
C           -------------------------------------------------------
C           .SUMRULES: 
C           -------------------------------------------------------
            SUM_RULES  = .true.
            ABSANALYZE = .false.
            GO TO 100

14          CONTINUE
C           ------------------------------------------------------
C           .REDMEM: use lanczos_drv2 to reduced memory calc of RE and LE
C           DEBUG STUFF
C           ------------------------------------------------------
            REDMEML = .true.
            GO TO 100

15          CONTINUE
C           ------------------------------------------------------
C           .DBGSYM: debug symmetry
C           ------------------------------------------------------
            Debug_sym = .true.
            GO TO 100

16          CONTINUE
C           ------------------------------------------------------
C           .NEWDRV: use new lanczos drive with alternative F constr
C           ------------------------------------------------------
!            newdriver = .true.
            GO TO 100

17          CONTINUE
C           ------------------------------------------------------
C           .OLDDRV: use old lanczos drive with old F constr
C           ------------------------------------------------------
!            olddriver = .true.
            GO TO 100

18          CONTINUE
C           ------------------------------------------------------
C           .NONORM: no csi and eta norms available on restart file
C           ------------------------------------------------------
!            nonorms = .true.
            GO TO 100
C
19          CONTINUE
C           ------------------------------------------------------
C           .OLDFMAT use dot function in QR F mat build instead
C           of my own batching strategy
C           ------------------------------------------------------
!            oldfmat = .true.
            GO TO 100

20          CONTINUE
C           ------------------------------------------------------
C           .SNSEEDS use alternative seeds. 
C            First case is exc state dependent seeds
C           ------------------------------------------------------
!            Sn_seeds = .true.
!            READ(LUCMD,*) StateNumber
!            READ(LUCMD,*) StateSymmetry
!            GO TO 100
21          CONTINUE
C           ------------------------------------------------------
C           .USE_RG use RG routine to diagonalize the T matrix
C           instead of DGEEV
C           ------------------------------------------------------
!            USE_RG = .true.
            GO TO 100
C           ------------------------------------------------------
C           .CVSEPAration remove valence only excitations
C           ------------------------------------------------------
22          CONTINUE
!            LCVSLCZ = .true.
!            WRITE(LUPRI,*)'LANCZOS_INPUT: CVS approx requested'
!            !how many per symmetry
!            READ(LUCMD,*) (NRHFCORE(I),I=1,MSYM)
!            !which ones
!            DO I = 1, MSYM
!               READ(LUCMD,*) (IRHFCORE(J,I),J=1,NRHFCORE(I))
!            END DO
!            WRITE(LUPRI,*)'Requested number of core orbs per sym'
!            write(lupri,*) (NRHFCORE(I),I=1,MSYM)
!            WRITE(LUPRI,*)'Indices of requested core orbs'
!            DO I = 1, MSYM
!               write(LUpri,*) (IRHFCORE(J,I),J=1,NRHFCORE(I))
!            END DO
            GO TO 100
C
C----------------------------------------------------------------
C           .IONISATION  ionisation trick
C----------------------------------------------------------------
C
          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF
200   CONTINUE
      IF (ABS_RANGE) THEN
            IF( FREQ_RANGE(1).GT.FREQ_RANGE(2) .OR.
     &         (FREQ_RANGE(2)-FREQ_RANGE(1)).LT.FREQ_RANGE(3)
     &         .OR. FREQ_RANGE(3).LE.0.0D0 ) THEN
               WRITE(LUPRI,'(/,A,/,A,/,/,A,/,A,/,/,A)')
     & ' Warning: Frequency range for damped calculation is not',
     & '          correctly specified according to:',
     & '                 .FREQ_RANGE',
     & '                 STARTFREQ ENDFREQ STEPSIZE',
     & '          This input will be ignored. Calculation continues.'
               ABS_RANGE = .FALSE.
               FREQ_RANGE(1)=D0
               FREQ_RANGE(2)=D0
               FREQ_RANGE(3)=D1
            END IF
      END IF
C-----------------------------------------------
C check, if operators and frequencies specified:
C-----------------------------------------------
C
!      IF (NBLRFR.EQ.0) THEN
!         NBLRFR   = 1
!         BLRFR(1) = 0.0D0
!      ENDIF
C
!      IF (ICHANG .NE. 0) THEN
!        IF (NLROP .EQ.0) WRITE(LUPRI,'(/A)')
!     &     '(*CCLR   input ignored, because no operators requested.)'
!      END IF
C----------------------------
C     Make wa frequency list. 
C----------------------------
!      DO IFREQ = 1, NBLRFR
!        ALRFR(IFREQ) = - BLRFR(IFREQ)
!      END DO
C

*---------------------------------------------------------------------*
* set CCLRLCZ flag and return:
*---------------------------------------------------------------------*
      CCLRLCZ  = .TRUE.

      RETURN
      END
C=====================================================================*
C                    END OF SUBROUTINE CC_LANCZOS_LRINP
C=====================================================================*
