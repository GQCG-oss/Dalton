*======================================================================*
       SUBROUTINE CC_MCDINP1(WORD)
*----------------------------------------------------------------------*
*
*    Written by: Sonia Coriani and Poul Joergensen (fall 1997)
*    Purpose: read input for CC magnetic circular dichroism
*
*    if (WORD .eq '*CCMCD  ') read & process input and set defaults, 
*    else set only defaults 
*    Use A,B for second order moment, C for first order moment
*
*    Relaxed/PDBS operators, Sonia Coriani february 2000
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccmcdinf.h>
#include <ccsdinp.h>
#include <ccsections.h>

* local parameters:
      CHARACTER SECNAM*(9)
      PARAMETER (SECNAM='CC_MCDINP')
      CHARACTER*(20) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_MCDINP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER NTABLE
      PARAMETER (NTABLE = 10)

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif

      PARAMETER (ZERO = 0.0d00)

* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER*8 LABELA, LABELB, LABELC
      CHARACTER*70 LABHELP
      CHARACTER TABLE(NTABLE)*(7)

      LOGICAL LARLX, LBRLX, LCRLX, LRELAX
      INTEGER IJUMP, IJ, ITOT
      INTEGER IDA(6), IDB(6), IDC(6), IDIP(3), IANG(3)
      INTEGER IXSYM , IXST
* data
      DATA SET /.FALSE./
      DATA TABLE /'.OPERAT','.MCD   ','.MCDLAO','.PRINT ','.NO2N+1',
     &            '.SELSTA','.RELAXE','.UNRELA','.USEPL1','.XXXXXX'/
      DATA IDA / 1, 2, 3, 1, 2, 3 /
      DATA IDB / 2, 3, 1, 3, 1, 2 /
      DATA IDC / 3, 1, 2, 2, 3, 1 /
* external function:
      INTEGER INDPRPCC

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
      CCMCD     = .FALSE.
      NMCDST    = 0                 ! # MCD (final) states
      NMCDOPER  = 0                 ! # MCD triples
      SELMCDST  = .FALSE.           ! Select MCD fin. state (default)
      LUSE2N1   = .TRUE.            ! 2N+1 rule (Mbar^f vects in LR, default)
      LUSEPL1   = .FALSE.           ! debug use of Left transformed vectors
      IPRMCD    = 0                 ! Print level (default)

      LARLX        = .FALSE.        !Relaxed A operator
      LBRLX        = .FALSE.        !Relaxed B operator
      LCRLX        = .FALSE.        !Relaxed C operator
      LRELAX       = .FALSE.        !Relaxation 
*
      ICHANG    = 0                   
*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCMCD  ') THEN

100   CONTINUE
       
* get new input line:

        READ (LUCMD,'(A7)') WORD
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
        END DO

        IF (WORD(1:1) .EQ. '.') THEN
*         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
*         jump to the appropriate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_MCDINP.')
      
C           -----------------------------------------------------------
C           .OPERAT :  manually select triples of operator labels A,B,C
C                      A,B  for second order moments
C                      C    for first  order moment 
C           -----------------------------------------------------------
1           CONTINUE
              READ (LUCMD,'(3A)') LABELA, LABELB, LABELC
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).EQ.'!') THEN
                   CONTINUE
                ELSE IF (LABELA(1:1).EQ.'(') THEN
                    LARLX = .FALSE.
                    LBRLX = .FALSE.
                    LCRLX = .FALSE.
                    IF (LABELA(1:7).EQ.'(RELAX)') LARLX = .TRUE.
                    IF (LABELB(1:7).EQ.'(RELAX)') LBRLX = .TRUE.
                    IF (LABELC(1:7).EQ.'(RELAX)') LCRLX = .TRUE.
                    IF (LARLX .OR. LBRLX .OR. LCRLX) THEN
                      KEEPAOTWO = MAX(KEEPAOTWO,1)
                    END IF                                    
                ELSE 
                  IF (NMCDOPER.LT.MXMCDOP) THEN
                    NMCDOPER = NMCDOPER + 1
                    IAMCDOP(NMCDOPER) = INDPRPCC(LABELA)
                    IBMCDOP(NMCDOPER) = INDPRPCC(LABELB)
                    ICMCDOP(NMCDOPER) = INDPRPCC(LABELC)
                    LAMCDRX(NMCDOPER) = LARLX
                    LBMCDRX(NMCDOPER) = LBRLX
                    LCMCDRX(NMCDOPER) = LCRLX     
                  ELSE
                    WRITE(LUPRI,'(/2A,I5,/)') 
     &               ' NB. OF OPERATOR TRIPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXMCDOP
                       CALL QUIT('TOO MANY OPERATOR-TRIPLETS IN CCMCD.')
                  END IF
                END IF
                READ (LUCMD,'(3A)') LABELA, LABELB, LABELC
              END DO
              BACKSPACE(LUCMD)
            GO TO 100

C           -------------------------------------------------------
C           .MCD : calculate full tensor (r x L) * r = 6 components
C                  all operators UNRELAXED
C           -------------------------------------------------------
2           CONTINUE
              IF (NMCDOPER+6 .GT. MXMCDOP) THEN
                WRITE(LUPRI,'(2(/A,I5)/)') 
     &          ' NB. OF OPERATOR TRIPLES SPECIFIED  : ',NMCDOPER+6,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXMCDOP
                CALL QUIT('TOO MANY OPERATOR TRIPLES IN CCMCD.')

              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              IANG(1) = INDPRPCC('XANGMOM ')
              IANG(2) = INDPRPCC('YANGMOM ')
              IANG(3) = INDPRPCC('ZANGMOM ')
              DO IJ = 1,6
                IAMCDOP(IJ+NMCDOPER) = IDIP(IDA(IJ))
                IBMCDOP(IJ+NMCDOPER) = IANG(IDB(IJ))
                ICMCDOP(IJ+NMCDOPER) = IDIP(IDC(IJ))
                LAMCDRX(IJ+NMCDOPER) = LRELAX
                LBMCDRX(IJ+NMCDOPER) = LRELAX
                LCMCDRX(IJ+NMCDOPER) = LRELAX       
              END DO
              NMCDOPER = NMCDOPER + 6
            GO TO 100
C           -------------------------------------------------------
C           .MCDLAO : calculate full tensor (r x L) * r = 6 compnts
C                     L operator is dh/dB
C             UNFINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C           -------------------------------------------------------
3           CONTINUE
              IF (NMCDOPER+6 .GT. MXMCDOP) THEN
                WRITE(LUPRI,'(2(/A,I5)/)')
     &          ' NB. OF OPERATOR TRIPLES SPECIFIED  : ',NMCDOPER+6,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXMCDOP
                CALL QUIT('TOO MANY OPERATOR TRIPLES IN CCMCD.')

              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              IANG(1) = INDPRPCC('dh/dBX  ')
              IANG(2) = INDPRPCC('dh/dBY  ')
              IANG(3) = INDPRPCC('dh/dBZ  ')
              DO IJ = 1,6
                IAMCDOP(IJ+NMCDOPER) = IDIP(IDA(IJ))
                IBMCDOP(IJ+NMCDOPER) = IANG(IDB(IJ))
                ICMCDOP(IJ+NMCDOPER) = IDIP(IDC(IJ))
                LAMCDRX(IJ+NMCDOPER) = LRELAX
                LBMCDRX(IJ+NMCDOPER) = LRELAX
                LCMCDRX(IJ+NMCDOPER) = LRELAX
              END DO
              NMCDOPER = NMCDOPER + 6
            GO TO 100


C           -------------------------------------------------------
C           .PRINT : set desired print level (default = 0)
C           -------------------------------------------------------
4           CONTINUE
              READ (LUCMD,*) IPRMCD
            GO TO 100

C           ------------------------------------------------------
C           .NO2N+1 : don't use the 2N+1 rule, ie don't use Mbar^f
C                     for the calculation of the one-photon moment 
C                     for the C operator
C           ------------------------------------------------------
5           CONTINUE
              LUSE2N1 = .FALSE.
            GO TO 100

C           ---------------------------------------------------------------
C           .SELSTA : Select (final) states (Bfrequency zero by default)
C                     Specify then symmetry (IXSYM) and state number (IXST)
C                     of the state(s) we wish to calculate the transition 
C                     moments (one line with IXSYM,IXST for each state)
C           ---------------------------------------------------------------
C
6           CONTINUE

              SELMCDST =.TRUE.             

              READ (LUCMD,'(A70)') LABHELP        !read buffer line from input
              DO WHILE (LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
                IF (LABHELP(1:1).NE.'!') THEN
                  !read sym/number fin.state (from buffer line)
                  READ(LABHELP,*) IXSYM,IXST
                  IF (NMCDST.LT.MXMCDST) THEN
                    NMCDST = NMCDST + 1         !count how many
                    !put state-sym in array IMCDSTSY(*)
                    IMCDSTSY(NMCDST) = IXSYM
                    !put state-nr  in array IMCDSTNR(*)
                    IMCDSTNR(NMCDST) = IXST
                  ELSE
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ' ,MXMCDST 
                    CALL QUIT(
     &                'TOO MANY STATES SPECIFIED BY .SELSTA in MCD')
                  END IF
                END IF
                READ (LUCMD,'(A70)') LABHELP
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
C           ----------------------------------------------------------
C           .RELAXE : switch to relaxed modus for all three operators:
C           ----------------------------------------------------------
7           CONTINUE
              ! LRELAX    = .TRUE.
              ! KEEPAOTWO = MAX(KEEPAOTWO,1)
              WRITE (LUPRI,*) 
     *            '.RELAXE keyword in *CCMCD section is disabled.'
            GO TO 100
C           ------------------------------------------------------------
C           .UNRELA : switch to unrelaxed modus for all three operators:
C           ------------------------------------------------------------
8           CONTINUE
              LRELAX = .FALSE.
            GO TO 100                         
C           -----------------------------------------------------
C           .USEPL1 : use left transformed contributions (debug)
C           -----------------------------------------------------
9           CONTINUE
              LUSEPL1 = .TRUE.
              WRITE (LUPRI,*) SECNAM, 
     *            ': Use PL1 and left A transformations'
            GO TO 100
C           -----------------------
C           .XXXXXX : unused labels
C           -----------------------
10          CONTINUE
              WRITE (LUPRI,*)SECNAM,': unused .XXXXXX label... ignored'
            GO TO 100
C

          ELSE
            WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
            CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
            CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A,/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE

*---------------------------------------------------------------------*
* check if any triple of operator labels has been specified:
* if not, use default: complete unrelaxed 
* {dipole x angmom * dipole} tensor
*---------------------------------------------------------------------*
      IF (NMCDOPER .EQ. 0) THEN 
         IF (NMCDOPER+6 .GT. MXMCDOP) THEN
            WRITE(LUPRI,'(2(/A,I5)/)') 
     &      ' NB. OF OPERATOR TRIPLES SPECIFIED  : ',NMCDOPER+6,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MXMCDOP
            CALL QUIT('TO MANY OPERATOR TRIPLES IN CCMCD.')
         END IF
         IDIP(1) = INDPRPCC('XDIPLEN ')
         IDIP(2) = INDPRPCC('YDIPLEN ')
         IDIP(3) = INDPRPCC('ZDIPLEN ')
         IANG(1) = INDPRPCC('XANGMOM ')
         IANG(2) = INDPRPCC('YANGMOM ')
         IANG(3) = INDPRPCC('ZANGMOM ')
         DO ITOT=1,6
            IAMCDOP(ITOT+NMCDOPER) = IDIP(IDA(ITOT))
            IBMCDOP(ITOT+NMCDOPER) = IDIP(IDB(ITOT))
            ICMCDOP(ITOT+NMCDOPER) = IDIP(IDC(ITOT))
            LAMCDRX(ITOT+NMCDOPER) = LRELAX
            LBMCDRX(ITOT+NMCDOPER) = LRELAX
            LCMCDRX(ITOT+NMCDOPER) = LRELAX       
         END DO
         NMCDOPER = NMCDOPER + 6
      END IF
*---------------------------------------------------------------------*
* set CCMCD = TRUE if we are to calculate anything at all
*---------------------------------------------------------------------*
C      WRITE (LUPRI,*) 'SET CCMCD TO TRUE'
      CCMCD  = .TRUE.

      RETURN
      END
*---------------------------------------------------------------------*
