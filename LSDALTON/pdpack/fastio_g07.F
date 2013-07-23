#ifdef UNDEF
===========================================================================
/* Comdeck notes */
C*  AUTHOR:   Per-Olof Widmark                                         *
C*            Department of Theoretical Chemistry                      *
C*            University of Lund                                       *
C*                                                                     *
C*  WRITTEN:  82-12-01                                                 *
C*  FPS-VER:  84-12-06                                                 *
C*  IBM-VER:  86-11-17                                                 *
C*  88-12-02  Hans Joergen Aa. Jensen -- Alliant version
C*            also changed IREC = IDISK/4 + 1 to = IDISK + 1
C*  93-10-05 hjaaj: LDISK byte storage where LDISK = 8/IRAT
C*                  Cray: LDISK = 8, i*4 machines: LDISK = 4
C*  94-05-20 hjaaj: DATEMP unformatted (was formatted !)
C*                  DARMOV,SQRMOV close and delete (remove)
C*  94-08-23 hjaaj: recl in *4 words for IRIX .. LRECL always in bytes
C*  97-04-09 hjaaj: renamed GETADR to SQGETA to avoid conflict with
C*           system routine GETADR; renamed SETADR to SQSETA for
C*           symmetry. Renamed GETERR to FIOERR to avoid future conflict.
===========================================================================
#endif
C  /* Deck fx8io */
C***********************************************************************
C*                                                                     *
C*       F A S T I O                                                   *
C*       ***********                                                   *
C*                                                                     *
C*  THIS ROUTINE PERFORMS EITHER SYNCRONOUS I/O OR ASYNCRONOUS I/O     *
C*  USING THE FILE$ ROUTINES.                                          *
C*                                                                     *
C*  THE ROUTINE HAS THE FOLLOWING ENTRY POINTS:                        *
C*  FASTIO - ENTRY TO SET PARAMETERS                                   *
C*  DAFILE - DIRECT ACCESS I/O                                         *
C*  DAFILF - DIRECT ACCESS I/O WITH EOF CHECK                          *
C*  SQFILE - SEQUENTIAL I/O                                            *
C*  SQFILF - SEQUENTIAL I/O WITH EOF CHECK                             *
C*  DAWAIT - WAIT FOR COMPLETION OF ASYNCRONOUS DIRECT ACCESS I/O      *
C*  SQWAIT - WAIT FOR COMPLETION OF ASYNCRONOUS SEQUENTIAL I/O         *
C*  DADONE - CHECK IF ASYNCRONOUS DIRECT ACCESS I/O IS COMPLETE        *
C*  SQDONE - CHECK IF ASYNCRONOUS SEQUENTIAL I/O IS COMPLETE           *
C*  DANAME - TO GIVE A FILE AN ALIAS NAME                              *
C*  SQNAME - TO GIVE A FILE AN ALIAS NAME                              *
C*  DATEMP - CREATE A TEMPORARY FILE                                   *
C*  SQTEMP - CREATE A TEMPORARY FILE                                   *
C*  GSLIST - CONSTRUCTS A LIST FOR GATHER/SCATTER I/O                  *
C*  SQSETA - TO SET ADDRESS FOR SEQUENTIAL I/O                         *
C*  SQGETA - TO GET ADDRESS FROM SEQUENTIAL I/O                        *
C*  FIOERR - TO GET ERROR CODE WHEN NOT USING HARD ERRFIX              *
C*                                                                     *
C*  GLOBAL VARIABLES:                                                  *
C*  NAME    LOCICAL UNIT NUMBER IN ASCII CHARACTERS.                   *
C*  ADDR    CURRENT DISK ADDRESS OF EACH LOGICAL UNIT.                 *
C*  OPEN    FLAG TO SIGNIFY THAT THE FILE IS OPENED.                   *
C*  ID      CHARACTER STRING IDENTIFYING THE ENTRY POINT USED.         *
C*  ASYNC   FLAG FOR ASYNCRONOUS I/O FOR EACH LOGICAL UNIT.            *
C*  LUPRI   LOGICAL UNIT NUMBER OF PRINT UNIT.                         *
C*  ERRFIX  CHARACTER STRING SPECIFYING ERROR HANDLING MODE,           *
C*          CAN BE: 'HARD', 'FIRM', 'SOFT' OR 'NONE'.                  *
C*  NERR    ERROR COUNTER (IF ERRFIX='FIRM').                          *
C*  MXERR   MAXIMUM NUMBER OF ERRORS (IF ERRFIX='FIRM').               *
C*                                                                     *
C*  Formal Parameters:                                                 *
C*  PARM1  - CHARACTER                                                 *
C*  LU     - INTEGER                                                   *
C*  IFUN   - INTEGER                                                   *
C*  IBUF   - INTEGER                                                   *
C*  LBUF   - INTEGER                                                   *
C*  IDISK  - INTEGER                                                   *
C*  EOF    - LOGICAL                                                   *
C*  ITEST  - INTEGER                                                   *
C*  ECODE  - INTEGER                                                   *
C*  LIST   - INTEGER                                                   *
C*  N      - INTEGER                                                   *
C*                                                                     *
C*  AUTHOR:   Per-Olof Widmark                                         *
C*            Department of Theoretical Chemistry                      *
C*            University of Lund                                       *
C*                                                                     *
C*  WRITTEN:  82-12-01                                                 *
C*  FPS-VER:  84-12-06                                                 *
C*  IBM-VER:  86-11-17                                                 *
C*  Alliant:  88-12-02 Hans Joergen Aa. Jensen                         *
C*                     also changed IREC = IDISK/4 + 1 to = IDISK + 1  *
C*            93-10-05 hjaaj: 8 byte storage for Cray                  *
C*                            (defined from IRAT)                      *
C*                                                                     *
C***********************************************************************
      SUBROUTINE FASTIO(PARM1)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*(*) PARM1
      DIMENSION IBUF(*),LIST(*)
C
#include <iratdef.h>
      PARAMETER ( LDISK = 8 / IRAT )
C
      PARAMETER (MXAUX=230*1024)
      LOGICAL*1 AUXBUF(MXAUX),REFBUF(2)
      LOGICAL   EOF, OLDDX
C
C  added by Robert Berger on 21.08.1998
C
C  REFBUF must be static, otherwise we'll loose our reference point
C  on subsequent calls
C
      SAVE REFBUF
C
      PARAMETER (MAXLU = 99)
      CHARACTER ID*6,ERRFIX*4,TRACE*3
      CHARACTER NAME(MAXLU)*14
      DIMENSION ADDR(MAXLU),OPEN(MAXLU),ASYNC(MAXLU)
      DIMENSION IFWC(MAXLU),LRECL(MAXLU)
      DIMENSION MXADDR(MAXLU)
#include <priunit.h>
C
      DATA NERR,MXERR/0,10/
      DATA ERRFIX,TRACE/'HARD','OFF'/
      DATA MXADDR /MAXLU*0/
      DATA ADDR   /MAXLU*0/
      DATA OPEN   /MAXLU*0/
      DATA ASYNC  /MAXLU*0/
      DATA IFWC   /MAXLU*0/
      DATA LRECL  /MAXLU*4096/
      DATA NAME   /MAXLU*'              '/
C
      INCADR(LU1)=(LDISK*IFWC(LU1)+LRECL(LU1)-1) / LRECL(LU1)
C     ... IFWC counts in LDISK byte integers (INTEGER*LDISK).
C         8 byte when IRAT = 1 (Cray); 4 byte when IRAT = 2
C
      IF(TRACE.EQ.'ON') THEN
         WRITE(LUPRI,*) 'FASTIO ENTERED AT ENTRY FASTIO'
         WRITE(LUPRI,*) 'PARAMETER:',PARM1
      END IF
C
      IF(PARM1.EQ.'ERRFIX=HARD') THEN
         ERRFIX='HARD'
      ELSE IF(PARM1.EQ.'ERRFIX=FIRM') THEN
         ERRFIX='FIRM'
      ELSE IF(PARM1.EQ.'ERRFIX=SOFT') THEN
         ERRFIX='SOFT'
      ELSE IF(PARM1.EQ.'ERRFIX=NONE') THEN
         ERRFIX='NONE'
      ELSE IF(PARM1.EQ.'TRACE=ON') THEN
         TRACE='ON'
      ELSE IF(PARM1.EQ.'TRACE=OFF') THEN
         TRACE='OFF'
      ELSE IF(PARM1.EQ.'STATUS') THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) '***************************************'
         WRITE(LUPRI,*) '******  FILE STATUS FROM FASTIO  ******'
         WRITE(LUPRI,*) '***************************************'
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) '        LU  NAME             XTENT     LRECL'
         DO 700 I=1,MAXLU
            IF(OPEN(I).EQ.1) THEN
               NOEXT=MXADDR(I)
               NOKB=LRECL(I)/1024
               WRITE(LUPRI,'(1X,I10,1X,A14,I9,I8,A)')
     &            I,NAME(I),NOEXT,NOKB,'KB'
            END IF
700      CONTINUE
      ELSE
         GOTO 910
      END IF
      RETURN
C$PAGE
C***********************************************************************
C*                                                                     *
C*  DAFILE:                                                            *
C*                                                                     *
C*  ENTRY POINT FOR DIRECT ACCESS I/O                                  *
C*  FORMAL PARAMETERS:                                                 *
C*    LU     LOGICAL UNIT NUMBER                                       *
C*    IFUN   TYPE OF I/O                                               *
C*      0    reserve space for later write, or skip logical record     *
C*      1    SYNCRONOUS WRITE                                          *
C*      2    SYNCRONOUS READ                                           *
C*      3    SYNCRONOUS GATHER WRITE                                   *
C*      4    SYNCRONOUS SCATTER READ                                   *
C*      5    REWIND                                                    *
C*      6    ASYNCRONOUS WRITE                                         *
C*      7    ASYNCRONOUS READ                                          *
C*      8    ASYNCRONOUS GATHER WRITE                                  *
C*      9    ASYNCRONOUS SCATTER READ                                  *
C*     10    REWIND                                                    *
C*    IBUF   BUFFER START ADDRESS                                      *
C*    LBUF   LENGTH OF BUFFER IN WORDS                                 *
C*    IDISK  DISK ADDRESS                                              *
C*                                                                     *
C*  SYNCRONOUS I/O:                                                    *
C*  IDISK SPECIFIES THE DISK ADRESS WHERE THE I/O IS TO START.         *
C*  IT IS UPDATED ON BOTH WRITE AND READ, ENBLING THE USER TO          *
C*  USE THE DAFILE ENTRY AS A SEQUENTIAL I/O                           *
C*                                                                     *
C*  ASYNCRONOUS I/O:                                                   *
C*  IDISK SPECIFIES THE DISK ADDRESS WHERE THE I/O IS TO START.        *
C*  IT IS UPDATED WHEN THE I/O IS COMPLETE, ON A CALL TO DADONE,       *
C*  OR DAWAIT.                                                         *
C*                                                                     *
C*  GATHER/SCATTER-I/O:                                                *
C*  IBUF IS THE LIST OF ACCESS WORDS CREATED IN ENTRY POINT GSLIST,    *
C*  AND LBUF IS NOT USED.                                              *
C*                                                                     *
C*                                                                     *
C***********************************************************************
      ENTRY DAFILE(LU,IFUN,IBUF,LBUF,IDISK)
      IF(TRACE.EQ.'ON') THEN
         WRITE(LUPRI,*) 'FASTIO ENTERED BY ENTRY DAFILE'
         WRITE(LUPRI,'(A,4(1X,Z8))') ' PARAMETERS:',
     &                           LU,IFUN,LBUF,IDISK
      END IF
C
      ID='DAFILE'
      IF(IDISK.EQ.0) IDISK=IABS(IDISK)
      IF(IDISK.LT.0) GOTO 901
      GOTO 100
C***********************************************************************
C*                                                                     *
C*  DAFILF:                                                            *
C*                                                                     *
C*  ENTRY POINT FOR DIRECT ACCESS I/O WITH EOF CHECK.                  *
C*                                                                     *
C***********************************************************************
      ENTRY DAFILF(LU,IFUN,IBUF,LBUF,IDISK,EOF)
      IF(TRACE.EQ.'ON') THEN
         WRITE(LUPRI,*) 'FASTIO ENTERED BY ENTRY DAFILF'
         WRITE(LUPRI,'(A,4(1X,Z8))') ' PARAMETERS:',
     &                           LU,IFUN,LBUF,IDISK
      END IF
      ID='DAFILF'
      EOF=.FALSE.
      IF(IDISK.EQ.0) IDISK=IABS(IDISK)
      IF(IDISK.LT.0) GOTO 901
      GOTO 100
C***********************************************************************
C*                                                                     *
C*  SQFILE:                                                            *
C*                                                                     *
C*  ENTRY POINT FOR SEQUENTIAL I/O.                                    *
C*                                                                     *
C***********************************************************************
      ENTRY SQFILE(LU,IFUN,IBUF,LBUF)
      IF(TRACE.EQ.'ON') THEN
         WRITE(LUPRI,*) 'FASTIO ENTERED BY ENTRY SQFILE'
         WRITE(LUPRI,'(A,3(1X,Z8))') ' PARAMETERS:',
     &                           LU,IFUN,LBUF
      END IF
      ID='SQFILE'
      GOTO 100
C***********************************************************************
C*                                                                     *
C*  SQFILF:                                                            *
C*                                                                     *
C*  ENTRY POINT FOR SEQUENTIAL I/O WITH EOF CHECK.                     *
C*                                                                     *
C***********************************************************************
      ENTRY SQFILF(LU,IFUN,IBUF,LBUF,EOF)
      IF(TRACE.EQ.'ON') THEN
         WRITE(LUPRI,*) 'FASTIO ENTERED BY ENTRY SQFILF'
         WRITE(LUPRI,'(A,3(1X,Z8))') ' PARAMETERS:',
     &                           LU,IFUN,LBUF
      END IF
      ID='SQFILF'
      EOF=.FALSE.
      GOTO 100
C$PAGE
C***********************************************************************
C*                                                                     *
C*  HERE THE I/O IS PERFORMED.                                         *
C*                                                                     *
C***********************************************************************
100   CONTINUE
      ISTAT = 0
      IF(.NOT.(0.LE.LU.AND.LU.LE.MAXLU)) GOTO 903
      IF(.NOT.(0.LE.IFUN.AND.IFUN.LE.10)) GOTO 904
      IF(IFUN.EQ.0) THEN
C     ... reserve space for later write, or skip logical record
         IFWC(LU) = LBUF
C        ... previous line inserted by hjaaj 881202
         IF(ID(1:5).EQ.'DAFIL') ADDR(LU)=IDISK
         ADDR(LU)=ADDR(LU)+INCADR(LU)
         IF(ID(1:5).EQ.'DAFIL') IDISK=ADDR(LU)
         MXADDR(LU)=MAX(MXADDR(LU),ADDR(LU))
         RETURN
      END IF
      IF(MOD(IFUN,5).EQ.0) THEN
C     ... rewind
         ADDR(LU)=0
         IF(ID(1:5).EQ.'DAFIL') IDISK=0
         RETURN
      END IF
      IF(MOD(IFUN,5).LE.2.AND.LBUF.LE.0) GOTO 905
      IF(OPEN(LU).EQ.0) THEN
         LURECL = LRECL(LU)
         CALL GPOPEN(LU,' ','UNKNOWN','DIRECT','UNFORMATTED',LURECL,
     &               OLDDX)
         OPEN(LU)=1
      END IF
      IF(ASYNC(LU).EQ.0) THEN
         IF(ID(1:5).EQ.'DAFIL') ADDR(LU)=IDISK
      ELSE
         IF(IFUN.NE.ASYNC(LU)) GOTO 906
         WRITE(LUPRI,*) 'INTERNAL ERROR 001 IN FASTIO'
         CALL FABORT
         IF(ID(1:5).NE.'DAFIL') THEN
            ADDR(LU)=ADDR(LU)+INCADR(LU)
            MXADDR(LU)=MAX(MXADDR(LU),ADDR(LU))
         ELSE
            IF(ADDR(LU).EQ.IDISK) THEN
               ADDR(LU)=ADDR(LU)+INCADR(LU)
               MXADDR(LU)=MAX(MXADDR(LU),ADDR(LU))
               IDISK=ADDR(LU)
            ELSE
               ADDR(LU)=IDISK
            END IF
         END IF
      END IF
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      IF(IFUN.EQ.1 .OR. IFUN.EQ.6) THEN
C     ... syncronous/asyncronous write
         NREC=ADDR(LU)+1
         DO 1000 ISTART=1,LBUF,LRECL(LU)/LDISK
            IEND=MIN(LBUF,ISTART+LRECL(LU)/LDISK-1)
               CALL WRITDX(LU,NREC,IEND-ISTART+1,IBUF(ISTART))
            NREC=NREC+1
1000     CONTINUE
         IFWC(LU)=LBUF
      ELSE IF(IFUN.EQ.2 .OR. IFUN.EQ.7) THEN
C     ... syncronous/asyncronous read
         NREC=ADDR(LU)+1
         DO 1100 ISTART=1,LBUF,LRECL(LU)/LDISK
            IEND=MIN(LBUF,ISTART+LRECL(LU)/LDISK-1)
            CALL READDX(LU,NREC,IEND-ISTART+1,IBUF(ISTART))
            NREC=NREC+1
1100     CONTINUE
         IFWC(LU)=LBUF
      ELSE IF(IFUN.EQ.3 .OR. IFUN.EQ.8) THEN
C     ... syncronous/asyncronous gather write
         NTOT=0
         DO 1200 NBUF=1,IBUF(1)
            DO 1210 I=1,IBUF(2*NBUF+2)
               AUXBUF(I+NTOT)=REFBUF(IBUF(2*NBUF+1)+I)
1210        CONTINUE
            NTOT=NTOT+IBUF(2*NBUF+2)
1200     CONTINUE
         IF(NTOT.NE.IBUF(2)) THEN
            WRITE(LUPRI,*) 'INTERNAL ERROR 002 IN FASTIO'
            WRITE(LUPRI,*) 'OR EXTERNAL GATHER SCATTER LIST'
            WRITE(LUPRI,*) 'DESTROYED BY USER PROGRAM.'
            CALL FTNWB
            CALL FABORT
         END IF
         NREC=ADDR(LU)+1
         DO 1220 ISTART=1,NTOT,LRECL(LU)
            IEND=MIN(NTOT,ISTART+LRECL(LU)-1)
            CALL WRITDX(LU,NREC,IEND-ISTART,AUXBUF(ISTART))
            NREC=NREC+1
1220     CONTINUE
         IFWC(LU)=NTOT/LDISK
      ELSE IF(IFUN.EQ.4 .OR. IFUN.EQ.9) THEN
C     ... syncronous/asyncronous gather read
         NTOT=IBUF(2)
         NREC=ADDR(LU)+1
         DO 1300 ISTART=1,NTOT,LRECL(LU)
            IEND=MIN(NTOT,ISTART+LRECL(LU)-1)
            CALL READDX(LU,NREC,IEND-ISTART,AUXBUF(ISTART))
            NREC=NREC+1
1300     CONTINUE
         IFWC(LU)=NTOT/LDISK
         NTOT=0
         DO 1310 NBUF=1,IBUF(1)
            DO 1320 I=1,IBUF(2*NBUF+2)
               REFBUF(IBUF(2*NBUF+1)+I)=AUXBUF(I+NTOT)
1320        CONTINUE
            NTOT=NTOT+IBUF(2*NBUF+2)
1310     CONTINUE
         IF(NTOT.NE.IBUF(2)) THEN
            WRITE(LUPRI,*) 'INTERNAL ERROR 003 IN FASTIO'
            WRITE(LUPRI,*) 'OR EXTERNAL GATHER SCATTER LIST'
            WRITE(LUPRI,*) 'DESTROYED BY USER PROGRAM.'
            CALL FTNWB
            CALL FABORT
         END IF
      ELSE
         ISTAT=1
         GOTO 911
      END IF
      IF(ASYNC(LU).EQ.0) THEN
         ADDR(LU)=ADDR(LU)+INCADR(LU)
         MXADDR(LU)=MAX(MXADDR(LU),ADDR(LU))
         IF(ID(1:5).EQ.'DAFIL') IDISK=ADDR(LU)
      END IF
      IF(ISTAT.NE.0) GOTO 900
      RETURN
C$PAGE
C***********************************************************************
C*                                                                     *
C*  DAWAIT:                                                            *
C*                                                                     *
C*  ENTRY POINT FOR WAIT ON ASYNCRONOUS DIRECT ACCESS I/O.             *
C*                                                                     *
C***********************************************************************
      ENTRY DAWAIT(LU,IDISK)
      ID='DAWAIT'
      GOTO 200
C***********************************************************************
C*                                                                     *
C*  SQWAIT:                                                            *
C*                                                                     *
C*  ENTRY POINT FOR WAIT ON ASYNCRONOUS SEQUENTIAL I/O                 *
C*                                                                     *
C***********************************************************************
      ENTRY SQWAIT(LU)
      ID='SQWAIT'
      GOTO 200
C***********************************************************************
C*                                                                     *
C*  WAIT FOR COMPLETON OF I/O                                          *
C*                                                                     *
C***********************************************************************
200   CONTINUE
      IF(.NOT.(1.LE.LU.AND.LU.LE.MAXLU)) GOTO 903
      IF(ASYNC(LU).EQ.0) RETURN
      ADDR(LU)=ADDR(LU)+INCADR(LU)
      MXADDR(LU)=MAX(MXADDR(LU),ADDR(LU))
      IF(ID(1:5).EQ.'DAFIL') IDISK=ADDR(LU)
      ASYNC(LU)=0
      RETURN
C$PAGE
C***********************************************************************
C*                                                                     *
C*  DADONE:                                                            *
C*                                                                     *
C*  TESTING FOR COMPLETION OF ASYNCRONOUS DIRECT ACCESS I/O            *
C*                                                                     *
C***********************************************************************
      ENTRY DADONE(LU,ITEST,IDISK)
      ID='DADONE'
      GOTO 300
C***********************************************************************
C*                                                                     *
C*  SQDONE:                                                            *
C*                                                                     *
C*  TESTING FOR COMPLETION OF ASYNCRONOUS SEQUENTIAL I/O               *
C*                                                                     *
C***********************************************************************
      ENTRY SQDONE(LU,ITEST)
      ID='SQDONE'
      GOTO 300
C***********************************************************************
C*                                                                     *
C*  CHECK IF THE I/O IS COMPLETE.                                      *
C*  IF IT IS COMPLETE SET ITEST=1 AND JUMP TO WAIT SECTION FOR         *
C*  UPDATE OF ADDRESS.                                                 *
C*  IF NOT SET ITEST=0 AND RETURN.                                     *
C*                                                                     *
C*  AT PRESENT EQUIVALENT TO DAWAIT/SQWAIT                             *
C***********************************************************************
300   CONTINUE
      IF(.NOT.(1.LE.LU.AND.LU.LE.MAXLU)) GOTO 903
      IF(ASYNC(LU).EQ.0) RETURN
      ADDR(LU)=ADDR(LU)+INCADR(LU)
      MXADDR(LU)=MAX(MXADDR(LU),ADDR(LU))
      IF(ID.EQ.'DADONE') IDISK=ADDR(LU)
      ASYNC(LU)=0
      RETURN
C***********************************************************************
C*                                                                     *
C*  DANAME:                                                            *
C*                                                                     *
C*  ENTRY POINT TO GIVE LU AN ALIAS NAME.                              *
C*                                                                     *
C***********************************************************************
      ENTRY DANAME(LU,PARM1)
      IF(TRACE.EQ.'ON') THEN
         WRITE(LUPRI,*) 'FASTIO ENTERED AT ENTRY DANAME'
         WRITE(LUPRI,*) ' PARAMETERS: ',LU,',',PARM1
         CALL FLSHFO(LUPRI)
      END IF
      ID='DANAME'
      IF(.NOT.(1.LE.LU.AND.LU.LE.MAXLU)) GOTO 903
      IF(OPEN(LU).NE.0) GOTO 912
      NAME(LU)=PARM1
      LURECL = LRECL(LU)
      CALL GPOPEN(LU,NAME(LU),'UNKNOWN','DIRECT','UNFORMATTED',LURECL,
     *            OLDDX)
      OPEN(LU)=1
      RETURN
C***********************************************************************
C*                                                                     *
C*  SQNAME:                                                            *
C*                                                                     *
C*  ENTRY POINT TO GIVE LU AN ALIAS NAME.                              *
C*                                                                     *
C***********************************************************************
      ENTRY SQNAME(LU,PARM1)
      ID='SQNAME'
      IF(.NOT.(1.LE.LU.AND.LU.LE.MAXLU)) GOTO 903
      IF(OPEN(LU).NE.0) GOTO 912
      NAME(LU)=PARM1
      LURECL = LRECL(LU)
      CALL GPOPEN(LU,NAME(LU),'UNKNOWN','DIRECT','UNFORMATTED',LURECL,
     &            OLDDX)
      OPEN(LU)=1
      RETURN
C***********************************************************************
C*                                                                     *
C*  DARECL:                                                            *
C*                                                                     *
C*  ENTRY POINT TO ASSIGN BUFFER SIZE TO FILES                         *
C*                                                                     *
C***********************************************************************
      ENTRY DARECL(LU,LREC)
      ID='DARECL'
      IF(.NOT.(1.LE.LU.AND.LU.LE.MAXLU)) GOTO 903
      IF(OPEN(LU).NE.0) GOTO 915
      LRECL(LU)=1024*((LREC+1023)/1024)
      RETURN
C***********************************************************************
C*                                                                     *
C*  SQRECL:                                                            *
C*                                                                     *
C*  ENTRY POINT TO ASSIGN BUFFER SIZE TO FILES                         *
C*                                                                     *
C***********************************************************************
      ENTRY SQRECL(LU,LREC)
      ID='SQRECL'
      IF(.NOT.(1.LE.LU.AND.LU.LE.MAXLU)) GOTO 903
      IF(OPEN(LU).NE.0) GOTO 915
      LRECL(LU)=1024*((LREC+1023)/1024)
      RETURN
C***********************************************************************
C*                                                                     *
C*  DATEMP:                                                            *
C*                                                                     *
C*  ENTRY POINT TO OPEN A TEMPORARY FILE.                              *
C*                                                                     *
C***********************************************************************
      ENTRY DATEMP(LU,SIZE)
      ID='DATEMP'
400   IF(.NOT.(1.LE.LU.AND.LU.LE.MAXLU)) GOTO 903
      IF(OPEN(LU).NE.0) GOTO 913
      LURECL = LRECL(LU)
      CALL GPOPEN(LU,' ','UNKNOWN','DIRECT','UNFORMATTED',LURECL,OLDDX)
      RETURN
C***********************************************************************
C*                                                                     *
C*  SQTEMP:                                                            *
C*                                                                     *
C*  ENTRY POINT TO OPEN A TEMPORARY FILE.                              *
C*                                                                     *
C***********************************************************************
      ENTRY SQTEMP(LU,SIZE)
      ID='SQTEMP'
      GOTO 400
C***********************************************************************
C*                                                                     *
C*  DACLOS:                                                            *
C*                                                                     *
C*  ENTRY POINT TO CLOSE A LOGICAL UNIT.                               *
C*                                                                     *
C***********************************************************************
      ENTRY DACLOS(LU)
      ID='DACLOS'
      IF(OPEN(LU).EQ.0) GOTO 914
      OPEN(LU)=0
      CALL GPCLOSE(LU,'KEEP')
      RETURN
C***********************************************************************
C*                                                                     *
C*  SQCLOS:                                                            *
C*                                                                     *
C*  ENTRY POINT TO CLOSE A LOGICAL UNIT.                               *
C*                                                                     *
C***********************************************************************
      ENTRY SQCLOS(LU)
      ID='SQCLOS'
      IF(OPEN(LU).EQ.0) GOTO 914
      OPEN(LU)=0
      CALL GPCLOSE(LU,'KEEP')
      RETURN
C***********************************************************************
C*                                                                     *
C*  DARMOV:                                                            *
C*                                                                     *
C*  ENTRY POINT TO CLOSE and delete A LOGICAL UNIT.                    *
C*                                                                     *
C***********************************************************************
      ENTRY DARMOV(LU)
      ID='DARMOV'
      IF(OPEN(LU).EQ.0) GOTO 914
      OPEN(LU)=0
      CALL GPCLOSE(LU,'DELETE')
      RETURN
C***********************************************************************
C*                                                                     *
C*  SQRMOV:                                                            *
C*                                                                     *
C*  ENTRY POINT TO CLOSE and delete A LOGICAL UNIT.                    *
C*                                                                     *
C***********************************************************************
      ENTRY SQRMOV(LU)
      ID='SQRMOV'
      IF(OPEN(LU).EQ.0) GOTO 914
      OPEN(LU)=0
      CALL GPCLOSE(LU,'DELETE')
      RETURN
C$PAGE
C***********************************************************************
C*                                                                     *
C*  SQSETA:                                                            *
C*                                                                     *
C*  THE ADDRESS IS SET BY PARAMETER                                    *
C*                                                                     *
C***********************************************************************
      ENTRY SQSETA(LU,IDISK)
      ID='SQSETA'
      IF(.NOT.(1.LE.LU .AND. LU.LE.MAXLU)) GOTO 903
      IF(IDISK.EQ.0) IDISK=IABS(IDISK)
      IF(IDISK.LT.0) GOTO 901
      ADDR(LU)=IDISK
      RETURN
C***********************************************************************
C*                                                                     *
C*  SQGETA:                                                            *
C*                                                                     *
C*  THE ADDRESS IS TRANSFERRED TO THE PARAMETER IDISK                  *
C*                                                                     *
C***********************************************************************
      ENTRY SQGETA(LU,IDISK)
      ID='SQGETA'
      IF(.NOT.(1.LE.LU .AND. LU.LE.MAXLU)) GOTO 903
      IDISK=ADDR(LU)
      RETURN
C***********************************************************************
C*                                                                     *
C*  FIOERR:                                                            *
C*                                                                     *
C*  THE I/O ERROR CODE IS TRANSFERRED TO THE CALLING PROGRAM.          *
C*                                                                     *
C***********************************************************************
      ENTRY FIOERR(ECODE)
      ID='FIOERR'
      ECODE=ISTAT
      ISTAT=0
      RETURN
C$PAGE
C***********************************************************************
C*                                                                     *
C*  GSLIST:                                                            *
C*                                                                     *
C*  ENTRY POINT FOR CONSTRUCTION OF LIST OF ACCESS WORDS USED          *
C*  IN GATHER/SCATTER-I/O.                                             *
C*                                                                     *
C***********************************************************************
      ENTRY GSLIST(LIST,N,
     & V01,N01,V02,N02,V03,N03,V04,N04,V05,N05,
     & V06,N06,V07,N07,V08,N08,V09,N09,V10,N10,
     & V11,N11,V12,N12,V13,N13,V14,N14,V15,N15,
     & V16,N16,V17,N17,V18,N18,V19,N19,V20,N20)
C    & V21,N21,V22,N22,V23,N23,V24,N24,V25,N25,
C    & V26,N26,V27,N27,V28,N28,V29,N29,V30,N30,
C    & V31,N31,V32,N32,V33,N33,V34,N34,V35,N35,
C    & V36,N36,V37,N37,V38,N38,V39,N39,V40,N40,
C    & V41,N41,V42,N42,V43,N43,V44,N44,V45,N45,
C    & V46,N46,V47,N47,V48,N48,V49,N49,V50,N50)
      IF(N.LT.1.OR.N.GT.20) GOTO 908
      NTOT=0
      GOTO( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20)
     &     ,N
C    &     21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
C    &     41,42,43,44,45,46,47,48,49,50),N
20    LIST( 42)=LDISK*N20
      LIST( 41)=LOC(V20)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N20
      IF(N20.LE.0) GOTO 909
19    LIST( 40)=LDISK*N19
      LIST( 39)=LOC(V19)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N19
      IF(N19.LE.0) GOTO 909
18    LIST( 38)=LDISK*N18
      LIST( 37)=LOC(V18)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N18
      IF(N18.LE.0) GOTO 909
17    LIST( 36)=LDISK*N17
      LIST( 35)=LOC(V17)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N17
      IF(N17.LE.0) GOTO 909
16    LIST( 34)=LDISK*N16
      LIST( 33)=LOC(V16)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N16
      IF(N16.LE.0) GOTO 909
15    LIST( 32)=LDISK*N15
      LIST( 31)=LOC(V15)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N15
      IF(N15.LE.0) GOTO 909
14    LIST( 30)=LDISK*N14
      LIST( 29)=LOC(V14)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N14
      IF(N14.LE.0) GOTO 909
13    LIST( 28)=LDISK*N13
      LIST( 27)=LOC(V13)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N13
      IF(N13.LE.0) GOTO 909
12    LIST( 26)=LDISK*N12
      LIST( 25)=LOC(V12)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N12
      IF(N12.LE.0) GOTO 909
11    LIST( 24)=LDISK*N11
      LIST( 23)=LOC(V11)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N11
      IF(N11.LE.0) GOTO 909
10    LIST( 22)=LDISK*N10
      LIST( 21)=LOC(V10)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N10
      IF(N10.LE.0) GOTO 909
 9    LIST( 20)=LDISK*N09
      LIST( 19)=LOC(V09)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N09
      IF(N09.LE.0) GOTO 909
 8    LIST( 18)=LDISK*N08
      LIST( 17)=LOC(V08)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N08
      IF(N08.LE.0) GOTO 909
 7    LIST( 16)=LDISK*N07
      LIST( 15)=LOC(V07)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N07
      IF(N07.LE.0) GOTO 909
 6    LIST( 14)=LDISK*N06
      LIST( 13)=LOC(V06)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N06
      IF(N06.LE.0) GOTO 909
 5    LIST( 12)=LDISK*N05
      LIST( 11)=LOC(V05)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N05
      IF(N05.LE.0) GOTO 909
 4    LIST( 10)=LDISK*N04
      LIST(  9)=LOC(V04)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N04
      IF(N04.LE.0) GOTO 909
 3    LIST(  8)=LDISK*N03
      LIST(  7)=LOC(V03)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N03
      IF(N03.LE.0) GOTO 909
 2    LIST(  6)=LDISK*N02
      LIST(  5)=LOC(V02)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N02
      IF(N02.LE.0) GOTO 909
 1    LIST(  4)=LDISK*N01
      LIST(  3)=LOC(V01)-LOC(REFBUF)
      NTOT=NTOT+LDISK*N01
      IF(N01.LE.0) GOTO 909
      LIST(1)=N
      LIST(2)=NTOT
      IF(NTOT.GT.MXAUX) THEN
        WRITE(LUPRI,*) 'AUXILIARY BUFFER IF FASTIO TOO SMALL'
        WRITE(LUPRI,'(A,I6)') ' SIZE REQUIRED: ',NTOT
        III=MXAUX
        WRITE(LUPRI,'(A,I6)') ' SIZE AVAILABLE:',III
        CALL FTNWB
        CALL FABORT
      END IF
      RETURN
C$PAGE
C***********************************************************************
C*                                                                     *
C*  ERROR HANDLING SECTION                                             *
C*                                                                     *
C***********************************************************************
900   CONTINUE
#if defined (SYS_AIX) || defined (SYS_LINUX) || defined (SYS_IRIX) || defined (SYS_HAL) || defined (SYS_HPUX)
      ISTEOF = -1
#endif
#if defined (SYS_CRAY) || defined (SYS_T90)
C     eof value for Cray from Trygve Sep.93 /hjaaj
      ISTEOF = 2
#endif
#if !defined (SYS_AIX) && !defined (SYS_LINUX) && !defined (SYS_IRIX) && !defined (SYS_HAL) && !defined (SYS_HPUX) && !defined (SYS_CRAY) && !defined (SYS_T90)
      CALL QUIT('FASTIO ERROR: ISTEOF has not been '//
     &          'defined for this computer')
#endif
      IF(ISTAT.EQ.ISTEOF .AND. (ID.EQ.'DAFILF'.OR.ID.EQ.'SQFILF')) THEN
         EOF=.TRUE.
         RETURN
      END IF
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,'(A,I3,2A)') ' I/O-ERROR ON UNIT:',LU,'   ENTRY:',ID
      IF(ISTAT.EQ.ISTEOF) THEN
         WRITE(LUPRI,*) 'END OF FILE'
      ELSE
         WRITE(LUPRI,*) 'ERROR NOT DECODED'
      END IF
      WRITE(LUPRI,'(A,I10)') ' STATUS=',ISTAT
      IF(ERRFIX.EQ.'FIRM') GOTO 980
      GOTO 990
C
901   ISTAT=901
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,'(A,I3)') ' NEGATIVE DISK ADDRESS ON UNIT:',LU
      IF(ERRFIX.EQ.'FIRM') GOTO 980
      GOTO 990
C
903   ISTAT=903
      IF(ERRFIX.EQ.'NONE') GOTO 960
      WRITE(LUPRI,'(A,I6,A)') ' UNIT NUMBER:',LU,' OUT OF RANGE'
      GOTO 990
C
904   ISTAT=904
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,'(A,I6,2A)') ' I/O-FUNCTION:',IFUN,' OUT OF ',
     & 'RANGE'
      WRITE(LUPRI,'(A,I3)') ' UNIT:',LU
      IF(ERRFIX.EQ.'FIRM') GOTO 980
      GOTO 990
C
905   ISTAT=905
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,*) ' ATTEMPT TO READ/WRITE NULL RECORD.'
      WRITE(LUPRI,'(A,I3)') ' UNIT:',LU
      IF(ERRFIX.EQ.'FIRM') GOTO 980
      GOTO 990
C
906   ISTAT=906
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,'(2A,I3)') ' ILLEGAL STACKING OF ASYNC. I/O ',
     & 'ON UNIT:',LU
      IF(ERRFIX.EQ.'FIRM') GOTO 980
      GOTO 990
C
C
908   ISTAT=908
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,'(A,I6,2A)') ' NUMBER OF SUBRECORDS',N,' OUT OF',
     & ' RANGE'
      WRITE(LUPRI,'(A,I3)') ' UNIT:',LU
      IF(ERRFIX.EQ.'FIRM') GOTO 980
      GOTO 990
C
909   ISTAT=909
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,*) ' ATTEMPT TO READ/WRITE NULL SUBRECORD IN ',
     & 'GATHER/SCATTER I/O'
      WRITE(LUPRI,'(A,I3)') ' UNIT:',LU
      IF(ERRFIX.EQ.'FIRM') GOTO 980
      GOTO 990
C
910   ISTAT=910
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,*) 'UNRECOGNIZABLE PARAMETER TO FASTIO!'
      IF(ERRFIX.EQ.'FIRM') GOTO 980
      GOTO 990
C
911   WRITE(LUPRI,'(A,I3.3,A)') ' INTERNAL ERROR ',ISTAT,' IN FASTIO'
      CALL FABORT
C
912   ISTAT=912
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,*) 'ATTEMPT TO NAME AN OPENED FILE TO:'
      WRITE(LUPRI,*) PARM1
      WRITE(LUPRI,'(A,I3)') ' UNIT:',LU
      IF(ERRFIX.EQ.'FIRM') GOTO 980
      GOTO 990
C
913   ISTAT=913
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,*) 'ATTEMPT TO OPEN A TEMPORARY FILE IN OPENED UNIT'
      WRITE(LUPRI,'(A,I3)') ' UNIT:',LU
      GOTO 990
C
914   ISTAT=914
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,*) 'ATTEMPT TO CLOSE A FILE THAT HAS NOT BEEN OPENED'
      WRITE(LUPRI,'(A,I3)') ' UNIT:',LU
      GOTO 990
C
915   ISTAT=915
      IF(ERRFIX.EQ.'NONE') GOTO 960
      IF(ERRFIX.EQ.'SOFT') GOTO 970
      WRITE(LUPRI,*) 'ATTEMPT TO CHANGE RECORDLENGTH OF OPENED FILE'
      WRITE(LUPRI,'(A,I3)') ' UNIT:',LU
      GOTO 990
C
960   RETURN
C
970   LU=ISTAT
      RETURN
C
980   NERR=NERR+1
      LU=ISTAT
      IF(NERR.GT.MXERR) THEN
         WRITE(LUPRI,*) 'ABORT DUE TO ERROR COUNT'
         GOTO 990
      END IF
      RETURN
C
990   CONTINUE
      CALL FTNWB
      CALL FABORT
C
      END
C  /* Deck fabort */
      SUBROUTINE FABORT
#include <priunit.h>
      WRITE(LUPRI,*) 'FASTIO.FABORT: abort'
      CALL QUIT(' ')
      END
C  /* Deck ftnwb */
      SUBROUTINE FTNWB
#include <priunit.h>
      WRITE(LUPRI,*) 'FASTIO.FTNWB: WALKBACK NOT PERFORMED'
      RETURN
      END
