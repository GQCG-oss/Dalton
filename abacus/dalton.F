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
#include <single.h>
C
#ifdef UNDEF
===========================================================================
/* Comdeck log */
010131-vebjornb: Old dalton.F split into dalton.F (program) and
                 dalgnr.F (subroutines).
951129-hjaaj: also CALL EXESOR in geo.opt. when RUNINTS
              force NEWTRA true if .INTSORT
              replaced ITMPMB with IMORDB in inpbak.h (IMORDB was not saved)
              removed CTMPMB in inpbak.h
              removed ITMPMB and CTMPMB references in SSAVE and SRESTR
              s/H1VORB/FCVORB/ and s/H1VRBB/FCVRBB/
951009-kenneth : Changed name from PSISTAR to DALTON
950621-vebjornb: First working optimization, GRAND changed to PSISTAR
950516-vebjornb: First effort to link all programs: GRAND
===========================================================================
#endif
C  /* Deck dalton */
      PROGRAM DALTON
C
C     NOTE
C
C     This is an experimental code for the evaluation of SCF, CC and
C     MCSCF molecular properties as energy derivatives or response
C     functions.  The authors accept no responsibility for the
C     performance of the code or for the correctness of the results.
C
C     The code (in whole or part) is not to be reproduced for
C     further distribution without the written permission of the authors.
C
C     If results obtained with this code are published, an appropriate
C     citation would be:
C
C     "Dalton, a molecular electronic structure program, Release 1.2
C     (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C     J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C     O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C     T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C     H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C     H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C     T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C     B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C     and O. Vahtras"
C
C     This note is not to be removed or altered.
C
#include <implicit.h>
#include <dummy.h>
#include <mtags.h>
#include <maxorb.h>
#include <infpar.h>
#include <priunit.h>
#include <infpri.h>
#include <gnrinf.h>
#include <inftap.h>
      EXTERNAL EXEDRV
#if defined (VAR_MPI)
      EXTERNAL NODSTR
      INCLUDE 'mpif.h'
      LOGICAL FINISH
#endif
      CHARACTER*12 FILNAM
      PARAMETER (LWORK = INSTALL_WRKMEM)
#if !defined (SYS_CONVEX) && !defined (SYS_HAL) 
      PARAMETER (IWOFF = 0)
#endif
#if defined (SYS_AIX) || defined (SYS_CRAY) || defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_CONVEX) || defined (SYS_T3D) || defined (SYS_SUN) || defined (SYS_LINUX) || defined (SYS_HPUX) || defined (SYS_HAL) || defined (SYS_T90) || defined (SYS_DEC)
      CHARACTER*20 WRKMEM
#endif
#if defined (SYS_CRAY) || defined (SYS_IRIX) || defined (SYS_T3D) || defined (SYS_SUN) || defined (SYS_HPUX) || defined (SYS_HAL) || defined (SYS_T90) || defined (SYS_DEC) || (defined (SYS_LINUX) && defined (VAR_PGF77))
      POINTER (WADDR,WORK(1))
#else
#if defined (SYS_CONVEX) || (defined (SYS_LINUX) && defined (VAR_G77))
      DIMENSION WORK(2)
#else
#if defined (SYS_NEC) || defined (SYS_AIX) || (defined (SYS_LINUX) && defined (VAR_IFC))
#if defined (NO_F90)
      POINTER (WADDR,WORK(1))
#else
      DOUBLE PRECISION , ALLOCATABLE , DIMENSION(:) :: WORK
      INTEGER                                       :: WADDR
#endif
#else
C
C     Static memory allocation
C
      DIMENSION WORK(LWORK+2)
#endif
#endif
#endif
      DATA WRKDLM/8H*WRKDLM*/
C
#if defined (VAR_MPI)
      CALL MPIXINIT
#else
      MYNUM = 0
#endif
      CALL GPIOIN
      IF (MYNUM .EQ. 0) THEN
      CALL GPOPEN(LUERR,'DALTON.ERR','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      ELSE
         LUERR = 0
      END IF
      LMWORK = 0
#if defined (SYS_AIX) || defined (SYS_CRAY) || defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_CONVEX) || defined (SYS_T3D) || defined (SYS_SUN) || defined (SYS_LINUX) || defined (SYS_HPUX) || defined (SYS_HAL) || defined (SYS_T90) || defined (SYS_DEC)
C
C     Dynamic memory allocation; check for WRKMEM in environment
C     The T3D thingie might work for DEC_ALPHA as well? Try it out, kr-dec96
C
      WRKMEM = '                    '
#if defined (SYS_T3D) || defined (SYS_T90)
      CALL PXFGETENV('WRKMEM',6,WRKMEM,LWRKMM,IERR)
#else
      CALL GETENV('WRKMEM',WRKMEM)
#endif
      READ(WRKMEM, '(I20)', ERR=10) LMWORK
      GOTO 20
   10    WRITE(LUERR,'(/3A)')
     &   ' DALTON: WRKMEM conversion error; WRKMEM = "',WRKMEM,'"'
         LMWORK = 0
   20 CONTINUE
#if defined (SYS_T3D) || defined (SYS_T90)
      CALL PXFGETENV ('BASDIR',6,BASDIR,LBASDR,IERR)
#else
      CALL GETENV ('BASDIR',BASDIR)
#endif
#endif
C
C     Either default memory size,
C
      IF (LMWORK .LE. 0) THEN
         WRITE (LUERR,'(/A)') ' DALTON: default work memory size used.'
         LMWORK = LWORK
      ELSE
C
C     ... or user specified memory size.
C
         WRITE(LUERR,'(/A/3A)')
     &      ' DALTON: user specified work memory size used,',
     &      '         environment variable WRKMEM = "',WRKMEM,'"'
      END IF
      NBYTES = (LMWORK+2) * 8
      XMBYTES = NBYTES
      XMBYTES = XMBYTES / (1024*1024)
      WRITE(LUERR,'(/A,I12,A,F8.2,A)') ' Work memory size (LMWORK):',
     &   LMWORK,' =',XMBYTES,' megabytes.'
C
C    Dynamic memory allocation
C
 1000 FORMAT(/' >>> DALTON: requested memory could not be allocated.'
     &      ,/'     Probable causes:',
     &       /'       (1) Insufficient free memory+page+swap space',
     &       /'       (2) "datasize" resource limit too small.')
C
#if defined (SYS_IRIX) || defined (SYS_AIX) || defined (SYS_NEC) || defined (SYS_SUN) || defined (SYS_HPUX) || defined (SYS_LINUX) || defined (SYS_HAL) || defined (SYS_DEC)
#if defined (SYS_IRIX) || defined (SYS_SUN) || defined (SYS_HPUX) || defined (SYS_HAL) || defined (SYS_DEC) || (defined (SYS_LINUX) && defined (VAR_PGF77))||defined(VAR_ABSOFT)
      WADDR = MALLOC(NBYTES)
#endif

#if defined (SYS_NEC) || defined (SYS_AIX)||(defined (SYS_LINUX) && defined (VAR_IFC))
#if defined (NO_F90)
      WADDR = MALLOC(NBYTES)
#else
      WADDR = 1
      ALLOCATE (WORK(LMWORK+2))
#endif
#endif
#if (defined (SYS_LINUX) && defined (VAR_G77))
      WADDR = 1
#endif
      IF (WADDR .EQ. 0) THEN
         WRITE(LUERR,1000)
         IF (LUERR .NE. 0) WRITE(0,1000)
         CALL QUIT(' DALTON: memory allocation failure')
      END IF
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
C
      IABORT = 0
      CALL HPALLOC(WADDR,LMWORK+2,IERR,IABORT)
      IF (IERR .LT. 0) THEN
         WRITE(LUERR,*) ' >>> Cray HPALLOC failure; IERR   = ', IERR
         IF (LUERR .NE. 0)
     &   WRITE(0,*) ' >>> Cray HPALLOC failure; IERR   = ', IERR
         IF (IERR .EQ. -2) THEN
            WRITE(LUERR,1000)
            IF (LUERR .NE. 0) WRITE(0,1000)
         END IF
         CALL QUIT(' DALTON: Core HPALLOC failure')
      END IF
#endif
C
      IF (BASDIR(1:1) .NE. '/') THEN
C
#if defined (VAR_ABSOFT)
         BASDIR="Insert basis-directory path here"
C        Example:
C        BASDIR="/Users/helgaker/Dalton/basis/"
#else
         BASDIR=
     &INSTALL_BASDIR
#endif
         WRITE(LUERR,'(/A,/8X,A)') ' Default basis set library used:',
     &                              BASDIR
      ELSE
         WRITE(LUERR,'(/A,/8X,A)') ' User supplied basis '//
     &                               'set directory :',BASDIR
      END IF
      LENBAS = 0
      DO 43 I = 1, 60
         IF (BASDIR(I:I) .EQ. ' ') GO TO 44
         LENBAS = LENBAS + 1
 43   CONTINUE
 44   CONTINUE
      IF (BASDIR(LENBAS:LENBAS) .NE. '/') THEN
         LENBAS = LENBAS + 1
         BASDIR(LENBAS:LENBAS) = '/'
      END IF
C
C     Set memory traps and initialize unit numbers
C
#if !(defined (SYS_LINUX) && defined (VAR_G77))
      WORK(1)        = WRKDLM
      WORK(2+LMWORK) = WRKDLM
#endif
C
C     Next 20+ lines commented out, does not work for general parallel
C     calculations.
C
c#if defined (VAR_MPI)
c      CALL MPIXINIT
c#if defined (SYS_CRAY)
c      CALL INFO = ISHELL('mkdir $SCRATCHDIR')
c      CALL INFO = ISHELL('cd $SCRATCHDIR')
c#else
c      CALL SYSTEM('mkdir $SCRATCHDIR')
c      CALL SYSTEM('cd $SCRATCHDIR')
c#endif
C
C     Allowing now for parallelization over different molecular geometries,
C     we must first check if this calculation indeed is of that kind.
C     We do this by letting the Master node read DALTON.INP, searching for
C     the keyword .PARNMD. If found, read DALTON/MOLECULE.INP and send them
C     to the slaves, who writes them out in their own files
C
c      IF (NODTOT .GT. 0) THEN
c         IF (MYNUM .EQ. 0) THEN
c            CALL PARIOT
c         ELSE
c            CALL PARION
c         END IF
c      END IF
c      CALL GPOPEN(LUPRI,'DALTON.OUT',' ',' ','FORMATTED',IDUMMY,.FALSE.)
c      CALL GPOPEN(LUCME,'DALTON.CM',' ',' ','FORMATTED',IDUMMY,.FALSE.)
c      CALL MEMINI(LUPRI,LUERR)
c#else
C      CALL GPOPEN(LUPRI,'DALTON.OUT',' ',' ','FORMATTED',IDUMMY,.FALSE.)
C      CALL GPOPEN(LUCME,'DALTON.CM',' ',' ','FORMATTED',IDUMMY,.FALSE.)
C      CALL MEMINI(LUPRI,LUERR)
C
C     This is not an MPI version, but we anyway set
C     variable MYNUM = 0 which for MPI identifies
C     this node as the master node for output below
C
C      MYNUM = 0
C#endif
C
C     Back to old version...   vb jan-01
C     Revised June 01 /hjaaj
C
      IF (MYNUM .EQ. 0) THEN
         CALL GPOPEN(LUPRI,'DALTON.OUT','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
         CALL GPOPEN(LUCME,'DALTON.CM','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      ELSE
         WRITE (FILNAM,'(A,I2.2)') 'DALTON.OUT',MYNUM
         CALL GPOPEN(LUPRI,FILNAM,'NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      END IF
      CALL MEMINI(LUPRI,LUERR)
C
C     *** End of old version ***
C
      IF (MYNUM .EQ. 0) THEN
C         
         CALL TITLER(' DALTON - An electronic structure program ',
     &               '*',111)
C
C        Changed format (WK/UniKA/04-11-2002).
         WRITE (LUPRI,'(T10,A,//,T26,A,/,6(/T10,A),//,T26,A,//,
     &        (T10,A))')
     * ' This is output from DALTON (Release 2.0, ? 2001)',
     *  ' Principal authors:',
     *  ' Trygve Helgaker,       University of Oslo,           Norway ',
     *  ' Hans Joergen Jensen,   SDU - Odense University,      Denmark',
     *  ' Poul Joergensen,       Aarhus University,            Denmark',
     *  ' Jeppe Olsen,           Aarhus University,            Denmark',
     *  ' Kenneth Ruud,          University of Tromsoe,        Norway ',
     *  ' Hans Aagren,           KTH Stockholm,                Sweden ',
     *  ' Contributors:',
     *  ' Alexander A. Auer,     University of Mainz,          Germany',
     *  ' Keld L. Bak,           UNI-C,                        Denmark',
     *  ' Vebjoern Bakken,       University of Oslo,           Norway ',
     *  ' Ove Christiansen,      University of Lund,           Sweden ',
     *  ' Sonia Coriani,         University of Trieste,        Italy  ',
     *  ' Paal Dahle,            University of Oslo,           Norway ',
     *  ' Erik K. Dalskov,       UNI-C,                        Denmark',
     *  ' Thomas Enevoldsen,     SDU - Odense University,      Denmark',
     *  ' Berta Fernandez,       U. of Santiago de Compostela, Spain  ',
     *  ' Christof Haettig,      Forschungszentrum Karlsruhe,  Germany',
     *  ' Kasper Hald,           Aarhus University,            Denmark',
     *  ' Asger Halkier,         Aarhus University,            Denmark',
     *  ' Hanne Heiberg,         University of Oslo,           Norway ',
     *  ' Hinne Hettema,         University of Auckland,       NZ     ',
     *  ' Dan Jonsson,           KTH Stockholm                 Sweden ',
     *  ' Sheela Kirpekar,       SDU - Odense University,      Denmark',
     *  ' Wim Klopper,           University of Karlsruhe,      Germany',
     *  ' Rika Kobayashi,        Cambridge University,         England',
     *  ' Henrik Koch,           University of Trondheim,      Norway ',
     *  ' Kurt V. Mikkelsen,     University of Copenhagen,     Denmark',
     *  ' Patrick Norman,        University of Linkoeping,     Sweden ',
     *  ' Thomas B. Pedersen,    U. of Santiago de Compostela, Spain  ',
     *  ' Martin J. Packer,      University of Sheffield,      UK     ',
     *  ' Torgeir A. Ruden,      University of Oslo,           Norway ',
     *  ' Pawel Salek,           KTH Stockholm,                Sweden ',
     *  ' Alfredo Sanchez,       University of Valencia,       Spain  ',
     *  ' Trond Saue,            University of Toulouse,       France ',
     *  ' Stephan P. A. Sauer,   University of Copenhagen,     Denmark',
     *  ' Bernd Schimmelpfennig, KTH Stockholm,                Sweden ',
     *  ' K. O. Sylvester-Hvid,  University of Copenhagen,     Denmark',
     *  ' Peter R. Taylor,       University of Warwick,        UK     ',
     *  ' Olav Vahtras,          PDC, Stockholm,               Sweden '
C
C        Changed format (WK/UniKA/04-11-2002).
         WRITE (LUPRI,'(/1X,69A1/)') ('-',I=1,69)
         WRITE (LUPRI,'(/,(5X,A))')
     &'NOTE:',
     &' ',
     &'This is an experimental code for the evaluation of molecular',
     &'properties using (MC)SCF wave functions. The authors accept no',
     &'responsibility for the performance of the code or for the',
     &'correctness of the results.',
     &' ',
     &'The code (in whole or part) is provided under a license and',
     &'is not to be reproduced for further distribution without',
     &'the written permission of the authors or their representatives.',
     &'See licence agreement for further information.',
     &' ',
     &'If results obtained with this code are published, an',
     &'appropriate citation would be:',
     &' ',
     &'"Dalton, a molecular electronic structure program, Release 1.2',
     &'(2001), written by T. Helgaker, H. J. Aa. Jensen, ' //
     &'P. Joergensen,',
     &'J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,',
     &'O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,',
     &'T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,',
     &'H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, W. Klopper,',
     &'R. Kobayashi, H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,'
     &,'T. B. Pedersen, T. A. Ruden, P. Salek, A. Sanchez, T. Saue, ',
     &'S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, ',
     &'P. R. Taylor and O. Vahtras"'
C
C     Stamp date and time and hostname to output
C
         CALL TSTAMP(' ',LUPRI)
C
#if defined (VAR_MPI)
         CALL GNRLIN
         IF (TESTIN) THEN
            CALL QUIT('*** End of input test for GENERAL ***')
         END IF
         CALL RECVNAMES
#if defined (VAR_G77)
         CALL CEXE(EXEDRV, LMWORK, WRKDLM,MASTER,MYNUM)
#else
         CALL EXEDRV(WORK, LMWORK, WRKDLM,MASTER,MYNUM)
#endif
         CALL FLSHFO(LUPRI)
C
C     Send ending signal to slaves
C
         IPRTYP = 0
         CALL MPIXBCAST(IPRTYP,1,'INTEGER',MYNUM)
         CALL GPCLOSE(LUCME,'KEEP')
         CALL GPCLOSE(LUERR,'KEEP')
         CALL GPCLOSE(LUMOLDEN,'KEEP')
         CALL GPCLOSE(LUPRI,'KEEP')
      ELSE
         IF (PARIO) THEN
            CALL GNRLIN
            IF (TESTIN) THEN
               CALL QUIT('*** End of input test for GENERAL ***')
            END IF
#if defined (VAR_G77)
            CALL CEXE(EXEDRV, LMWORK, WRKDLM,MASTER,MYNUM)
#else
            CALL EXEDRV(WORK, LMWORK,WRKDLM,MASTER,MYNUM)
#endif
C
C     Send ending signal to slaves
C
            FINISH = .TRUE.
            CALL MPI_BCAST(FINISH,1,MPI_LOGICAL,NODEID(MASTER),
     &           MPI_COMM_WORLD,IERR)
            CALL GPCLOSE(LUCME,'KEEP')
            CALL GPCLOSE(LUERR,'KEEP')
            CALL GPCLOSE(LUMOLDEN,'KEEP')
            CALL GPCLOSE(LUPRI,'KEEP')
         ELSE
            CALL SENDNAME(MASTER,MYNUM)
#if defined (VAR_G77)
            CALL CEXE(NODSTR,LMWORK,WRKDLM,MASTER,MYNUM)
#else
            CALL NODSTR(WORK,LMWORK,WRKDLM,MASTER,MYNUM)
#endif
         END IF
      END IF
      CALL MPIXFINALIZE
#else
      END IF
C
C     All below is not relevant for MPI as all proc. has stopped
C     with CALL MPI_FINALIZE(IERR)
C
      CALL GNRLIN
      IF (TESTIN) THEN
         CALL QUIT('*** End of input test for GENERAL ***')
      END IF
#endif
#if defined (VAR_PVM)
C
C     Spawn all PVM processes, as we now know how many there are
C
      IF (PARCAL) THEN
         CALL PVMXMYTID(MYTID)
         CALL PVMXSPAWN('node.x',NODTOT,NODEID(1))
         CALL RECVNAMES
         NODEID(0) = MYTID
      END IF
#endif
#if !defined (VAR_MPI)
#if defined (VAR_G77)
      CALL CEXE(EXEDRV, LMWORK, WRKDLM,MASTER,MYNUM)
#else
      CALL EXEDRV(WORK, LMWORK,WRKDLM,MASTER,MYNUM)
#endif
      CALL GPCLOSE(LUERR,'KEEP')
      CALL GPCLOSE(LUCME,'KEEP')
      CALL GPCLOSE(LUMOLDEN,'KEEP')
      CALL GPCLOSE(LUPRI,'KEEP')
#endif
#if defined (VAR_PVM)
      IF (PARCAL) THEN
C
C        Tell slaves to quit.
C        --------------------
C
         IPRTYP = 0
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(IPRTYP,1,'INTEGER')
         CALL PVMXMCAST(NODTOT,NODEID(1),MTAG1)
         CALL PVMXEXIT
      END IF
#endif
      STOP '*** End of DALTON calculation ***'
      END
