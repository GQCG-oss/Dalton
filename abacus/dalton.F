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
C===========================================================================
C old log:
C010131-vebjornb: Old dalton.F split into dalton.F (program) and
C                 dalgnr.F (subroutines).
C951009-kenneth : Changed name from PSISTAR to DALTON
C950621-vebjornb: First working optimization, GRAND changed to PSISTAR
C950516-vebjornb: First effort to link all programs: GRAND
C===========================================================================
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
      CHARACTER*14 FILNAM
      PARAMETER (LWORK_DEF = INSTALL_WRKMEM)
#if defined (SYS_CRAY) || defined (SYS_IRIX) || defined (SYS_T3D) || defined (SYS_SUN) || defined (SYS_HPUX) || defined (SYS_HAL) || defined (SYS_T90) || defined (SYS_DEC) || (defined (SYS_LINUX) && defined (VAR_PGF77))
      POINTER (WADDR,WORK(1))
#else
#if (defined (SYS_LINUX))
      DIMENSION WORK(2)
#else
#if defined (SYS_NEC) || defined (SYS_AIX)
#if defined (NO_F90)
      POINTER (WADDR,WORK(1))
#else
      DOUBLE PRECISION , ALLOCATABLE , DIMENSION(:) :: WORK
      INTEGER                           :: WADDR
#endif
#else
C
C     Static memory allocation
C
      DIMENSION WORK(LWORK_DEF+2)
#endif
#endif
#endif
      DATA WRKDLM/8H*WRKDLM*/
C
#if defined (VAR_MPI)
      CALL MPIXINIT
      PARIO = .FALSE.
#else
      MYNUM = 0
#endif
      CALL GPIOIN
      IF (MYNUM .EQ. 0) THEN
         CALL GPOPEN(LUERR,'DALTON.ERR','UNKNOWN',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
      ELSE
         LUERR = 0
      END IF
      LMWORK = 0
C
C     Read environment variables for LMWORK and BASDIR
C     If LMWORK .le. 0 use default value (LWORK_DEF)
C
      CALL GETMMBAS(LMWORK,BASDIR,IDUMMY)
      IF (LMWORK .LE. 0) LMWORK = LWORK_DEF
C
      NBYTES  = (LMWORK+2) * 8
      XMBYTES = (LMWORK+2) * 8.D0
      XMBYTES = XMBYTES / (1.024D3**2)
      GBYTES  = XMBYTES / 1.024D3
      IF (GBYTES .GE. 1.D0) THEN
        WRITE(LUERR,'(/A,I12,A,F7.3,A)')' Work memory size (LMWORK+2):',
     &        LMWORK+2,' =',GBYTES,' gigabytes.'
      ELSE
        WRITE(LUERR,'(/A,I12,A,F8.2,A)')' Work memory size (LMWORK+2):',
     &        LMWORK+2,' =',XMBYTES,' megabytes.'
      END IF
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

#if defined (SYS_NEC) || defined (SYS_AIX)
#if defined (NO_F90)
      WADDR = MALLOC(NBYTES)
#else
      WADDR = 1
      ALLOCATE (WORK(LMWORK+2),STAT=I)
      IF (I .NE. 0) THEN
         WRITE (LUERR,*) 'ALLOCATE command failed to allocate the '//
     &        'requested memory'
         WRITE (LUERR,*) 'Reduce the memory demands and be welcome back'
         CALL QUIT('Failed to allocate memory')
      END IF
#endif
#endif
#if (defined (SYS_LINUX))
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
#if !(defined (SYS_LINUX))
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
C        ... GPOPEN will change name from DALTON.OUT to $WRKDIR/$OUTFIL
C            if OUTFIL is an environment variable
         CALL GPOPEN(LUCME,'DALTON.CM','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      ELSE
         WRITE (FILNAM,'(A,I4.4)') 'DALTON.OUT',MYNUM
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
     * ' This is output from DALTON (Release 2.0 rev. 0, Feb. 2005)',
     *  ' Principal authors:',
     *  ' Trygve Helgaker,       University of Oslo,           Norway ',
     *  ' Hans Joergen Jensen,   SDU - Odense University,      Denmark',
     *  ' Poul Joergensen,       Aarhus University,            Denmark',
     *  ' Jeppe Olsen,           Aarhus University,            Denmark',
     *  ' Kenneth Ruud,          University of Tromsoe,        Norway ',
     *  ' Hans Aagren,           KTH Stockholm,                Sweden ',
     *  ' Contributors:',
     *  ' Celestino Angeli,      University of Ferrara,        Italy  ',
     *  ' Alexander A. Auer,     University of Mainz,          Germany',
     *  ' Keld L. Bak,           UNI-C,                        Denmark',
     *  ' Vebjoern Bakken,       University of Oslo,           Norway ',
     *  ' Ove Christiansen,      University of Lund,           Sweden ',
     *  ' Renzo Cimiraglia,      University of Ferrara,        Italy  ',
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
     *  ' Ola B. Lutnaes,        University of Oslo,           Norway ',
     *  ' Kurt V. Mikkelsen,     University of Copenhagen,     Denmark',
     *  ' Patrick Norman,        University of Linkoeping,     Sweden ',
     *  ' Thomas B. Pedersen,    U. of Santiago de Compostela, Spain  ',
     *  ' Martin J. Packer,      University of Sheffield,      UK     ',
     *  ' Zilvinas Rinkevicius   KTH Stockholm,                Sweden ',
     *  ' Torgeir A. Ruden,      University of Oslo,           Norway ',
     *  ' Pawel Salek,           KTH Stockholm,                Sweden ',
     *  ' Alfredo Sanchez,       University of Valencia,       Spain  ',
     *  ' Trond Saue,            University of Toulouse,       France ',
     *  ' Stephan P. A. Sauer,   University of Copenhagen,     Denmark',
     *  ' Bernd Schimmelpfennig, KTH Stockholm,                Sweden ',
     *  ' K. O. Sylvester-Hvid,  University of Copenhagen,     Denmark',
     *  ' Peter R. Taylor,       University of Warwick,        UK     ',
     *  ' Olav Vahtras,          KTH Stockholm,                Sweden ',
     *  ' David J. Wilson,       University of Oslo,           Norway '
C
C        Changed format (WK/UniKA/04-11-2002).
         WRITE (LUPRI,'(/1X,69A1/)') ('-',I=1,69)
         WRITE (LUPRI,'(5X,A)')
     &'NOTE:',
     &' ',
     &'This is an experimental code for the evaluation of molecular',
     &'properties using (MC)SCF and CC wave functions. The authors',
     &'accept no responsibility for the performance of the code or',
     &'for the correctness of the results.',
     &' ',
     &'The code (in whole or part) is provided under a licence and',
     &'is not to be reproduced for further distribution without',
     &'the written permission of the authors or their representatives.',
     &' ',
     &'See the home page "http://www.kjemi.uio.no/software/dalton"',
     &'for further information.',
     &' ',
     &'If results obtained with this code are published,',
     &'an appropriate citation would be:',
     &' ',
     &'"Dalton, a molecular electronic structure program, Release 2.0',
     &'(2005), see http://www.kjemi.uio.no/software/dalton/dalton.html"'
C
C     Stamp date and time and hostname to output
C
         CALL TSTAMP(' ',LUPRI)
C
C     Read general input
C
         CALL GNRLIN
         IF (TESTIN) THEN
            CALL QUIT('*** End of input test for GENERAL ***')
         END IF
#if !defined (VAR_MPI)
      END IF
#else /* VAR_MPI */
         CALL RECVNAMES
#if defined (SYS_LINUX)
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
#if defined (SYS_LINUX)
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
#if defined (SYS_LINUX)
            CALL CEXE(NODSTR,LMWORK,WRKDLM,MASTER,MYNUM)
#else
            CALL NODSTR(WORK,LMWORK,WRKDLM,MASTER,MYNUM)
#endif
         END IF
      END IF
      CALL MPIXFINALIZE
C
C     All below is not relevant for MPI as all proc. has stopped
C     with CALL MPI_FINALIZE(IERR)
C
#endif /* VAR_MPI */
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
#if defined (SYS_LINUX)
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
      END
      SUBROUTINE GETMMBAS(LMWORK,BASDIR,LBASDR)
#include <implicit.h>
#include <priunit.h>
      CHARACTER*70 BASDIR
#if defined(SYS_AIX)||defined(SYS_CRAY)||defined(SYS_IRIX)||defined(SYS_NEC)||defined(SYS_CONVEX)||defined(SYS_T3D)||defined(SYS_SUN)||defined (SYS_LINUX)||defined(SYS_HPUX)||defined(SYS_HAL)||defined(SYS_T90)||defined(SYS_DEC)
      CHARACTER*20 WRKMEM
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
C
C     Dynamic memory allocation; check for WRKMEM in environment
C     The T3D thingie might work for DEC_ALPHA as well? Try it out, kr-dec96
C
      WRKMEM = '                    '
#if defined (SYS_T3D) || defined (SYS_T90)
      CALL PXFGETENV('WRKMEM',6,WRKMEM,LWRKMM,IERR)
      CALL PXFGETENV ('BASDIR',6,BASDIR,LBASDR,IERR)
#else
      CALL GETENV('WRKMEM',WRKMEM)
      CALL GETENV ('BASDIR',BASDIR)
#endif
      READ(WRKMEM, '(I20)', ERR=10) LMWORK
      IF (LMWORK .GE. 0) GOTO 20
   10    WRITE(LUERR,'(/3A)')
     &   ' DALTON: WRKMEM conversion error; WRKMEM = "',WRKMEM,'"'
         LMWORK = 0
   20 CONTINUE
C
C     Either default memory size,
C
      IF (LMWORK .LE. 0) THEN
         WRITE (LUERR,'(/A)') ' DALTON: default work memory size used.'
      ELSE
C
C     ... or user specified memory size.
C
         DO I=20,1,-1
           IF(WRKMEM(I:I).NE.' ') GOTO 1312
         ENDDO
         CALL QUIT('Empty environment variable WRKMEM')
 1312    LWMEM = I
         WRITE(LUERR,'(/A/3A)')
     &   ' DALTON: user specified work memory size used,',
     &   '         environment variable WRKMEM = "',WRKMEM(1:LWMEM),'"'
      END IF
#if defined(VAR_MPI)
      CALL MPI_BCAST(LMWORK,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#endif /* defined(VAR_MPI) */
#endif
      END
