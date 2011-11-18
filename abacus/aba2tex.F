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
C  /* Deck twoinp */
      SUBROUTINE TWOINP(WORD)
C
C     TUH
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 20)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <abainf.h>
#include <cbitwo.h>
      DATA TABLE /'.SKIP  ', '.PRINT ', '.DIRTST', '.FIRST ', '.SECOND',
     *            '.PTRSKI', '.SORSKI', '.INTSKI',
     *            '.PTRPRI', '.SORPRI', '.INTPRI',
     *            '.NODC  ', '.NODV  ', '.NOPV  ', '.RETURN',
     *            'XXXXXXX', '.TIME  ', '.NOCONT', '.STOP  ','.PTRNOD'/
C
      CALL QENTER('TWOINP')
      NEWDEF = WORD .EQ. '*TWOEXP'
      MAXOLD = MAXDIF
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,
     &                      11,12,13,14,15,16,17,18,19,20),I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in TWOINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in TWOINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRALL
                  IF (IPRALL .EQ. IPRDEF) ICHANG = ICHANG - 1
                  IPRINT = IPRALL
                  IPRPTR = IPRALL
                  IPRSOR = IPRALL
               GO TO 100
    3             CONTINUE
                  DIRTST = .TRUE.
               GO TO 100
    4             MAXDIF = 1
                  IF (MAXDIF .EQ. MAXOLD) ICHANG = ICHANG - 1
               GO TO 100
    5             MAXDIF = 2
                  IF (MAXDIF .EQ. MAXOLD) ICHANG = ICHANG - 1
               GO TO 100
    6             RUNPTR = .FALSE.
               GO TO 100
    7             RUNSOR = .FALSE.
               GO TO 100
    8             RUNINT = .FALSE.
               GO TO 100
    9             READ (LUCMD,*) IPRPTR
                  IF (IPRPTR .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   10             READ (LUCMD,*) IPRSOR
                  IF (IPRSOR .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   11             READ (LUCMD,*) IPRINT, IPRNTA, IPRNTB,
     *                                          IPRNTC, IPRNTD
                  IPRSUM = IPRNTA + IPRNTB + IPRNTC + IPRNTD
                  IF (IPRINT .EQ. IPRDEF .AND. IPRSUM .EQ. 0) THEN
                     ICHANG = ICHANG - 1
                  END IF
               GO TO 100
   12             NODC = .TRUE.
               GO TO 100
   13             NODV = .TRUE.
               GO TO 100
   14             NOPV = .TRUE.
               GO TO 100
   15             RETUR = .TRUE.
               GO TO 100
   16             CONTINUE
C              GO TO 100
   17             TKTIME = .TRUE.
               GO TO 100
   18             NOCONT = .TRUE.
               GO TO 100
   19             CUT    = .TRUE.
               GO TO 100
   20             NODPTR = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in TWOINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in TWOINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .EQ. 0) RETURN
      IF (NEWDEF) THEN
         CALL HEADER('Changes of defaults for TWOEXP:',0)
      END IF
      IF (SKIP) THEN
         WRITE (LUPRI,'(A)') ' TWOEXP skipped in this run.'
      ELSE IF (NEWDEF) THEN
         IF (IPRALL .NE. IPRDEF) THEN
            WRITE (LUPRI,'(A,I5)') ' Print level in TWOEXP:',IPRALL
         END IF
         IF (RUNPTR) THEN
            IF (IPRPTR .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in PTRAN :',IPRPTR
            END IF
            IF (NODPTR) THEN
               WRITE (LUPRI,'(A)') ' PV matrix transformed to AO '
     &              //'basis using noddy routine for comparison.'
            END IF
         ELSE
            WRITE (LUPRI,'(A)') ' PTRAN skipped in this run.'
         END IF
         IF (RUNSOR) THEN
            IF (IPRSOR .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in PSORT :',IPRSOR
            END IF
         ELSE
            WRITE (LUPRI,'(A)') ' PSORT skipped in this run.'
         END IF
         IF (RUNINT) THEN
            IF (MAXDIF .NE. MAXOLD) THEN
               WRITE (LUPRI,'(A,I1)') ' Maximum differentiation: ',
     *                                MAXDIF
            END IF
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in TWOINT:',
     *                                IPRINT
            END IF
            IF (IPRNTA + IPRNTB + IPRNTC + IPRNTD .GT. 0) THEN
               WRITE(LUPRI,'(2A,4I3)')' Extra output for the following',
     *                    ' shells:', IPRNTA, IPRNTB, IPRNTC, IPRNTD
               IF (RETUR) WRITE (LUPRI,'(A)')
     *             ' Program will exit TWOINT after these shells.'
            END IF
            IF (NODC) WRITE (LUPRI,'(/,2A)') ' Inactive one-electron',
     *         ' density matrix neglected in TWOEXP.'
            IF (NODV) WRITE (LUPRI,'(/,2A)') ' Active one-electron',
     *         ' density matrix neglected in TWOEXP.'
            IF (NOPV) WRITE (LUPRI,'(/,2A)') ' Active two-electron',
     *         ' density matrix neglected in TWOEXP.'
            IF (NOCONT) WRITE (LUPRI,'(/,2A)')
     *         ' Derivative integrals will not be contracted.'
            IF (TKTIME) WRITE (LUPRI,'(/,2A)') ' Detailed timing for',
     *         ' integral calculation will be provided.'
            IF (DIRTST) WRITE (LUPRI,'(/A)') ' Direct calculation '//
     &         'of Fock matrices and integral distributions tested.'
         ELSE
            WRITE (LUPRI,'(A)') ' TWOINT skipped in this run.'
         END IF
         IF (CUT) THEN
            WRITE (LUPRI,'(/,A)') ' Program is stopped after TWOEXP.'
         END IF
      END IF
      CALL QEXIT('TWOINP')
      RETURN
      END
C  /* Deck twoini */
      SUBROUTINE TWOINI
C
C     Initialize /CBITWO/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbitwo.h>
C
      CALL QENTER('TWOINI')
      SKIP   = .FALSE.
      CUT    = .FALSE.
      RUNPTR = .TRUE.
      RUNSOR = .TRUE.
      RUNINT = .TRUE.
      IPRALL = IPRDEF
      IF (MOLHES .AND. .NOT. HELFEY) THEN
         MAXDIF = 2
      ELSE IF ((MOLGRD .OR. DIPDER .OR. QPGRAD) .AND. .NOT. HELFEY) THEN
         MAXDIF = 1
      ELSE IF (EXPGRD) THEN
         MAXDIF = 2
      ELSE
         SKIP = .TRUE.
      END IF
      IPRINT = IPRDEF
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      IPRPTR = IPRDEF
      IPRSOR = IPRDEF
      NODC   = .FALSE.
      NODV   = .FALSE.
      NOPV   = .FALSE.
      NOCONT = .FALSE.
      RETUR  = .FALSE.
      TKTIME = .FALSE.
      DIRTST = .FALSE.
      NODPTR = .FALSE.
      CALL QEXIT('TWOINI')
      RETURN
      END
C  /* Deck twoexp */
      SUBROUTINE TWOEXP(WORK,LWORK,PASS)
C
C     TUH
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <aovec.h>
      LOGICAL PASS, ANTI, PANTI, ANTIDM, DIA2SO, DOERI
      DIMENSION WORK(LWORK)
#include <abainf.h>
#include <exeinf.h>
#include <cbieri.h>
#include <cbitwo.h>
#include <cbisol.h>
#include <inforb.h>
#include <inftap.h>
#include <symmet.h>
#include <nuclei.h>
#include <gnrinf.h>
#include <infpar.h>
#include <infinp.h>
#include <taymol.h>
#include <taysol.h>
#include <energy.h>
#include <expopt.h>
#include <dftcom.h>
#include <rspprp.h>
#include <esg.h>
C
      PARAMETER (D0 = 0.0D0)
      PARAMETER (MXDMAT = 2*MXCOOR)
      LOGICAL DDFOCK, NOBLK, ZFS2EL, PARCAL_IS_OK, PARCAL_SAVE
      INTEGER IFCTYP(MXDMAT), ISYMDM(MXDMAT)
C
C     Control routine for calculation of expectation values
C     of differentiated two-electron integrals
C
      IF (SKIP) RETURN
      CALL QENTER('TWOEXP')
C
C     If not executed earlier, variables are saved, else restored.
C     Necessary for unified code.
C
      IF (FTWOXP) THEN
         NASTXP = NASHT
         SORTXP = RUNSOR
         PTRTXP = RUNPTR
         NPVTXP = NOPV
         FTWOXP = .FALSE.
      ELSE
         NASHT  = NASTXP
         RUNSOR = SORTXP
         RUNPTR = PTRTXP
         NOPV   = NPVTXP
      END IF
C
      I2TYP = 0
      IF (IPRALL .GE. 1) THEN
         CALL TIMER('START ',TIMEIN,TIMOUT)
         WRITE (LUPRI,'(/A/)')
     *     '  <<<<<<<<<< Output from TWOEXP >>>>>>>>>> '
      END IF
C
C     No active orbitals
C
      IF (NASHT .EQ. 0) NODV = .TRUE.
      IF (NASHT .LT. 2 .OR. HSROHF) THEN
         NOPV   = .TRUE.
         IF (HSROHF .AND. (DIPDER .OR. MOLHES)) THEN
            WRITE(LUPRI,'(/A)') ' Dipole gradients and molecular '//
     &           'not implmented for high-spin restricted open-shell HF'
            CALL QUIT('Option not available for HSROHF in TWOEXP')
         END IF 
      END IF
      IF (NOPV) THEN
         RUNPTR = .FALSE.
         RUNSOR = .FALSE.
      END IF
C
C     *******************************************************
C     ***** Set up COMMON /BLOCKS/ for PSORT and TWOINT *****
C     *******************************************************
C
      IF (RUNSOR .OR. RUNPTR) THEN
         NOBLK = .TRUE.
      ELSE
         NOBLK = .FALSE.
      END IF
      KJSTRS = 1
      KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KIORBS = KNCONT + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KJORBS = KIORBS + (MXSHEL*MXAOVC + 1)/IRAT
      KKORBS = KJORBS + (MXSHEL*MXAOVC + 1)/IRAT
      KLAST  = KKORBS + (MXSHEL*MXAOVC + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('TWOEXP','PAOVEC',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,NOBLK,IPRALL)
      KLAST = KJORBS
      LWRK = LWORK - KLAST + 1
C
C     **********************************************************
C     ***** Set up two-electron density matrix in AO basis *****
C     **********************************************************
C
      ANTI = .FALSE.
      PANTI = .FALSE.
      DIA2SO = .FALSE.
      ZFS2EL = .FALSE.
C
C     Transformation to AO basis
C
      IF (RUNPTR) THEN
         IF (IPRPTR .GE. 2) CALL TIMER('START ',TIMSTR,TIMEND)
         CALL PTRAN(NODPTR,WORK(KLAST),LWRK,IPRPTR,ANTI,PANTI,DIA2SO,
     &              ZFS2EL)
         IF (IPRPTR .GE. 2) CALL TIMER('PTRAN ',TIMSTR,TIMEND)
      END IF
C
C     Final sorting of integrals
C
      IF (RUNSOR) THEN
         IF (IPRPTR .GE. 2) CALL TIMER('START ',TIMSTR,TIMEND)
C
C     In order to avoid paging on machines with virtual memory, we
C     use at most 64 Mb on such machines as suggested by P.Taylor
C     K.Ruud, Aug.-98
C
#if defined (SYS_CRAY) || defined (SYS_NEC) || defined (SYS_T90)
         LWUSE = LWRK
#else
         CALL MEMSRT(MEMMIN,WORK(KNCONT))
C
C        1000000 is just an estimate of what the rest of the sorting 
C        program needs
C
         LWUSE = MEMMIN + 1000000
         IF (LWUSE .GT. LWRK)
     &        CALL STOPIT('TWOEXP','PSORG',LWUSE,LWRK)
#endif
         IF (LUPAO .LE. 0)
     &        CALL GPOPEN(LUPAO,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
         CALL PSORG(WORK(KLAST),WORK(KLAST),LWUSE,WORK(KNCONT),IPRSOR,
     &              ANTI,PANTI)
         IF (IPRPTR .GE. 2) CALL TIMER('PSORT ',TIMSTR,TIMEND)
      END IF
C
C     ****************************************
C     ***** Calculate Expectation Values *****
C     ****************************************
C
      IF (RUNINT) THEN
C
C        Zero EE gradient and Hessian, if calculated.
C        Allocate memory for HESSEE.
         IF (MOLGRD) CALL DZERO(GRADEE,MXCOOR)
         KHESEE = KLAST
         IF (MOLHES .OR. PARCAL) THEN
           KLAST  = KHESEE + MXCOOR*MXCOOR
           CALL DZERO(WORK(KHESEE),MXCOOR*MXCOOR)
         ENDIF
C
C        One-electron density matrices
C
         IF (NODV) THEN
            NDMAT = 1
         ELSE
            NDMAT = 2
         END IF
C
C  4 generalised density matrices are needed for excited state gradient
C  calculations

         IF (ESG) NDMAT=4
C
         KDMAT = KLAST
         IF (FCKDDR .AND. EXPFCK) THEN
            ITYPE = 6
CTROND
            DDFOCK = MAXDIF.GT.1
CTROND
CKR
CKR Think carefully through whether next test will ever become true
CKR
            IF ((DIPDER .OR. QPGRAD) .AND. .NOT. MOLHES) THEN
               ITYPE = -6
               DDFOCK = .TRUE.
            END IF
            NFMAT = 3*NUCDEP*NDMAT
            KFMAT = KDMAT + NDMAT*N2BASX
            KLAST = KFMAT + NFMAT*N2BASX
         ELSE
            ITYPE = 2
            DDFOCK = .FALSE.
            NFMAT = 0
            KLAST = KDMAT + NDMAT*N2BASX
            KFMAT = KDMAT
         END IF
C
         LWRK  = LWORK - KLAST + 1
         IF(KLAST.GT.LWORK) CALL STOPIT('TWOEXP','GETDMT',KLAST,LWORK)
         IF (NFMAT .GT. 0) THEN
            CALL DZERO(WORK(KFMAT),NFMAT*N2BASX)
C           dimension ISYMDM(MXDMAT), IFCTYP(MXDMAT)
            IF (NFMAT .GT. MXDMAT) THEN
               WRITE (LUPRI,'(/A,2I10)')
     &         ' MXDMAT too small in TWOEXP; NFMAT,MXDMAT=',NFMAT,MXDMAT
               CALL QUIT('MXDMAT too small in TWOEXP')
            END IF
            DO I = 1,NFMAT
               ISYMDM(I) = 0
               IFCTYP(I) = 13
            END DO
         END IF
C
C        Is it OK to do parallel calculation here ??
         PARCAL_IS_OK = PARCAL ! has user wanted PARCAL ?
         IF ((NASHT .GT. 0) .AND. (.NOT.(HSROHF.OR.DFTADD)))
     &      PARCAL_IS_OK = .FALSE. ! parallel TWOEXP only for SCF 
C
C        ERI can only be used for Hartree-Fock gradients
C
         DOERI = RUNERI .AND. .NOT. PARCAL_IS_OK .AND. NODV .AND. NOPV
     &                  .AND. ITYPE.EQ.2 .AND. MAXDIF .LT. 2
C
         CALL GETDMT(WORK(KDMAT),NDMAT,WORK(KLAST),LWRK,NODC,NODV,
     &               .NOT.DOERI,IPRINT)
C
C  For excited state gradient calculation we need to read 
C    3 more generalized density matrices from file
C
      IF (ESG) CALL GETESG_DENMAT(WORK(KDMAT),NDMAT,WORK(KLAST),
     &               LWRK)
C
#if defined (VAR_MPI) || defined (VAR_PVM)
         IF (PARCAL_IS_OK) THEN
            IF ((NASHT .GT. 0) .AND. (.NOT.(HSROHF.OR.DFTADD)))
     &         CALL QUIT(' ERROR in TWOEXP: PARCAL only'// 
     &                   ' implemented for SCF')
C
            KNSTAT = KLAST
            KLAST  = KNSTAT + (NODTOT + 1)/IRAT
            IF (KLAST .GT. LWRK)
     &         CALL STOPIT('TWOEXP','PARDRV',KLAST,LWRK)
            LWRK = LWRK - KLAST
C
            IATOM = 0
            IPRTYP = 1
            HFXM0 = HFXMU
            HFXMU = D0
            CALL PARDRV(WORK(KFMAT),WORK(KDMAT),NDMAT,ISYMDM,IFCTYP,
     &                  WORK(KLAST),WORK(KNSTAT),WORK(KHESEE),LWRK,
     &                  ITYPE,MAXDIF,IATOM,NODV,NOPV,NOCONT,TKTIME,
     &                  RETUR,IPRINT,IPRTYP,
     &                  IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,
     &                  DUMMY)
            HFXMU = HFXM0
            IF (HFXMU.NE.D0) THEN

               KFCTMP = KLAST
               KLAST  = KFCTMP + NFMAT*N2BASX
               IF (KLAST .GT. LWRK)
     &              CALL STOPIT('TWOEXP','PARDRV',KLAST,LWRK)
               LWRK = LWRK - NFMAT*N2BASX

               HFXFC0 = HFXFAC
               HFXFAC = HFXATT
               IATOM = 1
               IF (NFMAT .GT. 0) CALL DCOPY(NFMAT*N2BASX,WORK(KFMAT),1,
     &                                      WORK(KFCTMP),1)
               CALL PARDRV(WORK(KFMAT),WORK(KDMAT),NDMAT,ISYMDM,IFCTYP,
     &              WORK(KLAST),WORK(KNSTAT),WORK(KHESEE),LWRK,
     &              ITYPE,MAXDIF,IATOM,NODV,NOPV,NOCONT,TKTIME,
     &              RETUR,IPRINT,IPRTYP,
     &              IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,
     &              DUMMY)
               IF (NFMAT .GT. 0) CALL DAXPY(NFMAT*N2BASX,1D0,
     &                                     WORK(KFCTMP),1,WORK(KFMAT),1)
               HFXFAC = HFXFC0
            END IF
         ELSE
            PARCAL_SAVE = PARCAL
            PARCAL = .FALSE.
#endif
            IF (DOERI) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               KCCFBT = KLAST
               KINDXB = KCCFBT + MXPRIM*MXCONT
               KLAST  = KINDXB + 8*MXSHEL*MXCONT
               IF (KLAST.GT.LWRK) 
     &            CALL STOPIT('TWOEXP','ERIPRO',KLAST,LWRK)
               LWRK = LWRK - KLAST
               CALL ERIPRO(WORK(KFMAT),WORK(KDMAT),NDMAT,ISYMDM,
     &                     IFCTYP,.TRUE.,.FALSE.,IPRINT,WORK(KCCFBT),
     &                     WORK(KINDXB),WORK(KLAST),LWRK)
               CALL TIMER('ERIPRO',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            ELSE
               IF (IPRINT .GE. 2) CALL TIMER('START ',TIMSTR,TIMEND)
               HFXM0 = HFXMU
               HFXMU = D0
               CALL TWOINT(WORK(KLAST),LWRK,WORK(KHESEE),WORK(KFMAT),
     &                     WORK(KDMAT),NDMAT,ISYMDM,IFCTYP,DUMMY,IDUMMY,
     &                     NUMDIS,1,ITYPE,MAXDIF,0,NODV,NOPV,NOCONT,
     &                     TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,
     &                     RETUR,IDUMMY,I2TYP,WORK(KJSTRS),WORK(KNPRIM),
     &                     WORK(KNCONT),WORK(KIORBS),IDUMMY,IDUMMY,
     &                     DUMMY,DUMMY,DUMMY,DUMMY,.FALSE.)
               HFXMU = HFXM0
               IF (HFXMU.NE.D0) THEN
                  HFXFC0 = HFXFAC
                  HFXFAC = HFXATT
                  CALL TWOINT(WORK(KLAST),LWRK,WORK(KHESEE),WORK(KFMAT),
     &                     WORK(KDMAT),NDMAT,ISYMDM,IFCTYP,DUMMY,IDUMMY,
     &                     NUMDIS,1,ITYPE,MAXDIF,0,NODV,NOPV,NOCONT,
     &                     TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,
     &                     RETUR,IDUMMY,I2TYP,WORK(KJSTRS),WORK(KNPRIM),
     &                     WORK(KNCONT),WORK(KIORBS),IDUMMY,IDUMMY,
     &                     DUMMY,DUMMY,DUMMY,DUMMY,.FALSE.)
                  HFXFAC = HFXFC0
               END IF
               IF (IPRINT .GE. 2)
     &         CALL TIMER('TWOEXP/TWOINT',TIMSTR,TIMEND)
            END IF
#if defined (VAR_MPI) || defined (VAR_PVM)
            PARCAL = PARCAL_SAVE
         END IF
#endif
CTROND
         IF (.NOT.DOERI .AND..NOT.EXPGRA) THEN
            CALL SKLFCK(WORK(KFMAT),WORK(KHESEE),WORK(KLAST),LWRK,
     &                  IPRINT,.FALSE.,
     &                  DDFOCK,.TRUE.,.FALSE.,NODV,MAXDIF,.FALSE.,NDMAT,
     &                  ISYMDM,IFCTYP,0,.FALSE.)
         END IF
C
C     <<<<< Print Section - Gradient and Hessian Elements >>>>>
C
         IF (MOLHES) CALL ADDHES(WORK(KHESEE))
         KCSTRA = 1
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
         IF (KLAST .GT. LWORK) 
     &        CALL STOPIT('TWOLOP','PRINT ',KLAST,LWORK)
         IF (IPRINT .GT. 0) THEN
            IF (EXPGRA) THEN
               CALL HEADER('Alpha gradient',-1)
               WRITE (LUPRI,'(2X,5F12.8)') (ALPGRD(I),I=1,NPBAS)
            ELSE
               CALL HEADER('Two-electron integral gradient',-1)
               CALL PRIGRD(GRADEE,WORK(KCSTRA),WORK(KSCTRA))
               CALL HEADER
     &              ('Potential energy (NN + NE + EE) gradient',-1)
               CALL ZERGRD
               CALL ADDGRD(GRADNN)
               CALL ADDGRD(GRADNA)
               CALL ADDGRD(GRADEE)
               CALL PRIGRD(GRDMOL,WORK(KCSTRA),WORK(KSCTRA))
               CALL HEADER('Molecular gradient',-1)
               CALL ADDGRD(GRADFS)
               CALL ADDGRD(GRADKE)
               IF (SOLVNT) THEN
                  CALL ADDGRD(GSOLTT)
                  CALL ADDGRD(GSOLNN)
               END IF
               CALL PRIGRD(GRDMOL,WORK(KCSTRA),WORK(KSCTRA))
               NCDEP3 = 3*NUCDEP
               GRDNRM = DDOT(NCDEP3,GRDMOL,1,GRDMOL,1)
               GRDNRM = SQRT(GRDNRM)
               WRITE (LUPRI,'(/19X,A,1P,E10.2)')
     *            'Molecular gradient norm:', GRDNRM
               CALL ZERGRD
               IF (MOLHES) THEN
                  CALL HEADER('Two-electron integral Hessian',-1)
                  CALL PRIHES(WORK(KHESEE),'CENTERS',WORK(KCSTRA),
     &                        WORK(KSCTRA))
               END IF
            END IF
         END IF
CTROND
      END IF
C
      IF (IPRALL .GE. 1) CALL TIMER ('TWOEXP',TIMEIN,TIMOUT)
      IF (RUNSOR) CALL GPCLOSE(LUPAO,'DELETE')
      PASS = .TRUE.
C
C     Testing of direct integral calculation
C
      IF (DIRTST) THEN
         CALL FCKTES(WORK(KLAST),LWRK,MAXDIF,NODV,NOPV,NOCONT,TKTIME,
     &               IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR)
      END IF
      IF (CUT) THEN
         WRITE (LUPRI,'(/A)')
     &          ' Program stopped after TWOEXP as required.'
         WRITE (LUPRI,'(A)') ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (in TWOEXP) *****')
      END IF
      CALL QEXIT('TWOEXP')
      RETURN
      END
C  /* Deck getdmt */
      SUBROUTINE GETDMT(DMAT,NDMAT,WORK,LWORK,NODC,NODV,SYMUNQ,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>

#include <abainf.h>
#include <inforb.h>

      LOGICAL NODC, NODV, SYMUNQ
      DIMENSION DMAT(NBAST,NBAST,NDMAT), WORK(LWORK)

      CALL QENTER('GETDMT')
      IF (IPRINT .GT. 5) CALL TITLER('Output from GETDMT','*',103)
C
      KCMO  = 1
      KDV   = KCMO  + NCMOT
      KDTSQ = KDV   + NNASHX
      KDASQ = KDTSQ + N2BASX
      KLAST = KDASQ + N2BASX
      IF (KLAST.GT.LWORK) CALL STOPIT('GETDMT','ONEDSF',KLAST,LWORK)
C
      IF (SYMUNQ) THEN
         CALL ONEDSF(WORK(KCMO),WORK(KDV),WORK(KDTSQ),WORK(KDASQ),
     &               IPRINT,NODC,NODV)
         ISYMDM = 0
         CALL DSOTAO(WORK(KDTSQ),DMAT(1,1,1),NBAST,ISYMDM,IPRINT)
         IF (.NOT.NODV) CALL DSOTAO(WORK(KDASQ),DMAT(1,1,2),NBAST,
     &                              ISYMDM,IPRINT)
      ELSE
         CALL ONEDSF(WORK(KCMO),WORK(KDV),DMAT,WORK(KDASQ),
     &               IPRINT,NODC,NODV)
      END IF
C
      CALL QEXIT('GETDMT')
      RETURN
      END
C  /* Deck onedsf */
      SUBROUTINE ONEDSF(CMO,DV,DTSO,DASO,IPRINT,NODC,NODV)
C
C     This subroutine calculates the total and active one-electron
C     density matrices in SO basis (contravariant).  Input is
C     one-electron active density matrix.
C                                         880420  PRT
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <iratdef.h>
      PARAMETER (D0 = 0.0D0)
      INTEGER R, S, U, V, UV, UR, US, VR, VS
      LOGICAL NODC, NODV, FOUND
      DIMENSION CMO(*), DV(*), DTSO(NBAST,NBAST), DASO(NBAST,NBAST)
C
#include <abainf.h>
#include <inforb.h>
#include <inftap.h>
C
C     Read CMO and DV from SIRIFC
C
      CALL QENTER('ONEDSF')
      CALL RD_SIRIFC('CMO',FOUND,CMO,DTSO,N2BASX)
      IF (.NOT.FOUND) CALL QUIT('ONEDSF ERROR: CMO not found on SIRIFC')
      IF (NASHT .GT. 0) THEN
         CALL RD_SIRIFC('DV',FOUND,DV,DTSO,N2BASX)
         IF (.NOT.FOUND)
     &      CALL QUIT('ONEDSF ERROR: DV not found on SIRIFC')
      END IF
C
C     Print Section
C
      IF (IPRINT .GT. 03) THEN
         WRITE (LUPRI, '(//A/)') ' <<<<< Subroutine ONEDSF >>>>>'
         WRITE (LUPRI, '(A,8I5)') ' NISH ', (NISH(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NASH ', (NASH(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NOCC ', (NOCC(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NORB ', (NORB(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NBAS ', (NBAS(I),I = 1,NSYM)
         IF (IPRINT .GE. 05) THEN
            CALL HEADER('Occupied molecular orbitals',0)
            IEND = 0
            DO 1000 ISYM = 1,NSYM
               IF (NBAS(ISYM) .EQ. 0) GOTO 1000
               IF (NOCC(ISYM) .EQ. 0) GOTO 1100
               WRITE (LUPRI, '(//,A,I5,/)') ' Symmetry ', ISYM
               IENDI = 0
               DO 1200 I = 1, NOCC(ISYM)
                  WRITE (LUPRI, '(/,A,I5,/)')
     *                     ' Molecular orbital ', I
                  WRITE (LUPRI, '(6F12.6)')
     *               (CMO(IEND+IENDI+J), J = 1, NBAS(ISYM))
                  IENDI = IENDI + NBAS(ISYM)
1200           CONTINUE
1100           CONTINUE
               IEND = IEND + NORB(ISYM)*NBAS(ISYM)
1000        CONTINUE
            CALL HEADER('Active density matrix (MO basis)',-1)
            CALL OUTPAK(DV,NASHT,1,LUPRI)
         END IF
      END IF
C
C     ***** Construct contravariant SO matrices *****
C
      CALL DZERO(DTSO,N2BASX)
      CALL DZERO(DASO,N2BASX)
      ICEND = 0
      DO 110 ISYM = 1,NSYM
         DO 100 R = 1, NBAS(ISYM)
         DO 100 S = 1,R
C
            DTRS = D0
C
C           (I) Inactive contribution
C
            ICENDI = 0
            DO 300 I = 1, NISH(ISYM)
               DTRS = DTRS + CMO(ICEND+ICENDI+R)*CMO(ICEND+ICENDI+S)
               ICENDI = ICENDI + NBAS(ISYM)
  300       CONTINUE
            DTRS = DTRS + DTRS
            IF (NODC) DTRS = D0
C
C           (II) Active contribution
C
            DVRS = D0
            IF (.NOT. NODV) THEN
               IASHI = IASH(ISYM)
               UV = ((IASHI + 1)*(IASHI + 2))/2
               IDVEND = ICEND + NISH(ISYM)*NBAS(ISYM)
               ICENDU = IDVEND
               DO 400 U = 1,NASH(ISYM)
                  ICENDV = IDVEND
                  DO 410 V = 1, U
                     DUV = DV(UV)
                     IF (ABS(DUV) .GT. D0) THEN
                        TEMP = CMO(ICENDU+R)*CMO(ICENDV+S)
                        IF (U .NE. V) TEMP = TEMP
     *                       + CMO(ICENDU+S)*CMO(ICENDV+R)
                        DVRS = DVRS + DUV*TEMP
                     END IF
                     UV = UV + 1
                     ICENDV = ICENDV + NBAS(ISYM)
  410             CONTINUE
                  UV = UV + IASHI
                  ICENDU = ICENDU + NBAS(ISYM)
  400          CONTINUE
            END IF
            IR = IBAS(ISYM) + R
            IS = IBAS(ISYM) + S
            DTSO(IR,IS) = DTRS + DVRS
            DTSO(IS,IR) = DTRS + DVRS
            DASO(IR,IS) = DVRS
            DASO(IS,IR) = DVRS
  100    CONTINUE
         ICEND = ICEND + NORB(ISYM)*NBAS(ISYM)
110   CONTINUE
C
C     ***** Print Section *****
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('Total SO density matrix (not folded)',-1)
         CALL OUTPUT(DTSO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('Active SO density matrix (not folded)',-1)
         CALL OUTPUT(DASO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      CALL QEXIT('ONEDSF')
      RETURN
      END
! --- end of aba2tex.F ---
