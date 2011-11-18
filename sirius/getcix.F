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
#ifdef REV_LOG
c ===========================================================================
c Revision 1.2  2000/05/24 12:20:33  hjj
c 1) NCDETS new parameter in SETCI
c 2) GETDETS new subroutine for getting det. vector from csf vector
c ===========================================================================
#endif
C  /* Deck getcix */
      SUBROUTINE GETCIX(INDXCI,ICSYM,IHCSYM,WRK,LFRSAV,NOSYM)
C
C  20-Jan-1988 hjaaj+jo
C
C     GETCIX: get string information for determinant CI
C
C     NOTE: for density matrix routines
C         ILSYM = ICSYM  and KLOFF = KCOFF, KLDTAS = KCDTAS,
C                            KLCOOS = KICOOS, KLCSO = KICSO
C         IRSYM = IHCSYM and KROFF = KHCOFF, KRDTAS = KHDTAS,
C                            KRCOOS = KIHOOS, KRCSO = KIHCSO
C         This convention makes it possible to use same GETCIX call
C         for non-symmetric CISIGD and DENSID calls.  This convention
C         is the opposite of the one originally chosen by Jeppe Olsen
C         /890908-hjaaj
C
C Input:
C  none
C
C Output:
C  in INDXCI(*) : string information
C
#include <implicit.h>
      DIMENSION INDXCI(*), WRK(LFRSAV)
C
C Used from common blocks:
C   INFINP : NACTEL, ISPIN, FLAG(26), FLAG(27), LSYM
C   INFORB : MULD2H(8,8), NASHT,  NAS1(8), NAS2(8), NAS3(8)
C   INFPRI : IPRCIX
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
C
      LOGICAL     REFSPC
      CHARACTER*8 LAB123(3),LABSTR
      DATA LAB123 /'********','********','********'/
      DATA LABSTR /'CISTRING'/
      DATA IFIRST / 1 /

C     
      CALL QENTER('GETCIX')
      IF (IPRCIX .GT. 0) THEN
         WRITE (LUPRI,'(//A/)') ' >>>>> OUTPUT FROM GETCIX <<<<<'
         WRITE (LUPRI,'(A,3I5/)')
     *      ' ICSYM, IHCSYM, NOSYM :',ICSYM,IHCSYM,NOSYM
         DTIM = SECOND()
      END IF
      LFREE  = LFRSAV
      IF (MCTYPE .GT. 0) THEN
        LUSTR = -9999
        IF (FLAG(26) .AND. IFIRST .EQ. 1) THEN
          CALL GPOPEN(LUSTR,'SIRIUS.STRINGINFO','UNKNOWN',
     &                       ' ','UNFORMATTED',IDUMMY,.FALSE.)
          REWIND(LUSTR)
          CALL GETDAT(LAB123(2),LAB123(3))
          WRITE(LUSTR) LAB123,LABSTR
C         IFIRST = 0
C         ... write to LUSTR each call, CI space may change during
C             calculation in future. /hjaaj
        END IF
C
C980820-hjaaj: CSF's are only implented for the ref.sym. LSYM
C       in this version:
C
        IF (FLAG(27)) THEN
C        ... do not use CSF expansion
           ICSF = 0
        ELSE IF (ICSYM  .EQ. LSYM) THEN
           ICSF = -1
        ELSE IF (IHCSYM .EQ. LSYM) THEN
           ICSF = -2
        ELSE
C          none of ICSYM and IHCSYM are of reference symmetry
C          so no CSF's for this call.
           ICSF = 0
        END IF
        REFSPC = FLAG(28)
        CALL DETINF(NACTEL,NASHT,MULD2H,ISPIN,NCDET,NHCDET,ICSYM,IHCSYM,
     *              NAS1,NAS2,NAS3,IPRCIX,INDXCI,WRK,LFREE,LUSTR,
     *              NOSYM,ICSF,REFSPC)
C       CALL DETINF(NEL,NORB,SYMPRO,MULTS,NCDET,NHCDET,ICSYM,IHCSYM,
C    *              NRAS1,NRAS2,NRAS3,NTEST,XNDXCI,WORK,LFREE,LUSTR,
C    *              NOSYM,ICSF,REFSPC)
C
        IF (LUSTR.GT.0) CALL GPCLOSE(LUSTR,'KEEP')
      ELSE
C     ... RHF or ROHF
        IF (IPRCIX .GT. 0) WRITE (LUPRI,'(/A/A)')
     &     ' RHF with at most one open shell or high spin.',
     &     ' No determinant string information needed.'
      END IF
      IF (IPRCIX .GT. 0) THEN
         DTIM = SECOND() - DTIM
         WRITE (LUPRI,'(/A,F10.2,A)') ' Time used in GETCIX :',DTIM,'s'
      END IF
      CALL QEXIT('GETCIX')
C
      RETURN
      END
C  /* Deck getfij */
      SUBROUTINE GETFIJ(FIJ,FJI,H2AC)
C
C 20-Jan-1988 hjaaj :
C
C   Version for Determinant CI
C
C   extract FIJ and FJI (Coulumb and exchange integrals)
C   for CI diagonal routine
C
#include <implicit.h>
      DIMENSION FIJ(NASHT,NASHT),FJI(NASHT,NASHT),H2AC(NNASHX,NNASHX)
#include <iratdef.h>
C
C
C Used from common blocks:
C   INFORB : NASHT,NNASHX
C   INFIND : IROW(*)
C   INFTAP : LUH2AC
C   CBGETDIS : IADH2
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <inftap.h>
#include <cbgetdis.h>
C
C
C ***
C Construct FIJ(ij)  =  (ii/jj)
C           FJI(ij)  =  (ij/ij)
C
         DO 880 I = 1,NASHT
            IROWI = IROW(I)
            II    = IROW(I+1)
            DO 870 J = 1,I
               JJ    = IROW(J+1)
               IJ    = IROWI + J
               IF (IADH2 .GE. 0) THEN
                  CALL READDX(LUH2AC,IADH2+IJ,IRAT*NNASHX,H2AC)
                  IJX = 1
               ELSE
                  IJX = IJ
               END IF
               IF (I .EQ. J) THEN
                  DO 860 K = 1,NASHT
                     FIJ(K,I) = H2AC(IROW(K+1),IJX)
  860             CONTINUE
               END IF
               FJI(J,I) = H2AC(IJ,IJX)
               FJI(I,J) = FJI(J,I)
  870       CONTINUE
  880    CONTINUE
C
C     *** end of getfij ***
C
      RETURN
      END
C  /* Deck setci */
      SUBROUTINE SETCI(NCONF,NCDETS,ICSYM,WORK,LFREE,NOSYM)
C
C 8-Jan-1988 Hans Joergen Aa. Jensen
C 980820-hjaaj: extracted NEWCI check in separate function
C
C Set variables with CI information for
C reference wave fuction.
C
C for westa compatibility NOSYM has been added . if
C NOSYM is nonvanishing symmetry is neglected
C
#include <implicit.h>
      DIMENSION WORK(LFREE)
C
#include <iratdef.h>
C
C Used from common blocks:
C   INFINP : NACTEL,ISPIN, NELMN1, NELMX1, NELMN3, NELMX3
C   INFORB : NASHT,NAS1(*),NAS2(*),NAS3(*)
C   INFDIM : NCONMA,LCINDX,LACIMX,LBCIMX
C   INFPRI : IPRCIX
C   SCBRHF : IOPRHF
C   CIINFO : ICOMBI,IPSIGN, MXNDT, MXNCS
C   CBESPN : ISPIN1,ISPIN2
C   SPINFO : MULTS,MS2
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <scbrhf.h>
#include <ciinfo.h>
#include <spinfo.h>
#include <cbespn.h>
C
C
      LOGICAL  REFSPC, NEWCI, NOCSF, CHCKCI
      SAVE     REFSPC, NOCSF,
     &         LACISV, LBCISV,LCINSV
C
      CALL QENTER('SETCI ')
C
      IF (IPRCIX .GT. 0)
     &   WRITE (LUPRI,'(//A/)') ' >>>>> OUTPUT FROM SETCI <<<<<'
C
C CHCKCI returns true if new ci, either because it is first call
C or because specifications have changed.
C CHCKCI sets information in CIINFO and SPINFO
C Note that Hartree-Fock (i.e. nasht .le. 1) will not change
C CI information (NEWCI forced false).
C Note that change of symmetry will not change CI information.
C
      NEWCI = CHCKCI(NOSYM,IPRCIX,LUPRI)
      IF (NASHT .LE. 1 .OR. HSROHF) THEN
         NEWCI = .FALSE.
      ELSE
         NOCSF = FLAG(27)
         IF (ICSYM .NE. LSYM .AND. .NOT.NOCSF) THEN
C        980820-hjaaj: csf information is only ok for ref.sym. in
C        this version ...
            WRITE (LUPRI,'(//A/A,2I5)') ' Error in SETCI:'//
     &         ' SETCI called with ICSYM different from ref.sym. LSYM',
     &         ' ICSYM, LSYM:',ICSYM,LSYM
            CALL QUIT('SETCI called with ICSYM different from ref.sym.')
         END IF
C
         IF ( NOCSF ) THEN
C           ... do not use CSF expansion
            ICSF = 0
         ELSE
            ICSF = 1
         END IF
      END IF
C
C
      IF ( NEWCI ) THEN
C        890129/hjaaj: in this version only ICOMBI=0 and IPSIGN=1.
C        ICOMBI and IPSIGN are in CIINFO
         ICOMBI = 0
         IPSIGN = 1
C
         NCNSM  = 1
         REFSPC = FLAG(28)
         IF (REFSPC) THEN
            WRITE (LUPRI,'(//A)') ' Error in SETCI:'//
     &         ' REFSPC [= flag(28)] is not fully implemented yet'
            CALL QUIT('SETCI error: REFSPC is not implemented yet.')
CMAERKE-hjaaj-910730: code missing in NUMST2, NUMSTR, DETGN4, ?
         END IF
         LFREEL = LFREE
         IF (IPRCIX .GT. 0) THEN
            WRITE (LUPRI,'(A,I5)')
     *         ' Print level        :',IPRCIX,
     *         ' Reference symmetry :',ICSYM,
     *         ' CSF:1 or determ:0  :',ICSF,
     *         ' CSF spin multipl.  :',ISPIN,
     *         ' 2 * M_S value      :',MS2
            DTIM = SECOND()
         END IF
         CALL DETFO(NACTEL,NASHT,ISPIN,ICSYM,IPRCIX,
     *              NELMN1,NELMX3,NAS1,NAS2,NAS3,
     *              NCDET,NCCSF,LCINSV,LACISV,LBCISV,WORK,LFREEL,
     *              NELMX1,NELMN3,NOSYM,ICSF,NCNSM,REFSPC)
C        CALL DETFO(NEL,NORB,MULTS,ICSYM,NTEST,
C    *              NELMN1,NELMX3,NRAS1,NRAS2,NRAS3,
C    *              NCDET,NCCSF,LCINDX,LLOCA,LLOCB,WORK,LFREE,
C    *              MELMN3,NOSYM,ICSF,NCNSM,REFSPC)
C
C        LACISV and LBCISV contain information about core
C        allocation for lintrn/cilin routines (will be used
C        as LACIMX + NCSIM*LBCIMX).
C        LCINSV is length needed for core allocation to string
C        information.
C
C
C        Specify singlet-singlet 2-electron density matrix
C
         ISPIN1 = 0
         ISPIN2 = 0
C
         IF (IPRCIX .GT. 0) THEN
            WRITE (LUPRI,'(/A,3(/A,I10))')
     *      ' SETCI, core memory needed for CI:',
     *      ' LCINDX =',LCINSV,' LACIMX =',LACISV,' LBCIMX =',LBCISV
            WRITE (LUPRI,'(/A,I10)')
     *      ' Number of determinants:  ',NCDET,
     *      ' Number of configurations:',NCCSF
            DTIM = SECOND() - DTIM
            WRITE (LUPRI,'(/A,F10.2,A)')
     *         ' Time used in SETCI :',DTIM,'s'
         END IF
      END IF
C
      IF (NASHT .GT. 1 .AND. .NOT. HSROHF) THEN
         IF (ICSF .EQ. 0) THEN
C           ... do not use CSF expansion
            NCONF  = NDTASM(ICSYM)
            NCDETS = NCONF
            NCONMA = MXNDT
         ELSE
            NCONF  = NCSASM(ICSYM)
            NCDETS = NDTASM(ICSYM)
C           NCONMA = MXNCS
            NCONMA = MXNDT
C           ... we must have sufficient space
C               for determinant vectors of different
C               symmetries or triplet in ABACUS and RESPONSE
         END IF
         LACIMX = LACISV
         LBCIMX = LBCISV
         LCINDX = LCINSV
      ELSE
C     ... NASHT .le. 1 .or. HSROHF (Hartree-Fock)
         IF (HSROHF .OR. DOMC) THEN
C        ... high spin ROHF
            JCSYM = LSYM
         ELSE IF (NASHT .EQ. 1) THEN
C        ... one open shell ROHF
            JCSYM = IOPRHF
         ELSE 
            JCSYM = 1
         END IF
         IF (ICSYM .EQ. JCSYM) THEN
            NCONF = 1
            NCDETS = 1
         ELSE
            NCONF = 0
            NCDETS = 0
         END IF
         NCONMA = 1
         LACIMX = 0
         LBCIMX = 0
         LCINDX = 0
      END IF
C
      CALL QEXIT('SETCI ')
      RETURN
      END
C  /* Deck setci2 */
      SUBROUTINE SETCI2(NCONF,ICSYM,TRIPLET,NOSYM)
C
C 20-Aug-1998 Hans Joergen Aa. Jensen, based on SETCI
C
C Set variables with CI information for
C response equations (called from ABAVAR and RSPVAR)
C The cloning of SETCI2 from SETCI makes it possible
C to use CSF's in the reference symmetry while using
C determinants in other symmetries.
C
C for westa compatibility NOSYM has been added . if
C NOSYM is nonvanishing symmetry is neglected
C
#include <implicit.h>
C
      LOGICAL TRIPLET
C
#include <iratdef.h>
C
C Used from common blocks:
C   INFINP : NACTEL,ISPIN, NELMN1, NELMX1, NELMN3, NELMX3
C   INFORB : NASHT,NAS1(*),NAS2(*),NAS3(*)
C   INFDIM : NCONMA,LCINDX,LACIMX,LBCIMX
C   SCBRHF : IOPRHF
C   CIINFO : MXNDT, MXNCS
C   CBESPN : ISPIN1,ISPIN2
C   SPINFO : MULTS,MS2
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infdim.h>
#include <scbrhf.h>
#include <ciinfo.h>
#include <spinfo.h>
#include <cbespn.h>
#include <infpri.h>
C
C     local variables:
C
      LOGICAL  NEWCI, NOCSF, CHCKCI
C
      CALL QENTER('SETCI2')
C
C
C CHCKCI returns true if new ci, either because it is first call
C or because specifications have changed.
C CHCKCI sets information in CIINFO and SPINFO
C Note that CHCKCI true if Hartree-Fock (i.e. nasht .le. 1) follows
C MCSCF without call of SETCI first.
C Note that change of symmetry will not change CI information.
C
      NEWCI = CHCKCI(NOSYM,1,LUPRI)
C
C
      IF ( NEWCI ) THEN
C     this routine is only called from ABAVAR and RSPVAR to set
C     NCONF for perturbation symmetry. If NEWCI then something is
C     inconsistent ....
         WRITE (LUPRI,'(/A/A)')
     &      ' SETCI2: Fatal programming error, new CI definition for',
     &      ' response equations compared to reference state!'
         CALL QUIT('Fatal error in SETCI2: new CI definition')
      END IF
C
      NOCSF = FLAG(27)
C
      IF ( NOCSF .OR. ICSYM .NE. LSYM .OR. TRIPLET) THEN
C        ... do not use CSF expansion if NOCSF or if symmetry
C            is not reference symmetry (CSF's may only be used
C            for reference symmetry in this version). /980820-hjaaj
         ICSF = 0
      ELSE
         ICSF = 1
      END IF
C
      IF (NASHT .GT. 1 .AND. .NOT. HSROHF) THEN
         IF (ICSF .EQ. 0) THEN
C           ... do not use CSF expansion
            NCONF  = NDTASM(ICSYM)
            NCONMA = MXNDT
         ELSE
            NCONF  = NCSASM(ICSYM)
            IF (NSYM .GT. 1) THEN
               NCONMA = MXNDT
C              ... we must have sufficient space
C                  for determinant vectors of different
C                  symmetries in ABACUS and RESPONSE
            ELSE
               NCONMA = MXNCS
            END IF
         END IF
      ELSE
C     ... NASHT .le. 1 (Hartree-Fock)
         IF (NASHT .EQ. 1) THEN
            JCSYM = IOPRHF
         ELSE IF (HSROHF) THEN
            JCSYM = LSYM
         ELSE 
            JCSYM = 1
         END IF
         IF (ICSYM .EQ. JCSYM) THEN
            NCONF = 1
         ELSE
            NCONF = 0
         END IF
         NCONMA = 1
      END IF
C
      CALL QEXIT('SETCI2')
      RETURN
      END
C  /* Deck chckci */
      LOGICAL FUNCTION CHCKCI(NOSYM,JPRINT,LUWCIX)
C
C 20-Aug-1998 Hans Joergen Aa. Jensen, 
C CHCKCI test extracted from SETCI
C
C Check if the CI specification has changed since last call
C to SETCI/SETCI2.
C
#include <implicit.h>
C
C Used from common blocks:
C   INFINP : NACTEL,ISPIN, MCTYPE, NELMN1, NELMX1, NELMN3, NELMX3
C   INFORB : NASHT,NAS1(*),NAS2(*),NAS3(*)
C   CIINFO : ICOMBI,IPSIGN, MXNDT, MXNCS
C   SPINFO : MULTS,MS2
C
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <ciinfo.h>
#include <spinfo.h>
C
C
      LOGICAL  REFSPC, NEWCI, NOCSF
      INTEGER  MAS1(8), MAS2(8), MAS3(8)
      SAVE     REFSPC, NEWCI, NOCSF, NOSYMSV,
     &         MACTEL,MAS1,MAS2,MAS3,MASHT,
     &         MELMN1,MELMX1,MELMN3,MELMX3,MCTOLD
      DATA     NEWCI /.TRUE./, MASHT/-1/
C
      CALL QENTER('CHCKCI')
C
C If not first call (NEWCI false) check if CI specification has
C changed, if yes then reset NEWCI to true.
C Note that Hartree-Fock (i.e. nasht .le. 1) will not change
C CI information.
      IF (NASHT .LE. 1 .OR. HSROHF) THEN
         IF (NASHT .NE. MASHT) THEN
            CHCKCI = .TRUE.
C           force SETCI2 to quit if HF follows MCSCF,
C           without intervening call of SETCI /980820-hjaaj
            IF (JPRINT .GT. 0) WRITE (LUWCIX,'(/A/A,2I10)')
     &         ' CHCKCI: new CI definition (now RHF)!',
     &         ' old and new NASHT  :',MASHT,NASHT
         ELSE
            CHCKCI = .FALSE.
         END IF
         MASHT = NASHT
         GO TO 8000
      END IF
C
C Note that change of symmetry will not change CI information.
C
      IF (.NOT.NEWCI) THEN
         IF (MACTEL .NE. NACTEL) NEWCI = .TRUE.
         IF (MCTOLD .NE. MCTYPE) NEWCI = .TRUE.
         IF (MCTYPE .EQ. 2) THEN
C           if RAS
            IF (MELMN1 .NE. NELMN1) NEWCI = .TRUE.
            IF (MELMX1 .NE. NELMX1) NEWCI = .TRUE.
            IF (MELMN3 .NE. NELMN3) NEWCI = .TRUE.
            IF (MELMX3 .NE. NELMX3) NEWCI = .TRUE.
         END IF
         DO 100 ISYM = 1,NSYM
            IF (MAS1(ISYM) .NE. NAS1(ISYM)) NEWCI = .TRUE.
            IF (MAS2(ISYM) .NE. NAS2(ISYM)) NEWCI = .TRUE.
            IF (MAS3(ISYM) .NE. NAS3(ISYM)) NEWCI = .TRUE.
  100    CONTINUE
         IF (NOCSF  .NEQV. FLAG(27)) NEWCI = .TRUE.
         IF (MULTS  .NE.   ISPIN)    NEWCI = .TRUE.
         IF (REFSPC .NEQV. FLAG(28)) NEWCI = .TRUE.
         IF (NOSYMSV.NE.   NOSYM   ) NEWCI = .TRUE.
         IF (NEWCI .AND. JPRINT .GT. 0) THEN
            WRITE (LUWCIX,'(/A,4(/A,2I10),2(/A,2L10))')
     &      ' CHCKCI: new CI definition!',
     &      ' old and new NACTEL :',MACTEL,NACTEL,
     &      ' old and new MCTYPE :',MCTOLD,MCTYPE,
     &      ' old and new ISPIN  :',MULTS ,ISPIN ,
     &      ' old and new NOSYM  :',NOSYMSV,NOSYM,
     &      ' old and new NOCSF  :',NOCSF ,FLAG(27),
     &      ' old and new REFSPC :',REFSPC,FLAG(28)
            IF (MCTYPE .EQ. 2) THEN
               WRITE (LUWCIX,'(4(/A,2I10))')
     &         ' old and new NELMN1 :',MELMN1,NELMN1,
     &         ' old and new NELMX1 :',MELMX1,NELMX1,
     &         ' old and new NELMN3 :',MELMN3,NELMN3,
     &         ' old and new NELMX1 :',MELMX3,NELMX3
               WRITE (LUWCIX,'(/A,8I5)')
     &         ' old NAS1(i) :',(MAS1(ISYM),ISYM=1,NSYM)
               WRITE (LUWCIX,'(A,8I5)')
     &         ' new NAS1(i) :',(NAS1(ISYM),ISYM=1,NSYM)
               WRITE (LUWCIX,'(/A,8I5)')
     &         ' old NAS3(i) :',(MAS3(ISYM),ISYM=1,NSYM)
               WRITE (LUWCIX,'(A,8I5)')
     &         ' new NAS3(i) :',(NAS3(ISYM),ISYM=1,NSYM)
            END IF
            WRITE (LUWCIX,'(/A,8I5)')
     &      ' old NAS2(i) :',(MAS2(ISYM),ISYM=1,NSYM)
            WRITE (LUWCIX,'(A,8I5/)')
     &      ' new NAS2(i) :',(NAS2(ISYM),ISYM=1,NSYM)
         END IF
      END IF
C
      IF (NEWCI) THEN
         NOCSF = FLAG(27)
C. Ms component chosen as lowest possible component for CSF's
CMAERKE-901019: make input option for MS2
         MULTS = ISPIN
         IF (NOCSF) THEN
            MS2 = MULTS - 1
         ELSE
c--renzo
            IF (DONEVPT) THEN
               MS2 = MULTS - 1
            ELSE
               MS2 = MOD(MULTS-1,2) !original statement
            ENDIF
c--renzo
         END IF
C
C Input test:
C
         NERR = 0
         IF (NELMN1 .GT. NELMX1) NERR = NERR + 1
         IF (NELMN3 .GT. NELMX3) NERR = NERR + 1
         NAS1T = 0
         NAS3T = 0
         DO 200 I = 1,NSYM
            NAS1T = NAS1T + NAS1(I)
            NAS3T = NAS3T + NAS3(I)
  200    CONTINUE
         IF (NELMN1 .GT. 2*NAS1T) NERR = NERR + 1
         IF (NELMN3 .GT. 2*NAS3T) NERR = NERR + 1
         IF (NERR .GT. 0) THEN
            WRITE (LUWCIX,'(//A/A,I5/)')
     *         ' Error in RAS specification in CHCKCI',
     *         ' Number of errors found:',NERR
            WRITE (LUWCIX,'(A,2I8)') ' NELMN1, NELMX1',NELMN1,NELMX1,
     *                               ' NELMN3, NELMX3',NELMN3,NELMX3,
     *                               ' NAS1T,  NAS3T ',NAS1T ,NAS3T
            CALL QUIT('Error in RAS specification discovered in CHCKCI')
         END IF
C
         MACTEL = NACTEL
         MCTOLD = MCTYPE
         MELMN1 = NELMN1
         MELMX1 = NELMX1
         MELMN3 = NELMN3
         MELMX3 = NELMX3
         DO 300 ISYM = 1,NSYM
            MAS1(ISYM) = NAS1(ISYM)
            MAS2(ISYM) = NAS2(ISYM)
            MAS3(ISYM) = NAS3(ISYM)
  300    CONTINUE
C
C        890129/hjaaj: in this version only ICOMBI=0 and IPSIGN=1.
C        ICOMBI and IPSIGN are in CIINFO
         ICOMBI = 0
         IPSIGN = 1
C
         REFSPC  = FLAG(28)
         NOSYMSV = NOSYM
      END IF
C
      MASHT  = NASHT
      CHCKCI = NEWCI
      NEWCI  = .FALSE.
C
 8000 CALL QEXIT('CHCKCI')
      RETURN
      END
C  /* Deck getdets */
      SUBROUTINE GETDETS(ICSYM,NCSF,NDET,XNDXCI,CVEC,DVEC,WRK,LWRK)
C
C May 2000 Hans Joergen Aa. Jensen
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XNDXCI(*), WRK(LWRK)
C
C   INFPRI : IPRCIX
C
#include <csfbas.h>
#include <ciinfo.h>
#include <infpri.h>
C
      IF (NCSF .NE. NCSASM(ICSYM) .OR. NDET .NE. NDTASM(ICSYM)) THEN
         WRITE (LUPRI,'(/A/A,I4,2(/A,2I12))')
     &   ' GETDETS ERROR:  NCSF and/or NDET not consistent:',
     &   ' ICSYM :',ICSYM,
     &   ' NCSF, NCSASM(ICSYM) :',NCSF, NCSASM(ICSYM),
     &   ' NDET, NDTASM(ICSYM) :',NDET, NDTASM(ICSYM)
         CALL QUIT('GETDETS: NCSF and/or NDET not consistent')
      END IF
      IF (LWRK .LT. NDET) CALL ERRWRK('GETDETS',NDET,LWRK)
      CALL DCOPY(NCSF,CVEC,1,WRK,1)
C     ... CSDTVC destroys CSF vector
      CALL CSDTVC(WRK,DVEC,1,XNDXCI(KDTOC),
     &            XNDXCI(KICTS(1)),ICSYM,0,IPRCIX)
      RETURN
      END
