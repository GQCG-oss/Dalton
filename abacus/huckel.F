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
C   /* Deck hucdrv */
      SUBROUTINE HUCDRV(WORK,LWORK,CMO)
C
C     K.Ruud, July-96
C
#include <implicit.h>
#include <maxorb.h>
C
      DIMENSION WORK(LWORK), CMO(NCMOT)
#include <infdim.h>
#include <inforb.h>
#include <huckel.h>
C
      CALL QENTER('HUCDRV')
      NROW = NBAST + NHUCBA
      NELMNT = NROW*(NROW + 1)/2
      NHUCMX = 0
      DO 10 ISYM = 1, NSYM
         NHUCMX = MAX(NHUCMX,NHUCAO(ISYM))
 10   CONTINUE
      NNHUCMX = NHUCMX*(NHUCMX + 1)/2
      N2HUCMX = NHUCMX*NHUCMX
      N2HUCBA = NHUCBA*NHUCBA
      NTRDIM = MAX(NHUCMX*NBASMA,NBASMA*(NBASMA + 1)/2)
C
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KTOTOV,NELMNT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHUCMT,NNHUCMX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHUCOV,N2HUCBA,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVC,N2HUCMX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTRANM,NTRDIM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVL,NHUCBA,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVS,NHUCBA,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KOVERL,NBASMA*NBASMA,WORK,KFREE,LFREE)
      CALL HUCDR1(WORK(KFREE),LFREE,WORK(KTOTOV),NELMNT,NROW,
     &            WORK(KHUCMT),WORK(KHUCOV),NHUCMX,WORK(KEIGVC),
     &            WORK(KTRANM),WORK(KEIGVL),WORK(KEIGVS),WORK(KOVERL),
     &            CMO)
      CALL MEMREL('HUCDRV',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('HUCDRV')
      RETURN
      END
C   /* Deck hucdr1 */
      SUBROUTINE HUCDR1(WORK,LWORK,TOTOVL,NELMNT,NROW,
     &                  HUCMAT,HUCOVL,NHUCMX,
     &                  EIGVEC,TRANMT,EIGVAL,EIGVLS,OVERLP,CMO)
C
C     Driver routine for solving the extended Huckel problem, and 
C     projecting the eigenvectors onto the complete basis set
C     
C     K.Ruud, July-96
C
#include <implicit.h>
#include <dummy.h>
#include <maxorb.h>
      LOGICAL ANTSYM, FNDLB2
      CHARACTER*8 RTNLBL(2)
      DIMENSION WORK(LWORK)
      DIMENSION TOTOVL(NELMNT), HUCMAT(*), HUCOVL(*), TRANMT(*),
     &          EIGVEC(NHUCMX*NHUCMX), CMO(NCMOT),
     &          EIGVAL(NHUCBA), EIGVLS(NHUCBA), OVERLP(NBASMA*NBASMA)
C
#include <priunit.h>
#include <inftap.h>
#include <infdim.h>
#include <inforb.h>
#include <infpri.h>
#include <scbrhf.h>
#include <huckel.h>
C
D     p6flag(6) = .true.
      IF (P6FLAG(6) .OR. IPRI6.GE.0 .OR. IPRRHF.GE.0 .OR. EWMO)
     &   WRITE (LUPRI,'(/A/2(/5X,A,L10)/5X,A,8I5)')
     &      ' *** Output from Huckel module :',
     &      'Using EWMO model: ',EWMO,
     &      'Using EHT  model: ',.NOT.EWMO,
     &      'Number of Huckel orbitals each symmetry:',
     &      (NHUCAO(I),I=1,NSYM)
      IF (EWMO) THEN
         WRITE (LUPRI,'(/A/A/A)')
     &   ' EWMO - Energy Weighted Maximum Overlap - is a Huckel type'//
     &   ' method,',
     &   ' which normally is better than Extended Huckel Theory.',
     &   ' Reference: Linderberg and Ohrm, Propagators in Quantum'//
     &   ' Chemistry (Wiley, 1973)'
      END IF
C
      NNHUCBA = NHUCBA * (NHUCBA + 1) / 2
C
C     Read AO overlap matrix from AOPROPER file
C
      CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      IF (EWMO) THEN
C
C     ... read EWMO overlap matrix:
C
         REWIND LUPROP
         IF (FNDLB2('HUCKOVLP',RTNLBL,LUPROP)) THEN
            CALL READT(LUPROP,NELMNT,TOTOVL)
            IF (P6FLAG(6)) THEN
               CALL HEADER(
     &         'HUCDR1: HUCKOVLP matrix read from AOPROPER',-1)
               CALL OUTPAK(TOTOVL,NROW,1,LUPRI)
            END IF
         ELSE
            WRITE (LUPRI,'(//3A)') ' >>> HUCDR1 error: "HUCKOVLP"'//
     *         ' matrix not found on AOPROPER.'
            CALL QUIT('HUCDR1 error: HUCKOVLP not found on AOPROPER')
         END IF
C
C     ... extract EWMO overlap matrix from TOTOVL:
C
         CALL DZERO(HUCOVL,NNHUCBA)
         ISTRT  = 0
         NHUCT  = 0
         JPOS   = 1
C
         DO ISYM = 1, NSYM
            NHAOSY = NHUCAO(ISYM)
         IF (NHAOSY .GT. 0) THEN
            NBASI  = NBAS(ISYM)
            IBASTO = IBAS(ISYM) + NHUCT
            IHUCI  = IBASTO + NBASI
            NORBI  = NORB(ISYM)
            IF (NORBI .EQ. 0) THEN
               WRITE (LUPRI,'(/A,I3/I5,A)')
     &            'HUCDRV error in symmetry',ISYM,
     &            NHAOSY,' Huckel orbitals, but no orbitals!'
               CALL QUIT('HUCDRV error, NHAOSY.ne.0 .and. NORBI.eq.0')
            END IF
            ISTRT  = (IHUCI*(IHUCI + 1))/2
C
C           Copy the Huckel overlap matrix into a lower triangular matrix,
C           where the elements are sorted in consecutive order
C
            IPOS = JPOS
            DO I = 1, NHAOSY
               CALL DCOPY(I,TOTOVL(ISTRT + IHUCI + 1),1,HUCOVL(IPOS),1)
               IPOS = IPOS + I
               ISTRT = ISTRT + (IHUCI + I)
            END DO   ! I = 1,NHAOSY
            IF (P6FLAG(6) .OR. IPRRHF.GE.6) THEN
               WRITE (LUPRI,'(/A,I3,/A)')
     &            ' Extracted EWMO overlap matrix of symmetry :',ISYM,
     &            ' ------------------------------------------------'
               CALL OUTPAK(HUCOVL(JPOS),NHAOSY,1,LUPRI)
            END IF
            JPOS = IPOS
            NHUCT = NHUCT + NHAOSY
         END IF
         END DO   ! ISYM = 1,NSYM
      END IF
C
C     Read Huckel matrix, both for EWMO and EHT
C
      REWIND LUPROP
      IF (FNDLB2('HUCKEL  ',RTNLBL,LUPROP)) THEN
         CALL READT(LUPROP,NELMNT,TOTOVL)
         IF (P6FLAG(6)) THEN
            CALL HEADER(
     &      'HUCDR1: Huckel matrix read from AOPROPER',-1)
            CALL OUTPAK(TOTOVL,NROW,1,LUPRI)
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> HUCDR1 error: "HUCKEL  "'//
     *      ' matrix not found on AOPROPER.'
         CALL QUIT('HUCDR1 error: HUCKEL overlap not found on AOPROPER')
      END IF
      CALL GPCLOSE(LUPROP,'KEEP')
C
C     Save initial CMO for HUCORT
C
      KCMO1 = 1
      KWRK  = KCMO1 + NCMOT
      LWRK  = LWORK - (KWRK-1)
      CALL DCOPY(NCMOT,CMO,1,WORK(KCMO1),1)
C
C
C     Extract Huckel matrix, do one symmetry at a time
C
      ISTRT  = 0
      NHUCT  = 0
      JPOS   = 1
C
      DO 100 ISYM = 1, NSYM
         NHAOSY = NHUCAO(ISYM)
         IF (NHAOSY .EQ. 0) GO TO 100
         NBASI  = NBAS(ISYM)
         IBASTO = IBAS(ISYM) + NHUCT
         IHUCI  = IBASTO + NBASI
         NORBI  = NORB(ISYM)
         IF (NORBI .EQ. 0) THEN
            WRITE (LUPRI,'(/A,I3/I5,A)')
     &         'HUCDRV error in symmetry',ISYM,
     &         NHAOSY,' Huckel orbitals, but no orbitals!'
            CALL QUIT('HUCDRV error, NHAOSY.ne.0 .and. NORBI.eq.0')
         END IF
         ISTRT  = (IHUCI*(IHUCI + 1))/2
         ICSYM  = ICMO(ISYM) + 1
C
C     Extract ordinary overlap matrix
C
         IOVLST = (IBASTO*(IBASTO + 3))/2 + 1
         IPOS = 1
         DO 115 I = 1, NBASI
            CALL DCOPY(I,TOTOVL(IOVLST),1,TRANMT(IPOS),1)
            IPOS = IPOS + I
            IOVLST = IOVLST + (IBASTO + I)
 115     CONTINUE
         CALL DSPTSI(NBASI,TRANMT,OVERLP)
C
C     Copy the Huckel part into a lower triangular matrix, where the elements
C     are sorted in consecutive order
C
C     Also extract the transformation matrix from HUCKEL basis to full basis
C
         IPOS = 1
         DO 110 I = 1, NHAOSY
            CALL DCOPY(I,TOTOVL(ISTRT + IHUCI + 1),1,
     &                 HUCMAT(IPOS),1)
            CALL DCOPY(NBASI,TOTOVL(ISTRT + IBASTO + 1),1,
     &                 TRANMT((I - 1)*NBASI + 1),1)
            IPOS = IPOS + I
            ISTRT = ISTRT + (IHUCI + I)
 110     CONTINUE 
         IF (P6FLAG(6) .OR. IPRRHF.GE.8) THEN
            WRITE (LUPRI,'(/A,I3,/A)')
     &           ' Extracted Huckel matrix of symmetry :',ISYM,
     &           ' ------------------------------------------------'
            CALL OUTPAK(HUCMAT,NHAOSY,1,LUPRI)
         END IF
         IF (P6FLAG(6) .OR. IPRRHF.GE.15) THEN
            WRITE (LUPRI,'(/A,I3,/A)')
     &           ' Extracted transformation matrix of symmetry :',ISYM,
     &           ' ------------------------------------------------'
            CALL OUTPUT(TRANMT,1,NBASI,1,NHAOSY,NBASI,NHAOSY,1,LUPRI)
         END IF
C
C        EWMO or EHT ?
C
        IF (EWMO) THEN
C        ... call EISPACK routine RSG
         N2HAOSY = NHAOSY*NHAOSY
         IPOS = 1 + (NHUCT*(NHUCT + 1))/2 + NHUCT
         KHMAT = KWRK
         KSMAT = KHMAT+ N2HAOSY
         KFV1  = KSMAT+ N2HAOSY
         KFV2  = KFV1 + NHAOSY
         KLAST = KFV2 + NHAOSY
         IF (KLAST .GT. LWORK) CALL STOPIT('HUCDRV','RSG',KLAST,LWORK)
         CALL DSPTSI(NHAOSY,HUCMAT,WORK(KHMAT))
         CALL DSPTSI(NHAOSY,HUCOVL(JPOS),WORK(KSMAT))
         IF (P6FLAG(6) .OR. IPRRHF.GE.8) THEN
            WRITE (LUPRI,*) 'EWMO AO Hamiltonian matrix W'//
     &                      ' corresponding to separate atoms'
            CALL OUTPUT(WORK(KHMAT),1,NHAOSY,1,NHAOSY,
     &                  NHAOSY,NHAOSY,1,LUPRI)
         END IF
         CALL DGEMM('N','N',NHAOSY,NHAOSY,NHAOSY,1.D0,
     &              WORK(KHMAT),NHAOSY, WORK(KSMAT),NHAOSY,0.D0,
     &              EIGVEC,NHAOSY)
         CALL DGEMM('N','N',NHAOSY,NHAOSY,NHAOSY,1.D0,
     &              WORK(KSMAT),NHAOSY, EIGVEC,NHAOSY,0.D0,
     &              WORK(KHMAT),NHAOSY)
C        HMAT = SWS , SMAT = S
         IF (P6FLAG(6) .OR. IPRRHF.GE.4) THEN
            WRITE (LUPRI,*) 'EWMO Hamiltonian matrix SWS'
            CALL OUTPUT(WORK(KHMAT),1,NHAOSY,1,NHAOSY,
     &                  NHAOSY,NHAOSY,1,LUPRI)
            WRITE (LUPRI,*) 'EWMO overlap matrix S'
            CALL OUTPUT(WORK(KSMAT),1,NHAOSY,1,NHAOSY,
     &                  NHAOSY,NHAOSY,1,LUPRI)
         END IF
         JPOS = JPOS + ( NHAOSY * (NHAOSY+1) ) / 2
         CALL RSG(NHAOSY,NHAOSY,WORK(KHMAT),WORK(KSMAT),
     &            EIGVAL(NHUCT+1),1,EIGVEC,WORK(KFV1),WORK(KFV2),IERR)
         IF ( IERR.NE.0 ) THEN
            WRITE(LUPRI,'(/A,I10)')
     *      ' HUCDRV EWMO eigenvalue problem not converged, code: ',IERR
            CALL QUIT('HUCDRV: EWMO eigenvalue equation not converged')
         END IF
C           Save orbital symmetry in EIGVLS
            DO I = 1, NHAOSY
               EIGVLS(NHUCT + I) = ISYM
            END DO
        ELSE
C
C       EHT: Diagonalize Huckel matrix (Note: KLAST is not a misprint)
C
         KSCRA1 = KWRK
         KSCRA2 = KSCRA1 + NBASI
         KLAST  = KSCRA1 + NBASI*NHAOSY
         IF (KLAST .GT. LWORK) CALL STOPIT('HUCDRV','JACO',KLAST,LWORK)
         CALL DUNIT(EIGVEC,NHAOSY)
         CALL JACO(HUCMAT,EIGVEC,NHAOSY,NHAOSY,NHAOSY,WORK(KSCRA1),
     &             WORK(KSCRA2))
C
C     Sort the eigenvalues and eigenvectors
C     Save orbital symmetry in EIGVLS
C
         II = 0
         DISYM = ISYM
         DO 120 I = 1, NHAOSY
            II = II + I
            EIGVAL(NHUCT + I) = HUCMAT(II)
            EIGVLS(NHUCT + I) = DISYM
 120     CONTINUE 
         CALL ORDER(EIGVEC,EIGVAL(NHUCT + 1),NHAOSY,NHAOSY)
        END IF
C
C
C     We write the eigenvalues if high enough print level
C     and also the vectors if even higher print level
C
         IF (P6FLAG(6) .OR. IPRI6.GE.0 .OR. IPRRHF.GE.0) THEN
            IF (EWMO) THEN
               WRITE (LUPRI,'(/A,I3)') ' Huckel EWMO eigenvalues for'//
     &           ' symmetry :',ISYM
            ELSE
               WRITE (LUPRI,'(/A,I3)') ' Huckel EHT eigenvalues for'//
     &           ' symmetry :',ISYM
            END IF
            WRITE (LUPRI,'(5X,5F15.6)') (EIGVAL(NHUCT+I), I = 1, NHAOSY)
         END IF
         IF (P6FLAG(6) .OR. IPRRHF.GE.2) THEN
            IF (EWMO) THEN
               WRITE (LUPRI,'(/A,I3)') ' Huckel EWMO eigenvectors for'//
     &           ' symmetry :',ISYM
            ELSE
               WRITE (LUPRI,'(/A,I3)') ' Huckel EHT eigenvectors for'//
     &           ' symmetry :',ISYM
            END IF
            CALL OUTPUT(EIGVEC,1,NHAOSY,1,NHAOSY,NHAOSY,NHAOSY,1,LUPRI)
         END IF
         NHUCT = NHUCT + NHAOSY
C         
C     Project eigenvectors onto complete basis
C
C     We do this as (S_{NN})^{-1}*(S_{NO}*C_{O})
C     N and O refers to new (large) and old (small) basis respectively
C     For HUCKEL basis to HUCKEL basis this gives the identity matrix
C
C     Explicit inversion may take some time.....
C
         KSCRA1 = KWRK
         KSCRA2 = KSCRA1 + NBASI
         KLAST  = KSCRA1 + NBASI*NHAOSY
         IF (KLAST.GT.LWORK) CALL STOPIT('HUCDRV','DGEINV',KLAST,LWORK)
         CALL DGEINV(NBASI,OVERLP,OVERLP,WORK(KSCRA1),WORK(KSCRA2),INFO)
         IF (INFO .NE. 0) THEN
            WRITE (LUPRI,'(/A/A)') ' Singular overlap matrix HUCDRV',
     &           ' Huckel guess aborted'
            CALL QUIT('Singular overlap matrix when projecting '//
     &                'in Huckel')
         END IF
         CALL DGEMM('N','N',NBASI,NHAOSY,NBASI,1.D0,
     &              OVERLP,NBASI, TRANMT,NBASI,0.D0, WORK(KSCRA1),NBASI)
         CALL DGEMM('N','N',NBASI,NHAOSY,NHAOSY,1.D0,
     &              WORK(KSCRA1),NBASI, EIGVEC,NHAOSY,0.D0,
     &              CMO(ICSYM),NBASI)
 100  CONTINUE 
      IF (P6FLAG(6)) THEN
         WRITE (LUPRI,'(/A)')
     &   ' Huckel guess MOs before Gram-Schmidt orthogonalization:'
         CALL PRORB(CMO,.FALSE.,LUPRI)
      END IF
C
C     We do our own Gram-Schmidt Orthogonalization of the trial vectors
C
      CALL HUCORT(CMO,WORK(KCMO1),TOTOVL,WORK(KWRK),LWRK)
      IF (P6FLAG(6)) THEN
         WRITE (LUPRI,'(/A)')
     &   ' Huckel guess MOs after Gram-Schmidt orthogonalization:'
         CALL PRORB(CMO,.FALSE.,LUPRI)
      END IF
C
C     Determine Hartree-Fock occupation if requested
C
      IF (AUTOCC) THEN
         CALL ORDER (EIGVLS,EIGVAL,NHUCT,1)
         CALL IZERO(NISH,8)
         MOCC = NRHFEL/2
         DO 98 IK = 1, MOCC
            ISYM = NINT(EIGVLS(IK))
            NISH(ISYM) = NISH(ISYM) + 1
 98      CONTINUE
         CALL ICOPY(8,NISH,1,NOCC,1)
         CALL ICOPY(8,NISH,1,NRHF,1)
         IF (2*MOCC .NE. NRHFEL) THEN
            IOPRHF = NINT(EIGVLS(MOCC + 1))
            LSYM   = IOPRHF
            CALL IZERO(NASH,8)
            NASH(IOPRHF) = 1
            NOCC(IOPRHF) = NOCC(IOPRHF) + 1
         END IF
      END IF
C
C     We are done
C
      RETURN
      END
C   /* Deck hucfun */
      SUBROUTINE HUCFUN(NQ,IQCORE,IQM,ISHL1,HUCTER)
C
C     Ascribed orbital-ionization potential to a shell
C     Ionization potentials supplied from TURBOMOLE basis sets.
C
C     Input:
C     NQ     = atomic number
C     IQCORE = number of core electrons described with ECP
C     IQM    = l_quantum_number + 1
C     ISHL1  = orbital no. for this IQM
C
C     Output:
C     HUCTER = orbital energy
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
C
      DIMENSION EIGVL1(10,3), EIGVL2(11:18,5), EIGVL3(19:30,7),
     &          EIGVL4(31:36,8)
C
#include <huckel.h>
C
C           1s,2s,2p
            DATA ((EIGVL1(I,J),J=1,3),I=1,10) /
     &  -.49927730774995D0,0.0D0,0.0D0,
     &  -.91412243731124D0,0.0D0,0.0D0,
     &  -.24816746159666D+01,-.19440636753810D+00,0.0D0,
     &  -.47343786810965D+01,-.30862552384615D+00,0.0D0,
     &  -.76999577773441D+01,-.49234779827247D+00,-.30891241897079D+00,
     &  -.11330261219101D+02,-.70027915640211D+00,-.43091547364486D+00,
     &  -.15633247432629D+02,-.93724789614141D+00,-.56339961843608D+00,
     &  -.20673892113414D+02,-.12305199155485D+01,-.62371203325220D+00,
     &  -.26386066585568D+02,-.15540499409522D+01,-.71712342746009D+00,
     &  -.32778774575442D+02,-.19030293335186D+01,-.83418178868235D+00/
C           1s,2s,3s,2p,3p
      DATA ((EIGVL2(I,J),J=1,5),I=11,18) /
     &  -.40492097249337D+02,-.27894059577975D+01,-.18075821513267D+00,
     &  -.15202747041896D+01,0.0D0,
     &  -.49042962436482D+02,-.37522466129384D+01,-.25235717836570D+00,
     &  -.22826705448876D+01,0.0D0,
     &  -.58496434210500D+02,-.48998155658682D+01,-.38935511002256D+00,
     &  -.32095334365824D+01,-.20696626740742D+00,
     &  -.68812636968677D+02,-.61475051057187D+01,-.53622559159139D+00,
     &  -.42501446429514D+01,-.29430170913097D+00,
     &  -.79976261544148D+02,-.75043595509827D+01,-.69337872525626D+00,
     &  -.53987425814700D+01,-.38853865472178D+00,
     &  -.92011193529873D+02,-.89976454800258D+01,-.87537887942574D+00,
     &  -.66810815897629D+01,-.43238675379853D+00,
     &  -.10488880701021D+03,-.10599513103287D+02,-.10673569335135D+01,
     &  -.80700723698943D+01,-.49937394834322D+00,
     &  -.11860449689801D+03,-.12318825848607D+02,-.12712490599994D+01,
     &  -.95673899752542D+01,-.58224605996824D+00/
C           1s,2s,3s,4s,2p,3p,3d
      DATA ((EIGVL3(I,J),J=1,7),I=19,30) /
     &  -.13353069098600D+03,-.14485916490245D+02,-.17471301687431D+01,
     &  -.14708801737535D+00,-.11518191755577D+02,-.95336656131597D+00,
     &  0.0D0,
     &  -.14936143935251D+03,-.16818757548728D+02,-.22439704960253D+01,
     &  -.19524406319691D+00,-.13628535450378D+02,-.13399654223312D+01,
     &  0.0D0,
     &  -.16589683256107D+03,-.19075457745884D+02,-.25647622232670D+01,
     &  -.20941844395422D+00,-.15666732040035D+02,-.15728195654682D+01,
     &  -.34140154000383D+00,
     &  -.18327961760075D+03,-.21428535937406D+02,-.28783252991479D+01,
     &  -.22086095568251D+00,-.17800707638882D+02,-.18001814645575D+01,
     &  -.40228838660405D+00,
     &  -.20149843763423D+03,-.23867228240440D+02,-.31785445885957D+01,
     &  -.22932792514940D+00,-.20019592751106D+02,-.20157586107128D+01,
     &  -.50522452423977D+00,
     &  -.22058703568056D+03,-.26430356269507D+02,-.34932183673093D+01,
     &  -.23818021032623D+00,-.22361495552450D+02,-.22443021673275D+01,
     &  -.56349660940991D+00,
     &  -.24052800159493D+03,-.29099297724807D+02,-.38095159549645D+01,
     &  -.24603815015926D+00,-.24807928660070D+02,-.24739509857674D+01,
     &  -.63191681323148D+00,
     &  -.26136607277522D+03,-.31923241649121D+02,-.41603189279520D+01,
     &  -.25579970324891D+00,-.27407488016598D+02,-.27348524224893D+01,
     &  -.63779115886163D+00,
     &  -.28307236020005D+03,-.34870580572438D+02,-.45238526025025D+01,
     &  -.26527761764716D+00,-.30128927239222D+02,-.30070645013762D+01,
     &  -.64055221022396D+00,
     &  -.30561629584843D+03,-.37909005403997D+02,-.48795454632899D+01,
     &  -.27316053356929D+00,-.32940079186849D+02,-.32712974970763D+01,
     &  -.68290567072958D+00,
     &  -.32877476592812D+03,-.40792091501464D+02,-.49874527579018D+01,
     &  -.23137145006834D+00,-.35599071176379D+02,-.33030501793165D+01,
     &  -.46685257058320D+00,
     &  -.35329114905617D+03,-.44339943167914D+02,-.56197870069515D+01,
     &  -.28802215775704D+00,-.38911536374417D+02,-.38240801472867D+01,
     &  -.76393436766955D+00 /
C           1s,2s,3s,4s,2p,3p,4p,3d
      DATA ((EIGVL4(I,J),J=1,8),I=31,36) /
     &  -.37880682397662D+03,-.48148799963237D+02,-.63796010909107D+01,
     &  -.42021171960591D+00,-.42483270039189D+02,-.44701987106119D+01,
     &  -.20559461561373D+00,-.11773190246274D+01,
     &  -.40523389984477D+03,-.52131809639667D+02,-.71777951435091D+01,
     &  -.54906911734569D+00,-.46226995190672D+02,-.51513035355167D+01,
     &  -.28393862239525D+00,-.16204742000157D+01,
     &  -.43257598995361D+03,-.56291604153201D+02,-.80173521238562D+01,
     &  -.68144923597635D+00,-.50145406172146D+02,-.58713989315709D+01,
     &  -.36566953118867D+00,-.20990510281660D+01,
     &  -.46085658870671D+03,-.60649930780960D+02,-.89196551027896D+01,
     &  -.83199435553583D+00,-.54260396193709D+02,-.66521513180541D+01,
     &  -.39773888567469D+00,-.26357881872702D+01,
     &  -.49004878507123D+03,-.65180120390772D+02,-.98590440139025D+01,
     &  -.98640692417606D+00,-.58545396597353D+02,-.74685601227765D+01,
     &  -.45079853997024D+00,-.32059197938172D+01,
     &  -.52015312537666D+03,-.69882264742765D+02,-.10836105147812D+02,
     &  -.11458450497054D+01,-.63000574281761D+02,-.83214809593102D+01,
     &  -.51687272108885D+00,-.38104800526990D+01/
C
      CALL QENTER('HUCFUN')
      CALL ECPLCORE(IQCORE,IQM,ILOFF)
      ISHL = ISHL1 + ILOFF
C     ... calculate off-set for this l_quantum_number for ECP,
C         where the lowest orbitals are excluded.
C         (If not ECP, i.e. IQCORE .eq. 0, then ILOFF will be zero)
C
      IF (NQ .LE. 0) THEN
C           Output corrected (WK/UniKA/04-11-2002).
            WRITE (LUPRI,'(/A,I10/A)')
     &           ' WARNING: No Huckel guess for charge of',NQ,
     &           ' Hueckel start guess replaced'//
     &           ' with a diagonalization of H1'
            DOHUCKEL = .FALSE.
      ELSE IF (NQ .LE. 10) THEN
         I = 2*(IQM - 1) + ISHL 
C        I: 1s,2s,2p
         IF (I .LE. 3) THEN
            HUCTER = EIGVL1(NQ,I)
         ELSE
            WRITE (LUPRI,'(/A,I2,A,I3,/A)')
     &           ' WARNING: quantum number ',I,' not recognised for '//
     &           'nucleus :',NQ, ' Hueckel start guess replaced'//
     &           ' with a diagonalization of H1'
            DOHUCKEL = .FALSE.
         END IF
      ELSE IF (NQ .GT. 10 .AND. NQ .LE. 18) THEN
         IF (IQM .EQ. 1) THEN
            I = ISHL
         ELSE IF (IQM .EQ. 2) THEN
            I = 3 + ISHL
         ELSE
            I = 6
         END IF
C        I: 1s,2s,3s,2p,3p
         IF (I .LE. 5) THEN
            HUCTER = EIGVL2(NQ,I)
         ELSE
            WRITE (LUPRI,'(/A,I2,A,I3,/A)')
     &           ' WARNING: quantum number ',I,' not recognised for '//
     &           'nucleus :',NQ, ' Hueckel start guess replaced'//
     &           ' with a diagonalization of H1'
            DOHUCKEL = .FALSE.
         END IF
      ELSE IF (NQ .GT. 18 .AND. NQ .LE. 30) THEN
         IF (IQM .EQ. 1) THEN
            I = ISHL
         ELSE IF (IQM .EQ. 2) THEN
            I = 4 + ISHL
         ELSE IF (IQM .EQ. 3) THEN
            I = 6 + ISHL
         ELSE
            I = 8
         END IF
C        I: 1s,2s,3s,4s,2p,3p,3d
         IF (I .LE. 7) THEN
            HUCTER = EIGVL3(NQ,I)
         ELSE
            WRITE (LUPRI,'(/A,I2,A,I3,/A)')
     &           ' WARNING: quantum number ',I,' not recognised for '//
     &           'nucleus :',NQ, ' Hueckel start guess replaced'//
     &           ' with a diagonalization of H1'
            DOHUCKEL = .FALSE.
         END IF
      ELSE IF (NQ .GT. 30 .AND. NQ .LE. 36) THEN
         IF (IQM .EQ. 1) THEN
            I = ISHL
         ELSE IF (IQM .EQ. 2) THEN
            I = 4 + ISHL
         ELSE IF (IQM .EQ. 3) THEN
            I = 7 + ISHL
         ELSE
            I = 9
         END IF
C        I: 1s,2s,3s,4s,2p,3p,4p,3d
         IF (I .LE. 8) THEN
            HUCTER = EIGVL4(NQ,I)
         ELSE
            WRITE (LUPRI,'(/A,I2,A,I3,/A)')
     &           ' WARNING: quantum number ',I,' not recognised for '//
     &           'nucleus :',NQ, ' Hueckel start guess replaced'//
     &           ' with a diagonalization of H1'
            DOHUCKEL = .FALSE.
         END IF
      ELSE
         WRITE (LUPRI,'(/A,I3,/A)') ' No orbital ionization potential'//
     &        ' for nucleus with charge ', NQ, ' Hueckel start guess '//
     &        'replaced with a diagonlization of H1'
         DOHUCKEL = .FALSE.
      END IF
      CALL QEXIT('HUCFUN')
      RETURN
      END
C   /* Deck hucort */
      SUBROUTINE HUCORT(CMO,CMO1,S,WORK,LWORK)
C
C     Written by K.Ruud, aug-96
C     Version written on the basis of the ORTHO and NORM subroutines
C     of hjj
C     Revised by hjaaj May 2000 to handle round-off errors
C
#include <implicit.h>
#include <maxorb.h>
      PARAMETER (THNORM = 1.0D-6, D1 = 1.0D0, D0 = 0.0D0, THRRND=0.9D0)
      DIMENSION CMO(*), CMO1(*), S(*), WORK(LWORK)
C
#include <priunit.h>
#include <huckel.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
C
      CALL QENTER('HUCORT')
C
C     Read overlap matrix from LUONEL
C
      CALL RDONEL('OVERLAP ',.TRUE.,S,NNBAST)
      IF (P6FLAG(6)) THEN
         WRITE(LUPRI,'(/A)')
     *     ' (HUCORT) Overlap between AO basis functions :'
         CALL OUTPKB(S,NBAS,NSYM,-1,LUPRI)
      END IF
C
C     We orthogonalize orbitals symmetry by symmetry
C
      ISTS   = 1
      ISTC   = 1
      DO 100 ISYM = 1, NSYM
         NORBI = NORB(ISYM)
         NBASI = NBAS(ISYM)
         NHUCI = NHUCAO(ISYM)
         IF (NORBI .EQ. 0) GO TO 101
C
C     We orthogonalize the projected Huckel eigenvectors first
C
         CALL NORM(S(ISTS),CMO(ISTC),NBASI,NHUCI,WORK,THNORM,IRETUR)
         IF (IRETUR.NE.0) GO TO 5000
C
C     Now we add vectors that span the entire space and orthogonalize 
C     them to the previous vectors
Chj-may2000: modified to use the previous set of vectors in CMO1. The
C     original algorithm here could reinsert CMO vectors deleted by
C     DELMO or otherwise.
C     The initial vectors in CMO1 are generated by H1MO, so the NHUCI
C     first will probably have big overlap with the Huckel generated
C     orbitals. Therefore we begin with no. NHUCI+1 and wrap to 1
C     to end with NHUCI.
C
         IKROW1 = NHUCI
         DO 103 IKROW = NHUCI + 1, NORBI
            IVKC  = ISTC + (IKROW -1)*NBASI
 107        IKROW1 = IKROW1 + 1
            IF (IKROW1 .GT. NORBI) IKROW1 = 1
            IVKC1 = ISTC + (IKROW1-1)*NBASI
            CALL DCOPY(NBASI,CMO1(IVKC1),1,CMO(IVKC),1)
C
C           normalize initial vector C_k and find W_k = S*C_k
C
            ITURN = 0
 108        CONTINUE 
            CALL MPAPV(NBASI,S(ISTS),CMO(IVKC),WORK)
            TNORM = DDOT(NBASI,CMO(IVKC),1,WORK,1)
            IF (TNORM .LT. THNORM) THEN
Chj         ... zero vector on input
               IRETUR = -IKROW
               GO TO 5000
            END IF
Chj may2000: normalize input vector
C        (we ignore round-off errors as it is renormalized later)
            TNORM = D1 / SQRT(TNORM)
            CALL DSCAL(NBASI,TNORM,CMO(IVKC),1)
            CALL MPAPV(NBASI,S(ISTS),CMO(IVKC),WORK)
            TNORM = DDOT(NBASI,CMO(IVKC),1,WORK,1)
C
            IVJC = ISTC
            DO J = 1, IKROW-1
               T = DDOT(NBASI,CMO(IVJC),1,WORK,1)
               TNORM = TNORM - T*T
               IVJC = IVJC + NBASI
               WORK(NBASI+J) = -T
            END DO
            IF (TNORM .LT. THNORM) THEN
C
C              Zero norm,
C              we try the next mo vector from the initial mo's in CMO1
C
               IF (IKROW1 .EQ. NHUCI) THEN
C                 ... we have tried them all!
                  IRETUR = -10 000*IKROW
                  GO TO 5000
               END IF
               GOTO 107
            END IF
C
            IF (TNORM .LT. THRRND) THEN
               IF (ITURN .GT. 1) THEN
                  IRETUR = 10 000*IKROW
                  GO TO 5000
               END IF
               IROUND = 1
            ELSE
               IROUND = 0
            END IF
C
            TNORM = D1/SQRT(TNORM)
            DO J = 1, IKROW-1
               WORK(NBASI+J) = WORK(NBASI+J)*TNORM
            END DO
            WORK(NBASI+IKROW) = TNORM
            CALL DGEMM('N','N',NBASI,1,IKROW,1.D0,
     &                 CMO(ISTC),NBASI,
     &                 WORK(NBASI+1),IKROW,0.D0,
     &                 WORK,NBASI)
            CALL DCOPY(NBASI,WORK,1,CMO(IVKC),1)
C
            IF (IROUND .EQ. 1) THEN
               ITURN = ITURN + 1
               GO TO 108
            END IF
C
 103     CONTINUE 
 101     ISTS   = ISTS   + NBASI*(NBASI + 1)/2
         ISTC   = ISTC   + NORBI*NBASI
 100  CONTINUE 
      CALL QEXIT('HUCORT')
      RETURN
C
 5000 CONTINUE
      WRITE(LUW4,4010) ISYM, IRETUR
      IF (LUPRI.NE.LUW4)
     *   WRITE(LUPRI,4010) ISYM, IRETUR
      WRITE(LUERR,4010) ISYM, IRETUR
      CALL QUIT('*** ERROR *** FATAL ERROR IN HUCORT')
 4010 FORMAT(/' *** HUCORT-FATAL ERROR *** Linear dependency in',
     *        ' symmetry (Huckel part) =',I3/' CODE =',I10)
      END
