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
C  /* Deck rlxinp */
      SUBROUTINE RLXINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 5)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <abainf.h>
#include <cbirlx.h>
      DATA TABLE /'.SKIP  ', '.PRINT ','.SYMTES','.NOSELL','.STOP  '/
C
      NEWDEF = (WORD .EQ. '*RELAX')
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
                     GO TO (1,2,3,4,5), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in RLXINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in RLXINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  SYMTST = .TRUE.
               GO TO 100
    4             DOSELL = .FALSE.
               GO TO 100
    5             CUT    = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in RLXINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in RLXINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for RELAX:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' RELAX skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in RELAX :',IPRINT
            END IF
            IF (SYMTST) WRITE (LUPRI,'(2A)')
     *                    ' Calculation of all elements',
     *                    ' (i,j) and (j,i) in relaxation Hessian '
            IF (.NOT.DOSELL) WRITE (LUPRI,'(A,/,2A,/)')
     *           ' Sellers'' method for quadratic errors not used.',
     *           ' Errors in relaxation may be linear in threshold',
     *           ' for response equations.'
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after RELAX.'
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck rlxini */
      SUBROUTINE RLXINI
C
C     Initialize /CBIRLX/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbirlx.h>
C
      IPRINT = IPRDEF
      SKIP   = .NOT. (MOLHES .OR. DIPDER .OR. POLAR .OR. MAGSUS .OR.
     &                SHIELD .OR. SPNSPN .OR. VCD .OR. SPINRO .OR.
     &                MOLGFA .OR. QPGRAD)
      CUT    = .FALSE.
      SYMTST = .FALSE.
      DOSELL = .TRUE.
      RETURN
      END
C  /* Deck relax */
      SUBROUTINE RELAX(SPNPSO,SPNSD,SPNFC,SPSDFC,WORK,LWORK,PASS,DOTRPL)
C
C     Written 23-jan-1985 Hans Joergen Aa. Jensen
C     Modified 14-jun-1985 tuh
C     Rewritten to include dip. grad. and polarizabilities Jan 1990 tuh
C
C     Purpose: Calculate response contributions to molecular Hessian
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <iratdef.h>
#include <dummy.h>
C
      LOGICAL PASS, OLDDX, DOREAL, DOIMAG, DOSEL, DOTRPL, SELLER
      DIMENSION SPNPSO(*),SPNSD(*), SPNFC(*), SPSDFC(*), WORK(LWORK)

#ifdef VAR_G95
      real(8), allocatable :: fc_contribution(:, :)
!bugfix for the G95 compiler May 2011,
!should be removed again when G95 compiler fixed.
!  from the Dalton svn log:
!=========================================================
!r10285 | bast | 2011-05-14 16:28:15 +0200 (Sat, 14 May 2011) | 15 lines

!bugfix: spin-spin coupling using g95 with optimization

!runs crashed in the FC contribution
!reason: mystery
!more specifically due to invalid access to WORK array segment
!which looked correctly "allocated" but g95 optimized it into a crash

!changing a WORK "allocation" into a dynamic allocation solved
!the issue

!restores tests: prop_spinspin prop_spinspin2 prop_spinspin4
!prop_newtrasoppa
!prop_newtrasoppacc dft_properties_sym dft_properties_nosym
!vibave_spinspin
!hfreqfromg prop_spinspin5 prop_newtramcscf shield_spin
!=========================================================

#endif

C
C      Used from common blocks:
C       /ENERGY/: HESREL(*,*)
C       /NUCLEI/: NUCDEP
C       /INFTAP/: LUGDR,LURDR
C       /INFVAR/: NCONF,NWOPT,NVAR
C       /INFDIM/: NVARMA
C
#include <abainf.h>
#include <spnout.h>
#include <cbirlx.h>
C
#include <energy.h>
#include <moldip.h>
#include <dipole.h>
#include <suscpt.h>
#include <spinro.h>
#include <molgfa.h>
#include <sigma.h>
#include <difsec.h>
C
#include <nuclei.h>
#include <dorps.h>
#include <inftap.h>
#include <infvar.h>
#include <infdim.h>
#include <inforb.h>
#include <inflin.h>
#include <gdvec.h>
#include <dftcom.h>
C
      IF (SKIP) RETURN
      CALL QENTER('RELAX')
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPRINT .GT. 0 .OR. SYMTST) THEN
         CALL TIMER('START ',TIMEIN,TIMOUT)
         WRITE (LUPRI,'(A,/)')
     *    '  <<<<<<<<<< Output from RELAX >>>>>>>>>> '
      END IF
C
      DOREAL = (MOLHES .OR. DIPDER .OR. VCD .OR. POLAR .OR. QPGRAD)
     &          .AND. (.NOT. DOTRPL)
      DOIMAG = (SHIELD .OR. MAGSUS .OR. VCD .OR. SPNSPN .OR. MOLGFA .OR.
     &          SPINRO) .AND. (.NOT. DOTRPL)
C
C     **********************
C     ***** Open units *****
C     **********************
C
C     Open files for MCSCF right-hand sides and solutions
C
      IF (DOREAL) THEN
         CALL GPOPEN(LUGDR,ABAGDR,'OLD','DIRECT',' ',IRAT*NVARMA,OLDDX)
         CALL GPOPEN(LURDR,ABARDR,'OLD','DIRECT',' ',IRAT*NVARMA,OLDDX)
      END IF
      IF (DOIMAG) THEN
         CALL GPOPEN(LUGDI,ABAGDI,'OLD','DIRECT',' ',IRAT*NVARMA,OLDDX)
         CALL GPOPEN(LURDI,ABARDI,'OLD','DIRECT',' ',IRAT*NVARMA,OLDDX)
      END IF
      IF (DOTRPL) THEN
         CALL GPOPEN(LUGDT,ABAGDT,'OLD','DIRECT',' ',IRAT*NVARMA,OLDDX)
         CALL GPOPEN(LURDT,ABARDT,'OLD','DIRECT',' ',IRAT*NVARMA,OLDDX)
      END IF
C
C     Allocate extra work space if contribution from triplet
C     operators is to be calculated
C
      KTEST = 1
      KLAST = KTEST + 1
      WORK(KTEST) = -999.9D0
      IF (DOTRPL) THEN
         KWRKSD = KLAST
         KWRKFC = KWRKSD + 9*9*NUCDEP*NUCDEP
         KWRKSF = KWRKFC + NUCDEP*NUCDEP
         KLAST  = KWRKSF + 9*NUCDEP*NUCDEP
         LWRK1  = LWORK - KLAST + 1
      ELSE
         IF (MOLHES) THEN
            KHESRL = KLAST
#if 0
            IF (DFTRUN) THEN 
               KDFTSH = KHESRL + MXCOOR*MXCOOR
               KLAST  = KDFTSH + MXCOOR*MXCOOR
            ELSE 
               KDFTSH = KTEST
               KLAST  = KHESRL + MXCOOR*MXCOOR
            END IF
#else
            KDFTSH = KTEST
            KLAST  = KHESRL + MXCOOR*MXCOOR
#endif
         ELSE
            KHESRL = KTEST
            KDFTSH = KTEST
         END IF
      END IF
      LWRK1  = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ADDREL',' ',KLAST,LWORK)
C
C     *****************************************
C     ***** Calculate relaxation energies *****
C     ***** and orbital rotation matrix   *****
C     *****************************************
C
      IF (DOTRPL) THEN
         CALL HESZER(SPNSD)
         CALL HESZER(SPNFC)
         CALL HESZER(SPSDFC)
         CALL DZERO(WORK(KWRKSD),9*9*NUCDEP*NUCDEP)
#ifdef VAR_G95
!bugfix for the G95 compiler May 2011
         allocate(fc_contribution(nucdep, nucdep))
         fc_contribution = 0.0d0
#else
         CALL DZERO(WORK(KWRKFC),NUCDEP*NUCDEP)
#endif
         CALL DZERO(WORK(KWRKSF),9*NUCDEP*NUCDEP)
      ELSE
         IF (POLAR)  CALL DZERO(POLARS,9)
         IF (MOLHES) CALL HESZER(WORK(KHESRL))
         IF (DIPDER) CALL DZERO(DIP1,9*NUCDEP)
         IF (SHIELD) CALL DZERO(SIGMAR,9*NUCDEP)
         IF (SPINRO) CALL DZERO(ELSPRP,9*NUCDEP)
         IF (SPNSPN) CALL HESZER(SPNPSO)
         IF (QPGRAD) CALL DZERO(SEC1,27*NUCDEP)
      END IF
      
      DO 200 ITYPE = 1, 2
         IF (DOTRPL .AND. ITYPE .EQ. 1) THEN
            GOTO 200
         ELSE IF (DOTRPL .AND. ITYPE .EQ. 2) THEN
            LUGD = LUGDT
            LURD = LURDT
            DOSEL = DOSELL
         ELSE
            IF (ITYPE .EQ. 1) THEN
               IF (.NOT.DOREAL) GO TO 200
               LUGD = LUGDR
               LURD = LURDR
               DOSEL = DOSELL
            ELSE
               IF (.NOT.DOIMAG) GO TO 200
               LUGD = LUGDI
               LURD = LURDI
               DOSEL = DOSELL
            END IF
         END IF
         DO 210 ISYM = 1,NSYM
         IF (DOSYM(ISYM)) THEN
            IF (DOTRPL) THEN
               NVEC = NTRVEC(ISYM)
            ELSE
               NVEC = NGDVEC(ISYM,ITYPE)
            END IF
            IF (NVEC .GT. 0) THEN
               CALL ABAVAR(ISYM,DOTRPL,IPRINT,WORK(KLAST),LWRK1)
               IF (IPRINT .GT. 10) THEN
                  WRITE (LUPRI,'(1X,A,I10)') ' NCONST ', NCONST
                  WRITE (LUPRI,'(1X,A,I10)') ' NVARPT ', NVARPT
               END IF
            IF (NVARPT .EQ. 0) GO TO 210
               KREL = KLAST
               KWRK = KREL  + NVEC*NVEC
               LWRK = LWRK1 - KWRK + 1
               CALL DZERO(WORK(KREL),NVEC*NVEC)
               IF (ITYPE .EQ. 1) THEN
                  IF (MOLHES .AND. POLAR) THEN
                     SELLER = .TRUE.
                  ELSE
                     IF (DIPDER .OR. QPGRAD) THEN
                        SELLER = .FALSE.
                     ELSE
                        SELLER = .TRUE.
                     END IF
                  END IF
               ELSE
                  IF (DOTRPL) THEN
                     SELLER = .TRUE.
                     IF (ANISON .OR. NUCSPI .GT. 0) DOSEL = .FALSE.
                  ELSE
                     IF ((SHIELD .OR. SPINRO) .AND. .NOT.
     &                   (SPNSPN .AND. DOPSO)) THEN
                        SELLER = .FALSE.
                     ELSE
                        IF (NUCSPI .GT. 0) DOSEL = .FALSE.
                        SELLER = .TRUE.
                     END IF
                  END IF
               END IF
               IF (SELLER) THEN
                  CALL ABASEL(WORK(KREL),WORK(KWRK),LWRK,NVEC,DOSEL,
     &                        SYMTST,LUGD,LURD,ITYPE,
     &                        DOTRPL,IPRINT)
               ELSE
                  CALL PSOB(WORK(KREL),WORK(KWRK),LWRK,NVEC,LUGD,
     &                      LURD,ITYPE,IPRINT)
               END IF
               IF (DOTRPL) THEN
#ifdef VAR_G95
!bugfix for the G95 compiler May 2011
                  CALL ADDTRP(WORK(KREL),NVEC,IPRINT,WORK(KWRKSD),
     &                        fc_contribution,WORK(KWRKSF))
#else
                  CALL ADDTRP(WORK(KREL),NVEC,IPRINT,WORK(KWRKSD),
     &                        WORK(KWRKFC),WORK(KWRKSF))
#endif
               ELSE
                  CALL ADDREL(WORK(KREL),WORK(KHESRL),SPNPSO,NVEC,
     &                        ITYPE,IPRINT)
               END IF
            END IF
         END IF
 210     CONTINUE
 200  CONTINUE
C
C     ******************************************************
C     Calculate static contribution to dft molecular hessian
C     ******************************************************
C
      IF (MOLHES.AND.DFTRUN) 
     &           CALL DFTMOLHES(LURD,WORK(KLAST),LWRK1,IPRINT)
C
C     ******************************************************
C     Calculate static contribution to dft magnetic hessian
C     ******************************************************
C
      IF (DFTRUN.AND.(MAGSUS.OR.MOLGFA).AND.(.NOT.NOLOND))
     &           CALL DFTBHES(SUSDFT,WORK(KLAST),LWRK1,IPRINT)
C
C     ***********************
C     ***** Close units *****
C     ***********************
C
      IF (DOREAL) THEN
         IF (DOWALK .OR. VCD .OR. NACME) THEN
            CALL GPCLOSE(LUGDR,'KEEP')
            CALL GPCLOSE(LURDR,'KEEP')
         ELSE
            CALL GPCLOSE(LUGDR,'DELETE')
            CALL GPCLOSE(LURDR,'DELETE')
         END IF
      END IF
      IF (DOIMAG) THEN
         CALL GPCLOSE(LUGDI,'DELETE')
         IF (VCD) THEN
            CALL GPCLOSE(LURDI,'KEEP')
         ELSE
            CALL GPCLOSE(LURDI,'DELETE')
         END IF
      END IF
      IF (DOTRPL) THEN
         CALL GPCLOSE(LUGDT,'KEEP')
         CALL GPCLOSE(LURDT,'KEEP')
      END IF
C
C     ******************************************
C     ***** Add up SD and FC contributions *****
C     ******************************************
C
      IF (DOTRPL) THEN
#ifdef VAR_G95
!bugfix for the G95 compiler May 2011
          CALL TRPMOV(SPNSD,SPNFC,SPSDFC,WORK(KWRKSD),
     &                fc_contribution,WORK(KWRKSF))
         deallocate(fc_contribution)
#else
         CALL TRPMOV(SPNSD,SPNFC,SPSDFC,WORK(KWRKSD),
     &               WORK(KWRKFC),WORK(KWRKSF))
#endif
      END IF
      IF (MOLHES .AND. .NOT. DOTRPL) CALL ADDHES(WORK(KHESRL))
C
C     ************************************
C     ***** Print relaxation Hessian *****
C     ************************************
C
      IF (IPRINT .GT. 0) THEN
         KCSTRA = KLAST
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
         IF (KLAST .GT. LWRK1) 
     &        CALL STOPIT('RELAX ','TRACOR',KLAST,LWRK1)
         IF (MOLHES .AND. .NOT. DOTRPL) THEN
            CALL HEADER('Relaxation part of Hessian',-1)
            CALL PRIHES(WORK(KHESRL),'CENTERS',
     &                  WORK(KCSTRA),WORK(KSCTRA))
         END IF
         IF (DIPDER .AND. .NOT. DOTRPL) THEN
            CALL HEADER('Relaxation part of dipole gradient',-1)
            CALL FCPRI(DDIPR,'APT',WORK(KCSTRA),WORK(KSCTRA))
         END IF
         IF (QPGRAD .AND. .NOT. DOTRPL) THEN
            CALL HEADER('Relaxation part of second moment gradient',-1)
            CALL PRISEC(DSECR,'SECDER',WORK(KCSTRA),WORK(KSCTRA))
         END IF
         IF (POLAR .AND. .NOT. DOTRPL) THEN
            CALL HEADER('Polarizabilities',-1)
            CALL POLPRI(POLARS,'AU',1)
         END IF
         IF (MAGSUS .AND. .NOT. DOTRPL) THEN
            CALL HEADER('Relaxation part of susceptibilities',-1)
            CALL OUTPUT(SUSREL,1,3,1,3,3,3,1,LUPRI)
         END IF
         IF (MOLGFA .AND. .NOT. DOTRPL) THEN
            CALL HEADER('Electronic contribution to molecular'//
     &                  ' g-factor',-1)
            CALL OUTPUT(ELMGF,1,3,1,3,3,3,1,LUPRI)
         END IF
         IF (SHIELD .AND. .NOT. DOTRPL) THEN
            CALL HEADER('Relaxation part of shielding',-1)
            CALL FCPRI(SIGMAR,'SIGMA',WORK(KCSTRA),WORK(KSCTRA))
         END IF
         IF (SPINRO .AND. .NOT. DOTRPL) THEN
            CALL HEADER('Relaxation contribution to spin-rotation'//
     &                  ' constants',-1)
            CALL FCPRI(ELSPRP,'SIGMANO',WORK(KCSTRA),WORK(KSCTRA))
         END IF
         IF (SPNSPN .AND. .NOT. DOTRPL) THEN
            CALL HEADER('PSO part of spin-spin coupling tensors',-1)
            CALL PRIHES(SPNPSO,'SPNSPN',WORK(KCSTRA),WORK(KSCTRA))
         END IF
         IF (DOTRPL) THEN
            CALL HEADER('SD(+FC) part of spin-spin coupling tensors',-1)
            CALL PRIHES(SPNSD,'SPNSPN',WORK(KCSTRA),WORK(KSCTRA))
            CALL HEADER('FC part of spin-spin coupling tensors',-1)
            CALL PRIHES(SPNFC,'SPNSPN',WORK(KCSTRA),WORK(KSCTRA))
            CALL HEADER('SD-FC coupling part of spin-spin tensor',-1)
            CALL PRIHES(SPSDFC,'SPNSPN',WORK(KCSTRA),WORK(KSCTRA))
         END IF
      END IF
C
      IF (IPRINT .GT. 0 .OR. SYMTST) CALL TIMER ('RELAX ',TIMEIN,TIMOUT)
      IF (WORK(KTEST) .NE. -999.9D0) THEN
         CALL QUIT('WORK corrupted, WORK(KTEST) is wrong')
      END IF
C
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/A/A)')
     &      ' Program stopped after RELAX as required.',
     &      ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (in RELAX) *****')
      END IF
      CALL QEXIT('RELAX')
      RETURN
      END
C  /* Deck abasel */
      SUBROUTINE ABASEL(REL,WORK,LWORK,NOP,DOSEL,SYMTST,LUGD,LURD,
     &                  ITYPE,DOTRPL,IPRINT)
C
C     Calculate second-order properties.
C     Quadratic accuracy is obtained using Harrell Sellers' formula
C
C     REL(I,J) = GD(I)*SOL(J) + SOL(I)*RES(J)
C
C     GD(I)       : Right-hand side
C     SOL(I)      : Solution
C     RES(I)      : Residual
C     REL(I,J)    : Second-order property
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.D0, DP5 = 0.5D0)
C
      LOGICAL DOSEL, SYMTST, DOTRPL
      DIMENSION REL(NOP,NOP), WORK(LWORK)
C
#include <abainf.h>
#include <gdvec.h>
#include <inflin.h>
#include <spnout.h>
C
#include <chrnos.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from ABASEL','*',103)
      KIRHS = 1
      KISOL = KIRHS + NVARPT
      KJVEC = KISOL + NVARPT
      LAST  = KJVEC + NVARPT
      IF (LAST .GT. LWORK) CALL STOPIT('ABASEL',' ',LAST,LWORK)
C
      DO 100 IOP = 1,NOP
         IF (DOTRPL) THEN
            IREC = ITRREC(IOP,LSYMPT)
         ELSE
            IREC = IGDREC(IOP,LSYMPT,ITYPE)
         END IF
         CALL READDX(LUGD,IREC,IRAT*NVARPT,WORK(KIRHS))
         IF (DOSEL) THEN
            JST = 1
            IF (DOTRPL) THEN
               IREC = 2*ITRREC(IOP,LSYMPT) - 1
            ELSE
               IREC = 2*IGDREC(IOP,LSYMPT,ITYPE) - 1
            ENDIF
            CALL READDX(LURD,IREC,IRAT*NVARPT,WORK(KISOL))
         ELSE
            IF (SYMTST) THEN
              JST = 1
            ELSE
               IF (NUCSPI .GT. 0) THEN
                  JST = 1
               ELSE
                  JST = IOP
               END IF
            END IF
         END IF
         DO 200 JOP = JST,NOP
            IF (ITYPE .EQ. 1) THEN
               JREC = 2*IGDREC(JOP,LSYMPT,ITYPE) - 1
            ELSE
               IF (DOTRPL) THEN
                  JREC = 2*ITRREC(JOP,LSYMPT) - 1
                  JCOR = ITRCOR(JOP,LSYMPT)
               ELSE
                  JREC = 2*IGDREC(JOP,LSYMPT,ITYPE) - 1
               END IF
            END IF
            IF (DOTRPL) THEN
               IF (IDORCT(ITRREC(JOP,LSYMPT)).NE.0) GO TO 200
               IF (ANISON .AND. (JCOR .GT. 0)) GO TO 200
            ELSE IF (ITYPE .EQ. 2) THEN
               IF (IDORCI(IGDREC(JOP,LSYMPT,2),2).NE.0) GO TO 200
            END IF
            CALL READDX(LURD,JREC,IRAT*NVARPT,WORK(KJVEC))
            REL(IOP,JOP) = DDOT(NVARPT,WORK(KIRHS),1,WORK(KJVEC),1)
CHJ Sep 2001: calc. -<0| d2/dx2 |0> as <d0/dQ | d0/dQ>
            IF (DOD2DQ2 .AND. IOP.EQ.1 .AND. ITYPE.EQ.1) THEN
C           ... the JOP loop and not the IOP loop has solution vectors,
C               but we only want this once, thus the IOP.eq.1 test
C           ... MOLHES is for ITYPE.eq.1
               JCOOR = IGDCOR(JOP,LSYMPT,ITYPE)
               IF (JCOOR .GT. 0) THEN
C              ... this is a nuclear coordinate
                  D2DQ2 = DDOT(NVARPT,WORK(KJVEC),1,WORK(KJVEC),1)
                  WRITE (LUPRI,'(A,I5,A,1P,D15.7)')
     &              ' <d Psi/dQ| d Psi/dQ>, Q =',JCOOR,' :',D2DQ2
            END IF
            END IF
CHJ Sep 2001 end

            IF (IPRINT.GT.10) THEN
               WRITE (LUPRI,'(A,I5,A,I5,A,D15.7)')
     &          ' REL(',IOP,',',JOP,') (lin. error): ',REL(IOP,JOP)
            END IF
            IF (DOSEL) THEN
               IF (IOP .NE. JOP) THEN
                  IF (DOTRPL) THEN
                     JREC = 2*ITRREC(JOP,LSYMPT)
                  ELSE                  
                     JREC = 2*IGDREC(JOP,LSYMPT,ITYPE)
                  ENDIF
                  CALL READDX(LURD,JREC,IRAT*NVARPT,WORK(KJVEC))
                  REL(IOP,JOP) = REL(IOP,JOP)
     &               + DDOT(NVARPT,WORK(KISOL),1,WORK(KJVEC),1)
                  IF (IPRINT.GT.10) THEN
                     WRITE (LUPRI,'(A,I5,A,I5,A,D15.7)')
     &           ' REL(',IOP,',',JOP,') (Seller corr):',REL(IOP,JOP)
                  END IF
               END IF
            ELSE IF (.NOT.SYMTST) THEN
               REL(JOP,IOP) = REL(IOP,JOP)
            END IF
 200     CONTINUE
 100  CONTINUE
      IF (SYMTST .OR. IPRINT .GT. 5) THEN
         CALL HEADER('Relaxation matrix for sym. '//CHRNOS(LSYMPT),-1)
         CALL OUTPUT(REL,1,NOP,1,NOP,NOP,NOP,1,LUPRI)
      END IF
      IF (DOSEL .OR. SYMTST) THEN
         DIFMAX = D0
         DO 300 I = 1, NOP
            DO 400 J = 1, I - 1
               DIF = ABS(REL(I,J) - REL(J,I))
               IF (DIF .GT. DIFMAX) THEN
                  DIFMAX = DIF
                  IMX = I
                  JMX = J
               END IF
               AVER = DP5*(REL(I,J) + REL(J,I))
               REL(I,J) = AVER
               REL(J,I) = AVER
  400       CONTINUE
  300    CONTINUE
         IF (SYMTST .OR. (IPRINT .GT. 5)) THEN
            CALL HEADER('Relaxation matrix after averaging',-1)
            CALL OUTPUT(REL,1,NOP,1,NOP,NOP,NOP,1,LUPRI)
            WRITE (LUPRI,'(A,1P,D16.8)') ' Largest deviation:', DIFMAX
            WRITE (LUPRI,'(A,2I5)')      ' Element:          ', IMX,JMX
         END IF
      END IF
      RETURN
      END
C  /* Deck psob */
      SUBROUTINE PSOB(REL,WORK,LWORK,NOP,LUGD,LURD,ITYPE,IPRINT)
C
C     Calculate second-order properties for chemical shielding.
C
C     REL(I,J) = GDPSO(I)*SOLB(J)
C
C     GDPSO(I)       : Right-hand side for PSO
C     SOLB(I)        : Solution vector for magnetic field B
C     REL(I,J)       : Second-order property
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.D0, DP5 = 0.5D0)
C
      DIMENSION REL(NOP,NOP), WORK(LWORK)
C
#include <gdvec.h>
#include <abainf.h>
#include <spnout.h>
#include <inflin.h>
C
#include <chrnos.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from PSOB','*',103)
      KIRHS = 1
      KJVEC = KIRHS + NVARPT
      LAST  = KJVEC + NVARPT
      IF (LAST .GT. LWORK) CALL STOPIT('PSOB',' ',LAST,LWORK)
C
      DO 100 IOP = 1,NOP
         IREC  = IGDREC(IOP,LSYMPT,ITYPE)
         ICOOR = IGDCOR(IOP,LSYMPT,ITYPE)
         CALL READDX(LUGD,IREC,IRAT*NVARPT,WORK(KIRHS))
         DO 200 JOP = 1, NOP
            JREC  = IGDREC(JOP,LSYMPT,ITYPE)
            JCOOR = IGDCOR(JOP,LSYMPT,ITYPE)
            IF (ITYPE .EQ. 1) THEN
               IF (.NOT. ((MOLHES .AND. JCOOR .GT. 0) .OR. 
     &             (POLAR .AND. JCOOR .LT. 0) .OR. 
     &             (QPGRAD .AND. JCOOR .LT. -10))) GOTO 200
               JREC = 2*JREC - 1
            ELSE
               IF ((JCOOR .LT. 0 .AND. (SPNSPN .AND. DOPSO .AND. 
     &             .NOT. MAGSUS)) .OR. (JCOOR .GT. 0 .AND. .NOT. 
     &              (SPNSPN .AND. DOPSO))) GOTO 200
               JREC = 2*JREC - 1
            END IF
            IF (ITYPE .EQ. 2 .AND.
     &          IDORCI(IGDREC(JOP,LSYMPT,ITYPE),ITYPE).NE.0) GO TO 200
            CALL READDX(LURD,JREC,IRAT*NVARPT,WORK(KJVEC))
            REL(IOP,JOP) = DDOT(NVARPT,WORK(KIRHS),1,WORK(KJVEC),1)
            REL(JOP,IOP) = REL(IOP,JOP)
            IF (IPRINT.GT.30) THEN
               WRITE (LUPRI,'(A,I5,A,I5,A,D15.7)')
     &           ' REL(',IOP,',',JOP,'):              ',REL(IOP,JOP)
            END IF
 200     CONTINUE
 100  CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Relaxation matrix for sym. '//CHRNOS(LSYMPT),-1)
         CALL OUTPUT(REL,1,NOP,1,NOP,NOP,NOP,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck addrel */
      SUBROUTINE ADDREL(REL,HESREL,SPNPSO,NOP,ITYPE,IPRINT)
C
C     110190 tuh
C
C     This subroutine retrieves relaxation terms from REL and
C     adds these to the matrices HESREL, DDIPR, and POLARS
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0)
C
      DIMENSION REL(NOP,NOP), HESREL(MXCOOR,MXCOOR), 
     &          SPNPSO(MXCOOR,MXCOOR)
#include <abainf.h>
#include <spnout.h>
#include <energy.h>
#include <sigma.h>
#include <suscpt.h>
#include <molgfa.h>
#include <spinro.h>
#include <dipole.h>
#include <moldip.h>
#include <inflin.h>
#include <difsec.h>
#include <symmet.h>
#include <gdvec.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from ADDREL','*',103)
      DO 100 I = 1, NOP
         ICOOR = IGDCOR(I,LSYMPT,ITYPE)
         DO 200 J = 1, NOP
            JCOOR = IGDCOR(J,LSYMPT,ITYPE)
            RELIJ = REL(I,J)
            IF (IPRINT .GT. 5) WRITE (LUPRI,'(1X,A,2I5,F12.6)')
     &         ' Coordinates and element: ',ICOOR, JCOOR, RELIJ
C
C           Molecular Hessian
C
            IF ((ICOOR .GT. 0) .AND. (JCOOR .GT. 0)) THEN
               IF (ITYPE .EQ. 1 .AND. MOLHES) THEN
                  HESREL(ICOOR,JCOOR) = RELIJ
                  HESREL(JCOOR,ICOOR) = RELIJ
               ELSE
                  IF (SPNSPN .AND. DOPSO) THEN
                     SPNPSO(ICOOR,JCOOR) = -RELIJ
                     SPNPSO(JCOOR,ICOOR) = -RELIJ
                  END IF
               END IF
C
C           Dipole gradient
C
            ELSE IF ((ICOOR .GT. 0) .AND. (JCOOR .LT. 0)) THEN
               IF (ITYPE .EQ. 1 .AND. (JCOOR .GT. -4)) THEN
                  DDIPR (-JCOOR,ICOOR) = RELIJ
               ELSE IF (ITYPE .EQ. 2 .AND. (JCOOR .GT. -4)) THEN
                  SIGMAR(-JCOOR,ICOOR) = - RELIJ
                  ELSPRP(-JCOOR,ICOOR) = - RELIJ
C
C           Second moment gradient
C
               ELSE IF (ITYPE .EQ. 1 .AND. (JCOOR .LT. -10)) THEN
                  IX = - JCOOR / 10
                  IY = - MOD(JCOOR,10)
                  DSECR(IPTAX(IX,1),IPTAX(IY,1),ICOOR) = RELIJ
                  DSECR(IPTAX(IY,1),IPTAX(IX,1),ICOOR) = RELIJ
               END IF
C
C           Dipole gradient
C
            ELSE IF ((ICOOR .LT. 0) .AND. (JCOOR .GT. 0)) THEN
               IF (ITYPE .EQ. 1 .AND. (ICOOR .GT. -4)) THEN
                  DDIPR (-ICOOR,JCOOR) = RELIJ
               ELSE IF (ITYPE .EQ. 2 .AND. (ICOOR .GT. -4)) THEN
                  SIGMAR(-ICOOR,JCOOR) = - RELIJ
                  ELSPRP(-ICOOR,JCOOR) = - RELIJ
C
C           Second moment gradient
C
               ELSE IF (ITYPE .EQ. 1 .AND. (ICOOR .LT. -10)) THEN
                  IX = - ICOOR / 10
                  IY = - MOD(ICOOR,10)
                  DSECR(IPTAX(IX,1),IPTAX(IY,1),JCOOR) = RELIJ
                  DSECR(IPTAX(IY,1),IPTAX(IX,1),JCOOR) = RELIJ
               END IF
C
C           Polarizabilities and susceptibilities
C
            ELSE IF ((ICOOR .LT. 0) .AND. (JCOOR .LT. 0)) THEN
               IF (ITYPE .EQ. 1 .AND. ((ICOOR .GT. -4) .AND.
     &            (JCOOR .GT. -4))) THEN
                  POLARS(-ICOOR,-JCOOR) = - RELIJ
                  POLARS(-JCOOR,-ICOOR) = - RELIJ
               ELSE
                  IF ((ICOOR .GT. -4) .AND. (JCOOR .GT. -4)) THEN
                     SUSREL(-ICOOR,-JCOOR) = - RELIJ
                     SUSREL(-JCOOR,-ICOOR) = - RELIJ
                     ELMGF (-ICOOR,-JCOOR) = - RELIJ
                     ELMGF (-JCOOR,-ICOOR) = - RELIJ
                  END IF
               END IF
            END IF
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck addtrp */
      SUBROUTINE ADDTRP(REL,NOP,IPRINT,WRKSD,WRKFC,WRKSDF)
C
C     KR, Oct.-92, based on ADDREL, but modified for triplet
C     contributions to spin-spin-interactions
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
      DIMENSION REL(NOP,NOP), WRKSD(9*NUCDEP,9*NUCDEP),
     &          WRKFC(NUCDEP,NUCDEP), WRKSDF(9*NUCDEP,NUCDEP)
#include <abainf.h>
#include <energy.h>
#include <inflin.h>
#include <gdvec.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from ADDTRP','*',103)
      DO 100 I = 1, NOP
         ICOOR = ITRCOR(I,LSYMPT)
         DO 200 J = 1, I
            JCOOR = ITRCOR(J,LSYMPT)
            RELIJ = REL(I,J)
            IF (IPRINT .GT. 5) WRITE (LUPRI,'(1X,A,2I5,F12.6)')
     &           ' Coordinates and element: ',ICOOR,JCOOR,RELIJ
            IF ((ICOOR .GT. 0) .AND. (JCOOR .GT. 0)) THEN
               WRKSD(ICOOR,JCOOR) = RELIJ
               WRKSD(JCOOR,ICOOR) = RELIJ
            ELSE IF ((ICOOR .LT. 0) .AND. (JCOOR .LT. 0)) THEN
               WRKFC(-ICOOR,-JCOOR) = RELIJ
               WRKFC(-JCOOR,-ICOOR) = RELIJ
            ELSE IF ((ICOOR .LT. 0) .AND. (JCOOR .GT. 0)) THEN
               WRKSDF(JCOOR,-ICOOR) = RELIJ
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck trpmov */
      SUBROUTINE TRPMOV(SPNSD,SPNFC,SPSDFC,WRKSD,WRKFC,WRKSDF)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <nuclei.h>
      DIMENSION WRKSD(9*NUCDEP,9*NUCDEP), WRKFC(NUCDEP,NUCDEP),
     &          WRKSDF(9*NUCDEP,NUCDEP), SPNSD(MXCOOR,MXCOOR),
     &          SPNFC(MXCOOR,MXCOOR), SPSDFC(MXCOOR,MXCOOR)
#include <cbitrp.h>
#include <spnout.h>
#include <symmet.h>
#include <ibtfun.h>
C
      CALL DZERO(SPNSD,MXCOOR*MXCOOR)
      IF (DOSD .OR. DOSDFC) THEN
         JATOM1 = 0
         DO IREP = 0, MAXREP
         DO IATOM1 = 1, NUCIND
            IF (IBTAND(IREP,ISTBNU(IATOM1)).EQ.0) THEN
               JATOM1 = JATOM1 + 1
               JATOM2 = 0
               DO  IREP1 = 0, MAXREP
               DO  IATOM2 = 1, NUCIND
                  IF (IBTAND(IREP1,ISTBNU(IATOM2)).EQ.0) THEN
                     JATOM2 = JATOM2 + 1
                     DO ICOOR = 1, 3
                        IREP2 = IBTXOR(ISYMAX(ICOOR,2),IREP)
                        ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR,IREP2,2)
                        DO JCOOR = 1, 3
                           IREP3 = IBTXOR(ISYMAX(JCOOR,2),IREP1)
                           ISCOR2 = IPTCNT(3*(IATOM2-1)+JCOOR,IREP3,2)
                           IF (ISCOR1 .GT. 0 .AND. ISCOR2 .GT. 0) THEN
                              DO KCOOR = 1, 3
                                 IF (KCOOR .GT. ICOOR) THEN
                                   IREP4 = IBTXOR(ISYMAX(KCOOR,2),IREP)
                                   ISCOR3 = IPTCNT(3*(IATOM1-1)+KCOOR,
     &                                             IREP4,2)
                                   IADR1 = 3*(ISCOR3 - 1) + ICOOR
                                 ELSE
                                    IADR1 = 3*(ISCOR1 - 1) + KCOOR
                                 END IF
                                 IF (KCOOR .GT. JCOOR) THEN
                                   IREP4 = IBTXOR(ISYMAX(KCOOR,2),IREP1)
                                   ISCOR3 = IPTCNT(3*(IATOM2-1)+KCOOR,
     &                                             IREP4,2)
                                   IADR2 = 3*(ISCOR3 - 1) + JCOOR
                                 ELSE
                                    IADR2 = 3*(ISCOR2 - 1) + KCOOR
                                 END IF
                             SPNSD(ISCOR1,ISCOR2) = SPNSD(ISCOR1,ISCOR2)
     &                                            + WRKSD(IADR1,IADR2)
                              END DO
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END DO
            END IF
         END DO
         END DO
      END IF
      IF (DOFC) THEN
         JATOM1 = 0
         DO 400 IREP = 0, MAXREP
            KATOM1 = JATOM1
            DO 500 IATOM1 = 1, NUCIND
            IF (IBTAND(IREP,ISTBNU(IATOM1)).EQ.0) THEN
               JATOM1 = JATOM1 + 1
               JATOM2 = KATOM1
               DO 510 IATOM2 = 1, NUCIND
               IF (IBTAND(IREP,ISTBNU(IATOM2)).EQ.0) THEN
                  JATOM2 = JATOM2 + 1
                  DO 520 ICOOR = 1, 3
                     IREP2 = IBTXOR(ISYMAX(ICOOR,2),IREP)
                     ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR,IREP2,2)
                     ISCOR2 = IPTCNT(3*(IATOM2 - 1) + ICOOR,IREP2,2)
                     IF (ISCOR1 .GT. 0 .AND. ISCOR2 .GT. 0) THEN
                        SPNFC(ISCOR1,ISCOR2) = WRKFC(JATOM1,JATOM2)
                     END IF
  520             CONTINUE
               END IF
  510          CONTINUE
            END IF
  500       CONTINUE
  400    CONTINUE
      END IF
      IF (DOFC .AND. DOSD) THEN
         JATOM1 = 0
         DO 600 IREP = 0, MAXREP
         DO 610 IATOM1 = 1, NUCIND
            IF (IBTAND(IREP,ISTBNU(IATOM1)).EQ.0) THEN
               JATOM1 = JATOM1 + 1
               JATOM2 = 0
               DO 615 IREP1 = 0, MAXREP
               DO 620 IATOM2 = 1, NUCIND
                  IF (IBTAND(IREP1,ISTBNU(IATOM2)).EQ.0) THEN
                     JATOM2 = JATOM2 + 1
                     DO 630 ICOOR = 1, 3
                        IREP2 = IBTXOR(ISYMAX(ICOOR,2),IREP)
                        ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR,IREP2,2)
                        DO 640 JCOOR = 1, 3
                           IREP3 = IBTXOR(ISYMAX(JCOOR,2),IREP1)
                           ISCOR2 = IPTCNT(3*(IATOM2-1)+JCOOR,IREP3,2)
                           IF (ISCOR1 .GT. 0 .AND. ISCOR2 .GT. 0) THEN
                              IF (JCOOR .GT. ICOOR) THEN
                                 IREP4 = IBTXOR(ISYMAX(JCOOR,2),IREP)
                                 ISCOR3 = IPTCNT(3*(IATOM1-1)+JCOOR,
     &                                             IREP4,2)
                                 IADR1 = 3*(ISCOR3 - 1) + ICOOR
                                 IREP4 = IBTXOR(ISYMAX(JCOOR,2),IREP1)
                                 ISCOR3 = IPTCNT(3*(IATOM2-1)+JCOOR,
     &                                           IREP4,2)
                                 IADR2 = 3*(ISCOR3 - 1) + ICOOR
                              ELSE
                                 IADR1 = 3*(ISCOR1 - 1) + JCOOR
                                 IREP4 = IBTXOR(ISYMAX(ICOOR,2),IREP1)
                                 ISCOR3 = IPTCNT(3*(IATOM2-1)+ICOOR,
     &                                           IREP4,2)
                                 IADR2 = 3*(ISCOR3 - 1) + JCOOR
                              END IF
                            SPSDFC(ISCOR1,ISCOR2) = WRKSDF(IADR1,JATOM2)
     &                                            + WRKSDF(IADR2,JATOM1)
                           END IF
 640                    CONTINUE
 630                 CONTINUE
                  END IF
 620           CONTINUE
 615        CONTINUE
            END IF
 610     CONTINUE
 600     CONTINUE
      END IF
      RETURN
      END
