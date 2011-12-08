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
      SUBROUTINE CTOABA(WORK,LWORK)
C...
C...  This subroutine was written by Andrea Ligabue using
C...  the subroutine LNRABA as a model (1999)
C...
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <trkoor.h>
#include <sigma.h>
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <cbilnr.h>
#include <suscpt.h>
#include <infpri.h>
      LOGICAL CICLC, HFCLC, TRIPLE, EXECLC, FOUND
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL1, LABEL2, LABINT(4*MXCOOR+9)
      PARAMETER (D05=0.5D0,D025=0.25)
      LOGICAL TODOINT
C
#include <cbiexc.h>
#include <inflin.h>
#include <infvar.h>
#include <infdim.h>
#include <inforb.h>
#include <nuclei.h>
#include <inftap.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <maxmom.h>
#include <maxaqn.h>
#include <symmet.h>
#include <abainf.h>
#include <gnrinf.h>
#include <infsop.h>
#include <ibtfun.h>
C
      CALL QENTER('CTOABA')
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPRLNR .GT. 0) WRITE (LUPRI,'(A,/)')
     *    '1 <<<<<<<<<< Output from CTOABA >>>>>>>>>> '
C
      IPRRSP = IPRLNR
C
C     Get reference state
C     ===================
C
C     1. Work Allocations:
C
      IF (ABASOP) THEN
        LUDV   = NORBT * NORBT
        LPVX   = LPVMAT
      ELSE
        LUDV   = N2ASHX
        LPVX   = 0
      ENDIF
      KFREE  = 1
      LFREE  = LWORK
C      
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDV  ,LUDV  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVX  ,LPVX  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXINDX,LCINDX,WORK,KFREE,LFREE)
C
      KWORK1 = KFREE
      LWORK1 = LFREE
C
      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO),WORK(KFREE),LFREE)
      IF (.NOT.FOUND) CALL QUIT('CTOABA error: CMO not found on SIRIFC')
      IF (NASHT .GT. 0) THEN
         CALL RD_SIRIFC('DV',FOUND,WORK(KWORK1),WORK(KWORK1),LWORK1)
         IF (.NOT.FOUND)
     &      CALL QUIT('CTOABA error: DV not found on SIRIFC')
         CALL DSPTSI(NASHT,WORK(KWORK1),WORK(KUDV))
      END IF
C
      ISYM = 1
      CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK1),LWORK1)
C
      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LWORK1,0)
C
C     SOPPA :
C
      IF (ABASOP) THEN
C
C        Initialize XINDX
C
         CALL DZERO(WORK(KXINDX),LCINDX)
C
C        Find address array's for SOPPA calculation
C
         CALL SET2SOPPA(WORK(KXINDX+KABSAD-1),WORK(KXINDX+KABTAD-1),
     *                  WORK(KXINDX+KIJSAD-1),WORK(KXINDX+KIJTAD-1),
     *                  WORK(KXINDX+KIJ1AD-1),WORK(KXINDX+KIJ2AD-1),
     *                  WORK(KXINDX+KIJ3AD-1),WORK(KXINDX+KIADR1-1))
C
C
         REWIND (LUSIFC)
         IF (CCPPA) THEN
            CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
         ELSE
            CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
         ENDIF
C
C        reads the MP2 or CCSD correlation coefficients into PV
C
         CALL READT (LUSIFC,LPVMAT,WORK(KPVX))
C
         IF (IPRLNR.GT.10) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')' EXCIT1 : CCSD correlation ',
     &                           'coefficients'
            ELSE
               WRITE(LUPRI,'(/A,A)')' EXCIT1 :',
     &                              ' MP2 correlation coefficients'
            ENDIF
            CALL OUTPUT(WORK(KPVX),1,LPVMAT,1,1,LPVMAT,1,1,LUPRI)
         END IF
C
C        reads the MP2 or CCSD second order one particle density matrix 
C
         CALL READT (LUSIFC,NORBT*NORBT,WORK(KUDV))
C
C        UDV contains the MP2 one-density. Remove the diagonal
C        contribution from the zeroth order. (Added in MP2FAC)
C
         IF (IPRLNR.GT.10) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')' RSPMC : CCSD density'
            ELSE
               WRITE(LUPRI,'(/A)')' RSPMC : MP2 density'
            END IF
            CALL OUTPUT(WORK(KUDV),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,
     &                  LUPRI)
         END IF
C
         CALL SOPUDV(WORK(KUDV))
      END IF
C
C
C     Construct property-integrals and write to LUPROP
C     ================================================
C
C     2. Work Allocations:
C
      KIDSYM = KWORK1
      KIDADR = KIDSYM + 9*MXCENT
      KWORK2 = KIDADR + 9*MXCENT
      LWORK2 = LWORK  - KWORK2
C
      NLBTOT = 0
      NLBSHIS = 4
C
      NCOMP  = 0
      NPATOM = 0
cLig  <> added the TODOINT to see if the property was already in the file
      IF (TODOINT('DIPVEL  ',LUPROP)) THEN
        CALL GETLAB('DIPVEL',6,LABINT,WORK(KIDSYM),LUPROP)
      ELSE
        CALL GET1IN(DUMMY,'DIPVEL ',NCOMP,WORK(KWORK2),LWORK2,
     &              LABINT,WORK(KIDSYM),WORK(KIDADR),
     &              IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &              IPR1IN)
      ENDIF
      NLAB = 3
      CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      IF (MAGSUS) THEN
         NCOMP  = 0
         CALL GET1IN(DUMMY,'RANGMO ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &               IPR1IN)
         NLAB = 9
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
         NLBSHIS = 13
      ENDIF
cLig  <> added the TODOINT to see if the property was already in the file
      IF (SHIELD) THEN
         IF (TODOINT('RPSO    ',LUPROP)) THEN
           CALL GETLAB('RPSO',4,LABINT,WORK(KIDSYM),LUPROP)
         ELSE
            NCOMP  = 0
            CALL GET1IN(DUMMY,'RPSO   ',NCOMP,WORK(KWORK2),LWORK2,
     &                  LABINT,WORK(KIDSYM),WORK(KIDADR),
     &                  IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &                  IPR1IN)
         ENDIF
         NLAB = 3*NCOOR
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
         NLBSHIE = NLBTOT
         IF (TODOINT('PSO     ',LUPROP)) THEN
           CALL GETLAB('PSO',3,LABINT,WORK(KIDSYM),LUPROP)
         ELSE
            NCOMP  = 0
            CALL GET1IN(DUMMY,'PSO    ',NCOMP,WORK(KWORK2),LWORK2,
     &                  LABINT,WORK(KIDSYM),WORK(KIDADR),
     &                  IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &                  IPR1IN)
         ENDIF
         NLAB = NCOOR
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      ENDIF
C
C     Set variables and logicals
C
      CICLC  = .FALSE.
      HFCLC  = NASHT .LE. 1
      TRIPLE = .FALSE.
      EXECLC = .FALSE.
cLig
CLig  NABATY = 1 for real operators .. -1 for imm. op.
cLig  NABAOP is the number of right hand sides
      NABATY = -1
      NABAOP = 1
C
C     Zero the property tensors
C
      IF (MAGSUS) CALL DZERO(SUSDZD,9)
      IF (SHIELD) THEN
        CALL DZERO(SIGMADZ,9*NUCDEP)
        CALL DZERO(SIGMASFTP,9*NUCDEP)
        CALL DZERO(SIGMASFTM,9*NUCDEP)
      ENDIF
C   
C        Loop over the right operators which are the
C        the dipole velocity operators 
C        ===========================================
C   
      LUSOVE = -1
      LUGDVE = -1
      LUREVE = -1
      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      DO 300 IDIP = 1,3
         IF (IDIP.EQ.1) LABEL1 = 'XDIPVEL '
         IF (IDIP.EQ.2) LABEL1 = 'YDIPVEL '
         IF (IDIP.EQ.3) LABEL1 = 'ZDIPVEL '
         ISYM=ISYMAX(IDIP,1)+1
C         
         CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK2),LWORK2)
C
         IF (NFRVAL.GT.0) THEN
C
C           3. Work Allocations:
C
            KGD1   = KWORK1
            KWRKG1 = KGD1
            LWRKG1 = LWORK - KWRKG1
            KSLV   = KGD1 + 2*NVARPT
            KLAST  = KSLV + 2*NVARPT
            IF (KLAST.GT.LWORK) CALL STOPIT('CTOABA',' ',KLAST,LWORK)
            KWRK = KLAST
            LWRK = LWORK - KLAST + 1
C
C           Find right hand side for right operator and write to file
C           =========================================================
C
            KSYMOP = ISYM
            TRPLET = .FALSE.
C
            CALL GETGPV(LABEL1,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &           WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
            REWIND LUGDVE
            CALL WRITT(LUGDVE,2*NVARPT,WORK(KWRKG1))
            IF (IPRLNR.GT.3) THEN
               WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL1
               CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,2,1,LUPRI)
            ENDIF
C
C           Calculate eigenvector and write to file
C           =======================================
C
            CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     &            FRVAL,NFRVAL,NABATY,NABAOP,LABEL1,LUGDVE,LUSOVE,
     &            LUREVE,THCLNR,MAXITE,IPRRSP,MXRM,MXPHP,
     &            WORK(KWRK),LWRK)
C
C           Loop over the left side  property operators
C           ===========================================
C
            DO 200 IPRLBL = 4, NLBTOT
C
C              Find label and symmetry of the left side operator
C
               LABEL2 = LABAPP(IPRLBL)
               KSYM   = LABSYM(IPRLBL)

C
C              If symmetry of right operator equals symmetry of
C              the left operator, that is if ISYM = KSYM, then
C              ================================================
C              (otherwise 2. order property SNDPRP is zero)
C
               IF (KSYM.EQ.ISYM) THEN
                  KSYMOP = ISYM
                  TRPLET = .FALSE.
C
C                 Find right hand side for left operator
C                 ========================================
C
                  CALL GETGPV(LABEL2,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &                       WORK(KPVX),WORK(KXINDX),ANTSYM,
     &                       WORK(KWRKG1),LWRKG1)
C
                  IF (IPRLNR.GT.3) THEN
                     WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL2
                     CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,2,1,
     &                           LUPRI)
                  ENDIF
C
C                 Form second order properties SNDPRP
C                 ===================================
C
                  REWIND LUSOVE
                  CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
                  IF (IPRLNR.GT.3) THEN
                     WRITE (LUPRI,'(2A)') 'Solution Vector, label: ',
     &                                    LABEL1
                     CALL OUTPUT(WORK(KSLV),1,NVARPT,1,2,NVARPT,2,1,
     &                           LUPRI)
                  ENDIF
C
                  SNDPRP = DDOT(2*NVARPT,WORK(KSLV),1,WORK(KGD1),1)
C
                  IF (IPRLNR.GT.2) THEN
                     WRITE (LUPRI,'(4A,F15.8)')
     &               ' Second order property for ',LABEL2,LABEL1,
     &               ' = ',SNDPRP
                  ENDIF
C
C                 Write properties into the various property matrices
C                 ===================================================
C
C                 Magnetizability
C                 ----------------
                  IF (IPRLBL.LT.NLBSHIS) THEN 
                     IF (LABEL2(1:2).EQ.'XX') THEN  
                        IF (LABEL1(1:1).EQ.'Z') THEN  
                           SUSDZD(IPTAX(1,2),IPTAX(2,2))= 
     &                        SUSDZD(IPTAX(1,2),IPTAX(2,2))+ SNDPRP
                        ELSE IF (LABEL1(1:1).EQ.'Y') THEN 
                           SUSDZD(IPTAX(1,2),IPTAX(3,2))= 
     &                        SUSDZD(IPTAX(1,2),IPTAX(3,2))- SNDPRP
                        ENDIF
                     ELSE IF (LABEL2(1:2).EQ.'XY') THEN
                        IF (LABEL1(1:1).EQ.'Z') THEN
                           SUSDZD(IPTAX(2,2),IPTAX(2,2))=
     &                        SUSDZD(IPTAX(2,2),IPTAX(2,2))+ SNDPRP
                        ELSE IF (LABEL1(1:1).EQ.'Y') THEN
                           SUSDZD(IPTAX(2,2),IPTAX(3,2))=
     &                        SUSDZD(IPTAX(2,2),IPTAX(3,2))- SNDPRP
                        ENDIF
                     ELSE IF (LABEL2(1:2).EQ.'XZ') THEN
                        IF (LABEL1(1:1).EQ.'Z') THEN
                           SUSDZD(IPTAX(3,2),IPTAX(2,2))=
     &                        SUSDZD(IPTAX(3,2),IPTAX(2,2))+ SNDPRP
                        ELSE IF (LABEL1(1:1).EQ.'Y') THEN
                           SUSDZD(IPTAX(3,2),IPTAX(3,2))=
     &                        SUSDZD(IPTAX(3,2),IPTAX(3,2))- SNDPRP
                        ENDIF
                     ELSE IF (LABEL2(1:2).EQ.'YX') THEN
                        IF (LABEL1(1:1).EQ.'X') THEN
                           SUSDZD(IPTAX(1,2),IPTAX(3,2))=
     &                        SUSDZD(IPTAX(1,2),IPTAX(3,2))+ SNDPRP
                        ELSE IF (LABEL1(1:1).EQ.'Z') THEN
                           SUSDZD(IPTAX(1,2),IPTAX(1,2))=
     &                        SUSDZD(IPTAX(1,2),IPTAX(1,2))- SNDPRP
                        ENDIF
                     ELSE IF (LABEL2(1:2).EQ.'YY') THEN
                        IF (LABEL1(1:1).EQ.'X') THEN
                           SUSDZD(IPTAX(2,2),IPTAX(3,2))=
     &                        SUSDZD(IPTAX(2,2),IPTAX(3,2))+ SNDPRP
                        ELSE IF (LABEL1(1:1).EQ.'Z') THEN
                           SUSDZD(IPTAX(2,2),IPTAX(1,2))=
     &                        SUSDZD(IPTAX(2,2),IPTAX(1,2))- SNDPRP
                        ENDIF
                     ELSE IF (LABEL2(1:2).EQ.'YZ') THEN
                        IF (LABEL1(1:1).EQ.'X') THEN
                           SUSDZD(IPTAX(3,2),IPTAX(3,2))=
     &                        SUSDZD(IPTAX(3,2),IPTAX(3,2))+ SNDPRP
                        ELSE IF (LABEL1(1:1).EQ.'Z') THEN
                           SUSDZD(IPTAX(3,2),IPTAX(1,2))=
     &                        SUSDZD(IPTAX(3,2),IPTAX(1,2))- SNDPRP
                        ENDIF
                     ELSE IF (LABEL2(1:2).EQ.'ZX') THEN
                        IF (LABEL1(1:1).EQ.'Y') THEN
                           SUSDZD(IPTAX(1,2),IPTAX(1,2))=
     &                        SUSDZD(IPTAX(1,2),IPTAX(1,2))+ SNDPRP
                        ELSE IF (LABEL1(1:1).EQ.'X') THEN
                           SUSDZD(IPTAX(1,2),IPTAX(2,2))=
     &                        SUSDZD(IPTAX(1,2),IPTAX(2,2))- SNDPRP
                        ENDIF
                     ELSE IF (LABEL2(1:2).EQ.'ZY') THEN
                        IF (LABEL1(1:1).EQ.'Y') THEN
                           SUSDZD(IPTAX(2,2),IPTAX(1,2))=
     &                        SUSDZD(IPTAX(2,2),IPTAX(1,2))+ SNDPRP
                        ELSE IF (LABEL1(1:1).EQ.'X') THEN
                           SUSDZD(IPTAX(2,2),IPTAX(2,2))=
     &                        SUSDZD(IPTAX(2,2),IPTAX(2,2))- SNDPRP
                        ENDIF
                     ELSE IF (LABEL2(1:2).EQ.'ZZ') THEN
                        IF (LABEL1(1:1).EQ.'Y') THEN
                           SUSDZD(IPTAX(3,2),IPTAX(1,2))=
     &                        SUSDZD(IPTAX(3,2),IPTAX(1,2))+ SNDPRP
                        ELSE IF (LABEL1(1:1).EQ.'X') THEN
                           SUSDZD(IPTAX(3,2),IPTAX(2,2))=
     &                        SUSDZD(IPTAX(3,2),IPTAX(2,2))- SNDPRP
                        ENDIF
                     ENDIF
                  ENDIF
C
C                 Nuclear Shieldings
C                 ------------------
                  IF ((IPRLBL.GE.NLBSHIS).AND.
     &                (IPRLBL.LE.NLBSHIE))  THEN
                     INUCLEO = (IPRLBL-NLBSHIS)/3+1
                     IF (LABEL1.EQ.'XDIPVEL') THEN
                        IF (LABEL2(8:8).EQ.'Y') THEN
                           SIGMADZ(IPTAX(3,2),INUCLEO)=
     &                        SIGMADZ(IPTAX(3,2),INUCLEO)+SNDPRP
                        ELSE IF (LABEL2(8:8).EQ.'Z') THEN
                           SIGMADZ(IPTAX(2,2),INUCLEO)=
     &                        SIGMADZ(IPTAX(2,2),INUCLEO)-SNDPRP
                        ENDIF
                     ELSE IF (LABEL1.EQ.'YDIPVEL') THEN
                        IF (LABEL2(8:8).EQ.'Z') THEN
                           SIGMADZ(IPTAX(1,2),INUCLEO)=
     &                        SIGMADZ(IPTAX(1,2),INUCLEO)+SNDPRP
                        ELSE IF (LABEL2(8:8).EQ.'X') THEN
                           SIGMADZ(IPTAX(3,2),INUCLEO)=
     &                        SIGMADZ(IPTAX(3,2),INUCLEO)-SNDPRP
                        ENDIF
                     ELSE IF (LABEL1.EQ.'ZDIPVEL') THEN
                        IF (LABEL2(8:8).EQ.'X') THEN
                           SIGMADZ(IPTAX(2,2),INUCLEO)=
     &                        SIGMADZ(IPTAX(2,2),INUCLEO)+SNDPRP
                        ELSE IF (LABEL2(8:8).EQ.'Y') THEN
                           SIGMADZ(IPTAX(1,2),INUCLEO)=
     &                        SIGMADZ(IPTAX(1,2),INUCLEO)-SNDPRP
                        ENDIF
                     ENDIF
                  ENDIF
C
C                 Nuclear Shielding Corrections
C                 -----------------------------
                  IF (IPRLBL.GT.NLBSHIE) THEN
                     INUCLEO = (IPRLBL-NLBSHIE)
                     IF (LABEL1.EQ.'XDIPVEL') THEN
                        SIGMASFTP(IPTAX(3,2),INUCLEO)=
     &                     SIGMASFTP(IPTAX(3,2),INUCLEO)+SNDPRP
                        SIGMASFTM(IPTAX(2,2),INUCLEO)=
     &                     SIGMASFTM(IPTAX(2,2),INUCLEO)+SNDPRP
                     ELSE IF (LABEL1.EQ.'YDIPVEL') THEN
                        SIGMASFTP(IPTAX(1,2),INUCLEO)=
     &                     SIGMASFTP(IPTAX(1,2),INUCLEO)+SNDPRP
                        SIGMASFTM(IPTAX(3,2),INUCLEO)=
     &                     SIGMASFTM(IPTAX(3,2),INUCLEO)+SNDPRP
                     ELSE IF (LABEL1.EQ.'ZDIPVEL') THEN
                        SIGMASFTP(IPTAX(2,2),INUCLEO)=
     &                     SIGMASFTP(IPTAX(2,2),INUCLEO)+SNDPRP
                        SIGMASFTM(IPTAX(1,2),INUCLEO)=
     &                     SIGMASFTM(IPTAX(1,2),INUCLEO)+SNDPRP
                    ENDIF
                  ENDIF
               ENDIF
  200       CONTINUE
         END IF
  300 CONTINUE
      CALL GPCLOSE(LUSOVE,'DELETE')
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')
C   
C     for magnetizability conversion factor -0.25
C   
      IF (MAGSUS) THEN
        Do I=1,3
          Do J=1,3
            SUSDZD(I,J)=-D025*SUSDZD(I,J)
cDEBUG	
c	    write(LUPRI,*) 'Dentro ABACTODC'
c	    write(LUPRI,*) 'SUSDZD= ',SUSDZD
cDEBUG
          ENDDO
        ENDDO
      ENDIF
      IF (SHIELD) THEN
C   
C     for shieldings conversion factor -0.5
C   
        DO J=1,3
          Do I=1,NCOOR
            SIGMADZ(J,I)=-D05*SIGMADZ(J,I)
            SIGMASFTP(J,I)=D05*SIGMASFTP(J,I)
            SIGMASFTM(J,I)=D05*SIGMASFTM(J,I)
          ENDDO
        ENDDO
      ENDIF
C
      CALL TIMER ('CTOABA',TIMEIN,TIMOUT)
C
      CALL QEXIT('CTOABA')
      RETURN
      END
C...
c    ---------------------------------------
C  /* Deck donsshift */
      SUBROUTINE DONSSHIFT(SFTP,SFTM)
C...
C...
C     This subroutine compute the shift 
C     of the nuclear shielding tensor
C     to have the result of the proprerty with the origin
C     on the atom under examinations and store it in
C     SFTP
C   
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <nuclei.h>
#include <symmet.h>
#include <orgcom.h>
      DIMENSION SFTP(3,3,MXCENT),SFTM(3,3,MXCENT),
     &          SSFT(3,MXCENT)
#include <ibtfun.h>
C     ... compute the (R_cm - R_I) for each atom ...
      JATOM = 0
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         IF (MULT(MULCNT) .EQ. 1) THEN
           JATOM=JATOM+1
           DO i=1,3
             SSFT(I,JATOM) = GAGORG(I) - CORD(I,ICENT)
           END DO
         ELSE
            DO 200 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                  JATOM = JATOM + 1
                 DO i=1,3
                   SSFT(I,JATOM) = GAGORG(I) - 
     &             PT(IBTAND(ISYMAX(I,1),ISYMOP))*CORD(I,ICENT)
                 END DO
               END IF
  200       CONTINUE
         END IF
  100 CONTINUE
      DO 300 INUCLEO = 1,NUCDEP
        DO 300 I =1,3
          SFTP(1,I,INUCLEO) = SFTP(1,I,INUCLEO)*SSFT(3,INUCLEO) 
     &                       -SFTM(1,I,INUCLEO)*SSFT(2,INUCLEO)
          SFTP(2,I,INUCLEO) = SFTP(2,I,INUCLEO)*SSFT(1,INUCLEO) 
     &                       -SFTM(2,I,INUCLEO)*SSFT(3,INUCLEO)
          SFTP(3,I,INUCLEO) = SFTP(3,I,INUCLEO)*SSFT(2,INUCLEO) 
     &                       -SFTM(3,I,INUCLEO)*SSFT(1,INUCLEO)
  300 CONTINUE
      RETURN
      END
c    ---------------------------------------
C  /* Deck getlab */
      SUBROUTINE GETLAB(LABEL,LUNG,LABINT,INTREP,IFILE)
C     07-06-2000 ALig
C
C     This subroutine look in the file IFILE for the property indicated
C     by LABEL and returns the array LABINT and INTREP with the labels 
C     and the symmetry of the desired property. It works only for the 
C     DIPLEN, DIPVEL, ANGMOM, RPSO and PSO but is not difficult to 
C     extend it to other properties. The LUNG is an integer that 
C     contains the real lenght of the LABEL.
C
C     At present it look always in the AOPROPER files but it will be 
C     easy to generalize this subroutine.
C
#include <implicit.h>
#include <mxcent.h>
#include <trkoor.h>
#include <chrxyz.h>
#include <chrnos.h>
cLig DEBUG <>
#include <priunit.h>
      CHARACTER*8 LABEL,LABEL1
      CHARACTER LABINT(*)*8
      DIMENSION INTREP(*)
C
      CALL GPOPEN(IFILE,'AOPROPER','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
C
      REWIND IFILE
      IF ((LABEL(1:LUNG).EQ.'DIPLEN').OR.(LABEL(1:LUNG).EQ.'DIPVEL').OR.
     &   (LABEL(1:LUNG).EQ.'ANGMOM'))  THEN
        DO I=1,3
          LABEL1=CHRXYZ(I)//LABEL(1:LUNG)//' '
          CALL FNDLB3(LABEL1,IVALORE,IFILE)
          INTREP(I) = IVALORE -1
          LABINT(I)=LABEL1
        END DO
      ELSE IF (LABEL(1:LUNG).EQ.'RPSO') THEN
        IL=0
        DO I=1,NCOOR
          DO J=1,3
            LABEL1=CHRNOS(I/10)//CHRNOS(MOD(I,10))//'RPSO '//CHRXYZ(J)
            CALL FNDLB3(LABEL1,IVALORE,IFILE)
            IL=IL+1
            INTREP(IL) = IVALORE -1
            LABINT(IL) = LABEL1
          END DO
        END DO
      ELSE IF (LABEL(1:LUNG).EQ.'PSO') THEN
        IL=0
        DO I=1,NCOOR
          LABEL1='PSO '//CHRNOS(I/10)//CHRNOS(MOD(I,10))//'  '
          CALL FNDLB3(LABEL1,IVALORE,IFILE)
          IL=IL+1
          INTREP(IL) = IVALORE  -1
          LABINT(IL) = LABEL1
        END DO
      END IF
      CALL GPCLOSE(IFILE,'KEEP')
      END
C
c    ---------------------------------------
C  /* Deck todoint */
      LOGICAL FUNCTION TODOINT(LABEL,IFILE)
C     07-06-2000 ALig
C
C     is a non generic function to check if some property (LABEL) are 
C     already stored in file (IFILE).
C     ... it search always in the AOPROPER file
C
C     It could be necessary to correct the PSO part for molecules with 
C     more then 33 atoms
C
#include <implicit.h>
#include <mxcent.h>
#include <trkoor.h>
#include <chrnos.h>
cLig DEBUG <>
#include <priunit.h>
      CHARACTER*8 LABEL,LABEL1
      LOGICAL FNDLAB
C
      TODOINT = .FALSE. 
      NOINT = 0
      CALL GPOPEN(IFILE,'AOPROPER','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND IFILE
      IF (LABEL.EQ.'DIPLEN  ') THEN
         IF (FNDLAB('XDIPLEN ',IFILE)) NOINT=NOINT+1
         IF (FNDLAB('YDIPLEN ',IFILE)) NOINT=NOINT+1
         IF (FNDLAB('ZDIPLEN ',IFILE)) NOINT=NOINT+1
         IF (NOINT.EQ.3) TODOINT=.TRUE.
      ELSE IF (LABEL.EQ.'DIPVEL  ') THEN
         IF (FNDLAB('XDIPVEL ',IFILE)) NOINT=NOINT+1
         IF (FNDLAB('YDIPVEL ',IFILE)) NOINT=NOINT+1
         IF (FNDLAB('ZDIPVEL ',IFILE)) NOINT=NOINT+1
         IF (NOINT.EQ.3) TODOINT=.TRUE.
      ELSE IF (LABEL.EQ.'ANGMOM  ') THEN
         IF (FNDLAB('XANGMOM ',IFILE)) NOINT=NOINT+1
         IF (FNDLAB('YANGMOM ',IFILE)) NOINT=NOINT+1
         IF (FNDLAB('ZANGMOM ',IFILE)) NOINT=NOINT+1
         IF (NOINT.EQ.3) TODOINT=.TRUE.
      ELSE IF (LABEL.EQ.'RPSO    ') THEN
         IF (FNDLAB('01RPSO X',IFILE))  TODOINT = .TRUE.
      ELSE IF (LABEL.EQ.'PSO     ') THEN
         DO I=1,NCOOR
            LABEL1='PSO '//CHRNOS(I/10)//CHRNOS(MOD(I,10))//'  '
            IF (FNDLAB(LABEL1,IFILE)) NOINT=NOINT+1
         END DO  
         IF (NOINT.EQ.NCOOR) TODOINT=.TRUE.
      ENDIF   
      CALL GPCLOSE(IFILE,'KEEP')
      RETURN  
      END     
