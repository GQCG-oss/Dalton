C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
      SUBROUTINE LRSCLIN(WORK,LWORK)
C based on CTOABA 
C...
C...  This subroutine was written by Juan Ignacio Melo using
C...  the subroutine ABACTOCD  as a model (2012)
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
#include "chrxyz.h"
#include "chrnos.h"
C
cxu  sacar lo de abajo
      IPRLNR =14
cxu
      CALL QENTER('LRSCLIN')
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPRLNR .GT. 0) WRITE (LUPRI,'(A,/)')
     *    ' <<<<<<<<<< Output from CTOABA jim >>>>>>>>>> '
C
      IPRRSP = IPRLNR
C
C     Get reference state
C     ===================
C
C     1. Work Allocations:
C
      LUDV   = N2ASHX
      LPVX   = 0
      KFREE  = 1
      LFREE  = LWORK
      write(lupri,*) ' alocando 1 : '
      write(lupri,*) ' ----------------------------'
      write(lupri,*) '     KFREE = 1' 
      write(lupri,*) '     LFREE = LWORK :         ', LWORK 
      write(lupri,*) '                             '     
C      
      write(lupri,*) ' COMMON VARIABLES on LINEAR '
      write(lupri,*) ' ----------------------------'
      write(lupri,*) '     N2ASHX : ni idea ' , N2ASHX
      write(lupri,*) '     NASHT # Active Orbitals = 0 ? :', NASHT
      write(lupri,*) '     LABINT (4*MXCOOR+9) = ', 4*MXCOOR+9
      write(lupri,*) '     MXCOOR              = ', MXCOOR
      write(lupri,*) '     NCMOT = NORB * NORB = ', NCMOT
      write(lupri,*) '   '     
      write(lupri,*) ' memget....  '     

      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KUDV  ,LUDV  ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KPVX  ,LPVX  ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KXINDX,LCINDX,WORK ,KFREE ,LFREE)
c                  TYPE, KBASE, LENGTH, WORK, KFREE, LFREE
c            dimensiona work(KCMO, KCMO+NCMOT)
C
      KWORK1 = KFREE
      WORK1 = LFREE
      write(lupri,*) '   '     
      write(lupri,*) ' AFTER MEMGET  '
      write(lupri,*) ' ----------------------------'
      write(lupri,*) '        KCMO, NCMOT     =  ', KCMO,NCMOT
      write(lupri,*) '        KUDV, LUDV      =  ', KUDV,LUDV
      write(lupri,*) '        KPVX, LPVX      =  ', KPVX,LPVX
      write(lupri,*) '        KPXINDX, LCINDX =  ', KXINDX,KXINDX
      write(lupri,*) '        KWORK1 = KFREE :   ', KFREE 
      write(lupri,*) '        WORK1  = LFREE :   ', LFREE 
      write(lupri,*) '   '     

      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO),WORK(KFREE),LFREE)
C          RD_SIRIFC( KEY ,FOUND,   AMAT   ,  WRK      ,LWRK)
      IF (.NOT.FOUND) CALL QUIT('LRSCLIN error: CMO not on SIRIFC')
Cjim      write(lupri,*)' CMOS :  dan bien !!!' 
Cjim      CALL OUTPUT(WORK(KCMO),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,LUPRI)
      write(lupri,*) '   '     

cx ACA es para alguna capa activa
cx      IF (NASHT .GT. 0) THEN
cx         CALL RD_SIRIFC('DV',FOUND,WORK(KWORK1),WORK(KWORK1),LWORK1)
cx         write(lupri,*)'jim  DV found on RD_SIFC '
cx         IF (.NOT.FOUND)
cx     &      CALL QUIT('ROUTINE error: DV not found on SIRIFC')
cx         CALL DSPTSI(NASHT,WORK(KWORK1),WORK(KUDV))
cx      END IF
C
      ISYM = 1
      write(lupri,*) '   ' 
      write(lupri,*) ' about to call LNRVAR' 
      write(lupri,*) ' ----------------------------'
      write(lupri,*) '       ISYM     : ', ISYM 
      write(lupri,*) '       KWORK1   : ', KWORK1
      write(lupri,*) '       LWORK1   : ', LWORK1
      write(lupri,*) '   '     
      CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK1),LWORK1)
C
c tirar esto
      IPRCIX = 10
c
cd      write(lupri,*) '   ' 
cd      write(lupri,*) ' output of getcix            ' 
      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LWORK1,0)
cd      write(lupri,*) ' ----------------------------'
cd      write(lupri,*) ' INDXCI(*)  :  ' , WORK(KXINDX)
cd      write(lupri,*)  '  '
C
C     SOPPA :
C
cdx      IF (ABASOP) THEN
C
C        Initialize XINDX
C
cdx         CALL DZERO(WORK(KXINDX),LCINDX)
C
C        Find address array's for SOPPA calculation
C
cdx         CALL SET2SOPPA(WORK(KXINDX+KABSAD-1),WORK(KXINDX+KABTAD-1),
cdx     *                  WORK(KXINDX+KIJSAD-1),WORK(KXINDX+KIJTAD-1),
cdx     *                  WORK(KXINDX+KIJ1AD-1),WORK(KXINDX+KIJ2AD-1),
cdx     *                  WORK(KXINDX+KIJ3AD-1),WORK(KXINDX+KIADR1-1))
C
C
cdx         REWIND (LUSIFC)
cdx         IF (CCPPA) THEN
cdx            CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
cdx         ELSE
cdx            CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
cdx         ENDIF
C
C        reads the MP2 or CCSD correlation coefficients into PV
C
cdx         CALL READT (LUSIFC,LPVMAT,WORK(KPVX))
C
cdx         IF (IPRLNR.GT.10) THEN
cdx            IF (CCPPA) THEN
cdx               WRITE(LUPRI,'(/A)')' EXCIT1 : CCSD correlation ',
cdx     &                           'coefficients'
cdx            ELSE
cdx               WRITE(LUPRI,'(/A,A)')' EXCIT1 :',
cdx     &                              ' MP2 correlation coefficients'
cdx            ENDIF
cdx            CALL OUTPUT(WORK(KPVX),1,LPVMAT,1,1,LPVMAT,1,1,LUPRI)
cdx         END IF
C
C        reads the MP2 or CCSD second order one particle density matrix 
C
cdx         CALL READT (LUSIFC,NORBT*NORBT,WORK(KUDV))
C
C        UDV contains the MP2 one-density. Remove the diagonal
C        contribution from the zeroth order. (Added in MP2FAC)
C
cdx         IF (IPRLNR.GT.10) THEN
cdx            IF (CCPPA) THEN
cdx               WRITE(LUPRI,'(/A)')' RSPMC : CCSD density'
cdx            ELSE
cdx               WRITE(LUPRI,'(/A)')' RSPMC : MP2 density'
cdx            END IF
cdx            CALL OUTPUT(WORK(KUDV),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,
cdx     &                  LUPRI)
cdx         END IF
C
cdx         CALL SOPUDV(WORK(KUDV))
cdx      END IF
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
      write(lupri,*) '        '
      write(lupri,*) ' sigo alocando : '
      write(lupri,*) ' ----------------------------'
      write(lupri,*) ' KIDSYM = KWORK1           : ' , KIDSYM
      write(lupri,*) ' KIDADR = KIDSYM + 9MXCENT : ' , KIDADR
      write(lupri,*) ' KWORK2 = KIDADR + 9MXCENT : ' , KWORK2
      write(lupri,*) ' LWORK2 = LWORK - KWORK2   : ' , LWORK2
C
C  Look for NSCO LABELS 
C
      write(lupri,*)
      LABEL1='NSCO'
      write(lupri,*) '@Linear setting  LABEL1 : NSCO '
      IJ = 1
      DO I=1,NCOOR
        DO J=1,3
          LABINT(IJ)= CHRNOS(I/10)//CHRNOS(I/10)//CHRNOS(MOD(I,10))//
     &     'NSCO'//CHRXYZ(J)
          LABSYM(IJ)=1
          IF (LABEL1.EQ.'001NSCOX') write(lupri,*) 'equal'
          write(lupri,*) '   LABEL1 : ', LABINT(IJ), LABSYM(IJ)
          IJ = IJ + 1
        END DO
      END DO
      NLAB = 3*NCOOR
c      NLBTOT = 0
c.      CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
c                             LAB1   ISYM1 -->   LAB2  ISYM2
c     NLAB   = local  numero de labels?
c     NLBTOT = local  
c     LABINT = LABINT(4*MXCOOR+9) character *8
c     WORK
c     LABAPP = common cbiexe, cbiqr, LABAPP(MAXPP) character *8
c              MAXAPP = parameter = 80
c     LABSYM =  common cbiexe, cbiqr,   LABSYM(MAXPP)
C
C  Look for DARWIN LABELS 
C
      write(lupri,*)
      LABEL2='DARWIN  '
      write(lupri,*) '@Linear setting  LABEL2 :', LABEL2
      write(lupri,*) '   LABEL2 : ', LABEL2

cf take as example whe need to calculate integrals.
cf      NCOMP  = 0
cf      NPATOM = 0
cdx      IF (SHIELD) THEN
cdx         IF (TODOINT('RPSO    ',LUPROP)) THEN
cdx           CALL GETLAB('RPSO',4,LABINT,WORK(KIDSYM),LUPROP)
cdx         ELSE
cdx           CALL GET1IN(DUMMY,'RPSO   ',NCOMP,WORK(KWORK2),LWORK2,
cdx     &                 LABINT,WORK(KIDSYM),WORK(KIDADR),
cdx     &                 IDUMMY,.TRUE.,NPATOM,.TRUE.,IPR1IN)
cdx         ENDIF
cdx         NLAB = 3*NCOOR
cdx         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
cdx         write(lupri,*)'RPSO : LABAPP ' , LABAPP
cdx         write(lupri,*)'RPSO : LABSYM ' , LABSYM
cdx         NLBSHIE = NLBTOT
cdx      ENDIF
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
      NABATY = 1
      NABAOP = 1
C
C     Zero the property tensors
C

cdx      IF (MAGSUS) CALL DZERO(SUSDZD,9)

cdx      IF (SHIELD) THEN
cdx        CALL DZERO(SIGMADZ,9*NUCDEP)
cdx        CALL DZERO(SIGMASFTP,9*NUCDEP)
cdx        CALL DZERO(SIGMASFTM,9*NUCDEP)
cdx      ENDIF

C   
C        Loop over the right operators which are the
C        the dipole velocity operators 
C        ===========================================
C   
      LUSOVE = 456
      LUGDVE = 457
      LUREVE = 458
      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
cxv      DO 300 IDIP = 1,3
           IDIP = 1
           LABEL1 = LABINT(1)
c           LABEL1 = '001NSCOX'
           write(lupri,*)' label1 : ', LABEL1 ,'|',LABINT(1) 
C
cx       ISYM=ISYMAX(IDIP,1)+1
         ISYM =1          
C     
          write(lupri,*) ' isymax', ISYMAX(1,1), ISYMAX(2,1)
          write(lupri,*) ' isymax', ISYMAX(2,1), ISYMAX(2,2)
          write(lupri,*) ' isymax', ISYMAX(3,1), ISYMAX(2,1)
          write(lupri,*) ' isymaxxxxx', ISYMAX
          write(lupri,*)' ISYM    :              ',ISYM
c        set variables for response module
         write(lupri,*) ' about to call LNRVAR' 
         write(lupri,*) ' ----------------------------'
         CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK2),LWORK2)
C
         IF (NFRVAL.GT.0) THEN
            write(lupri,*)' antes NFRVAL:',NFRVAL 
            write(lupri,*)'       common de cbilnr.h ' 
C
C           3. Work Allocations:
C
            KGD1   = KWORK1
            KWRKG1 = KGD1
            LWRKG1 = LWORK - KWRKG1
            KSLV   = KGD1 + 2*NVARPT
            KLAST  = KSLV + 2*NVARPT
            write(lupri,*) ' KLAST ',KLAST
            write(lupri,*) ' LWORK ',LWORK
            write(lupri,*) ' KLAST GT LWORK ! error '
            IF (KLAST.GT.LWORK) CALL STOPIT('LINEARjim',' ',KLAST,LWORK)
            KWRK = KLAST
            LWRK = LWORK - KLAST + 1
            write(lupri,*)' NVARPT :', NVARPT 

C
C           Find right hand side for right operator and write to file
C           =========================================================
C
            KSYMOP = ISYM
            TRPLET = .FALSE.
            IPRSP = 15 
            ANTSYM = 0
C   
            write(lupri,*) ' about to call GETPV'  
            write(*,*) '@linearLRESC  ANTSYM', ANTSYM  
C jim sym to be read on luprop
C ANTSYM : matrix symmetry of PRPMO matrix
C          (1: symmetric, -1: antisymmetric, 0: unknown)

            CALL GETGPV(LABEL1,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &           WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
            REWIND LUGDVE
            CALL WRITT(LUGDVE,2*NVARPT,WORK(KWRKG1))
cvx            IF (IPRLNR.GT.-1) THEN
               WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL1
               CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,2,1,LUPRI)
cvx            ENDIF
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
C  en este caso todos los NSCO 
C
cx            DO 200 IPRLBL = 4, NLBTOT
C
C              Find label and symmetry of the left side operator
C
         write(lupri,*)' antes del lio22222.?' 
cx               LABEL2 = LABAPP(IPRLBL)
               LABEL2 = 'DARWIN  ' 
cd               LABEL2 = 'DARWIN '

               write(lupri,*)' LABAPP IPRLBL ' , LABEL2
             write(lupri,*) '  KSYM   = LABSYM(IPRLBL) ', LABSYM(IPRLBL)
cx               KSYM   = LABSYM(IPRLBL)
               KSYM   = 1

C
C              If symmetry of right operator equals symmetry of
C              the left operator, that is if ISYM = KSYM, then
C              ================================================
C              (otherwise 2. order property SNDPRP is zero)
C
cx               IF (KSYM.EQ.ISYM) THEN
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
                  IF (IPRLNR.GT.-1) THEN
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
                  IF (IPRLNR.GT.-1) THEN
                     WRITE (LUPRI,'(2A)') 'Solution Vector, label: ',
     &                                    LABEL1
                     CALL OUTPUT(WORK(KSLV),1,NVARPT,1,2,NVARPT,2,1,
     &                           LUPRI)
                  ENDIF
C
                  SNDPRP = DDOT(2*NVARPT,WORK(KSLV),1,WORK(KGD1),1)
C
                  IF (IPRLNR.GT.-1) THEN
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
cx                  IF (IPRLBL.LT.NLBSHIS) THEN 
cx                     IF (LABEL2(1:2).EQ.'XX') THEN  
cx                        IF (LABEL1(1:1).EQ.'Z') THEN  
cx                           SUSDZD(IPTAX(1,2),IPTAX(2,2))= 
cx     &                        SUSDZD(IPTAX(1,2),IPTAX(2,2))+ SNDPRP
cx                        ELSE IF (LABEL1(1:1).EQ.'Y') THEN 
cx                           SUSDZD(IPTAX(1,2),IPTAX(3,2))= 
cx     &                        SUSDZD(IPTAX(1,2),IPTAX(3,2))- SNDPRP
cx                        ENDIF
cx                     ELSE IF (LABEL2(1:2).EQ.'XY') THEN
cx                        IF (LABEL1(1:1).EQ.'Z') THEN
cx                           SUSDZD(IPTAX(2,2),IPTAX(2,2))=
cx     &                        SUSDZD(IPTAX(2,2),IPTAX(2,2))+ SNDPRP
cx                        ELSE IF (LABEL1(1:1).EQ.'Y') THEN
cx                           SUSDZD(IPTAX(2,2),IPTAX(3,2))=
cx     &                        SUSDZD(IPTAX(2,2),IPTAX(3,2))- SNDPRP
cx                        ENDIF
cx                     ELSE IF (LABEL2(1:2).EQ.'XZ') THEN
cx                        IF (LABEL1(1:1).EQ.'Z') THEN
cx                           SUSDZD(IPTAX(3,2),IPTAX(2,2))=
cx     &                        SUSDZD(IPTAX(3,2),IPTAX(2,2))+ SNDPRP
cx                        ELSE IF (LABEL1(1:1).EQ.'Y') THEN
cx                           SUSDZD(IPTAX(3,2),IPTAX(3,2))=
cx     &                        SUSDZD(IPTAX(3,2),IPTAX(3,2))- SNDPRP
cx                        ENDIF
cx                     ELSE IF (LABEL2(1:2).EQ.'YX') THEN
cx                        IF (LABEL1(1:1).EQ.'X') THEN
cx                           SUSDZD(IPTAX(1,2),IPTAX(3,2))=
cx     &                        SUSDZD(IPTAX(1,2),IPTAX(3,2))+ SNDPRP
cx                        ELSE IF (LABEL1(1:1).EQ.'Z') THEN
cx                           SUSDZD(IPTAX(1,2),IPTAX(1,2))=
cx     &                        SUSDZD(IPTAX(1,2),IPTAX(1,2))- SNDPRP
cx                        ENDIF
cx                     ELSE IF (LABEL2(1:2).EQ.'YY') THEN
cx                        IF (LABEL1(1:1).EQ.'X') THEN
cx                           SUSDZD(IPTAX(2,2),IPTAX(3,2))=
cx     &                        SUSDZD(IPTAX(2,2),IPTAX(3,2))+ SNDPRP
cx                        ELSE IF (LABEL1(1:1).EQ.'Z') THEN
cx                           SUSDZD(IPTAX(2,2),IPTAX(1,2))=
cx     &                        SUSDZD(IPTAX(2,2),IPTAX(1,2))- SNDPRP
cx                        ENDIF
cx                     ELSE IF (LABEL2(1:2).EQ.'YZ') THEN
cx                        IF (LABEL1(1:1).EQ.'X') THEN
cx                           SUSDZD(IPTAX(3,2),IPTAX(3,2))=
cx     &                        SUSDZD(IPTAX(3,2),IPTAX(3,2))+ SNDPRP
cx                        ELSE IF (LABEL1(1:1).EQ.'Z') THEN
cx                           SUSDZD(IPTAX(3,2),IPTAX(1,2))=
cx     &                        SUSDZD(IPTAX(3,2),IPTAX(1,2))- SNDPRP
cx                        ENDIF
cx                     ELSE IF (LABEL2(1:2).EQ.'ZX') THEN
cx                        IF (LABEL1(1:1).EQ.'Y') THEN
cx                           SUSDZD(IPTAX(1,2),IPTAX(1,2))=
cx     &                        SUSDZD(IPTAX(1,2),IPTAX(1,2))+ SNDPRP
cx                        ELSE IF (LABEL1(1:1).EQ.'X') THEN
cx                           SUSDZD(IPTAX(1,2),IPTAX(2,2))=
cx     &                        SUSDZD(IPTAX(1,2),IPTAX(2,2))- SNDPRP
cx                        ENDIF
cx                     ELSE IF (LABEL2(1:2).EQ.'ZY') THEN
cx                        IF (LABEL1(1:1).EQ.'Y') THEN
cx                           SUSDZD(IPTAX(2,2),IPTAX(1,2))=
cx     &                        SUSDZD(IPTAX(2,2),IPTAX(1,2))+ SNDPRP
cx                        ELSE IF (LABEL1(1:1).EQ.'X') THEN
cx                           SUSDZD(IPTAX(2,2),IPTAX(2,2))=
cx     &                        SUSDZD(IPTAX(2,2),IPTAX(2,2))- SNDPRP
cx                        ENDIF
cx                     ELSE IF (LABEL2(1:2).EQ.'ZZ') THEN
cx                        IF (LABEL1(1:1).EQ.'Y') THEN
cx                           SUSDZD(IPTAX(3,2),IPTAX(1,2))=
cx     &                        SUSDZD(IPTAX(3,2),IPTAX(1,2))+ SNDPRP
cx                        ELSE IF (LABEL1(1:1).EQ.'X') THEN
cx                           SUSDZD(IPTAX(3,2),IPTAX(2,2))=
cx     &                        SUSDZD(IPTAX(3,2),IPTAX(2,2))- SNDPRP
cx                        ENDIF
cx                     ENDIF
cx                  ENDIF
C
C                 Nuclear Shieldings
C                 ------------------
cx                  IF ((IPRLBL.GE.NLBSHIS).AND.
cx     &                (IPRLBL.LE.NLBSHIE))  THEN
cx                     INUCLEO = (IPRLBL-NLBSHIS)/3+1
cx                     IF (LABEL1.EQ.'XDIPVEL') THEN
cx                        IF (LABEL2(8:8).EQ.'Y') THEN
cx                           SIGMADZ(IPTAX(3,2),INUCLEO)=
cx     &                        SIGMADZ(IPTAX(3,2),INUCLEO)+SNDPRP
cx                        ELSE IF (LABEL2(8:8).EQ.'Z') THEN
cx                           SIGMADZ(IPTAX(2,2),INUCLEO)=
cx     &                        SIGMADZ(IPTAX(2,2),INUCLEO)-SNDPRP
cx                        ENDIF
cx                     ELSE IF (LABEL1.EQ.'YDIPVEL') THEN
cx                        IF (LABEL2(8:8).EQ.'Z') THEN
cx                           SIGMADZ(IPTAX(1,2),INUCLEO)=
cx     &                        SIGMADZ(IPTAX(1,2),INUCLEO)+SNDPRP
cx                        ELSE IF (LABEL2(8:8).EQ.'X') THEN
cx                           SIGMADZ(IPTAX(3,2),INUCLEO)=
cx     &                        SIGMADZ(IPTAX(3,2),INUCLEO)-SNDPRP
cx                        ENDIF
cx                     ELSE IF (LABEL1.EQ.'ZDIPVEL') THEN
cx                        IF (LABEL2(8:8).EQ.'X') THEN
cx                           SIGMADZ(IPTAX(2,2),INUCLEO)=
cx     &                        SIGMADZ(IPTAX(2,2),INUCLEO)+SNDPRP
cx                        ELSE IF (LABEL2(8:8).EQ.'Y') THEN
cx                           SIGMADZ(IPTAX(1,2),INUCLEO)=
cx     &                        SIGMADZ(IPTAX(1,2),INUCLEO)-SNDPRP
cx                        ENDIF
cx                     ENDIF
cx                  ENDIF
C
C                 Nuclear Shielding Corrections
C                 -----------------------------
cx                  IF (IPRLBL.GT.NLBSHIE) THEN
cx                     INUCLEO = (IPRLBL-NLBSHIE)
cx                     IF (LABEL1.EQ.'XDIPVEL') THEN
cx                        SIGMASFTP(IPTAX(3,2),INUCLEO)=
cx     &                     SIGMASFTP(IPTAX(3,2),INUCLEO)+SNDPRP
cx                        SIGMASFTM(IPTAX(2,2),INUCLEO)=
cx     &                     SIGMASFTM(IPTAX(2,2),INUCLEO)+SNDPRP
cx                     ELSE IF (LABEL1.EQ.'YDIPVEL') THEN
cx                        SIGMASFTP(IPTAX(1,2),INUCLEO)=
cx     &                     SIGMASFTP(IPTAX(1,2),INUCLEO)+SNDPRP
cx                        SIGMASFTM(IPTAX(3,2),INUCLEO)=
cx     &                     SIGMASFTM(IPTAX(3,2),INUCLEO)+SNDPRP
cx                     ELSE IF (LABEL1.EQ.'ZDIPVEL') THEN
cx                        SIGMASFTP(IPTAX(2,2),INUCLEO)=
cx     &                     SIGMASFTP(IPTAX(2,2),INUCLEO)+SNDPRP
cx                        SIGMASFTM(IPTAX(1,2),INUCLEO)=
cx     &                     SIGMASFTM(IPTAX(1,2),INUCLEO)+SNDPRP
cx                    ENDIF
cx                  ENDIF
cx               ENDIF
  200       CONTINUE
         END IF
cxv  300 CONTINUE
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
      CALL TIMER ('LRSCLIN',TIMEIN,TIMOUT)
C
      CALL QEXIT('LRSCLIN')
      RETURN
      END
C...
c    ---------------------------------------
C  /* Deck donsshift */
cx      SUBROUTINE DONSSHIFT(SFTP,SFTM)
C...
C...
C     This subroutine compute the shift 
C     of the nuclear shielding tensor
C     to have the result of the proprerty with the origin
C     on the atom under examinations and store it in
C     SFTP
C   
cx#include <implicit.h>
cx#include <mxcent.h>
cx#include <maxaqn.h>
cx#include <maxorb.h>
cx#include <nuclei.h>
cx#include <symmet.h>
cx#include <orgcom.h>
cx      DIMENSION SFTP(3,3,MXCENT),SFTM(3,3,MXCENT),
cx     &          SSFT(3,MXCENT)
cx#include <ibtfun.h>
C     ... compute the (R_cm - R_I) for each atom ...
cx      JATOM = 0
cx      DO 100 ICENT = 1, NUCIND
cx         MULCNT = ISTBNU(ICENT)
cx         IF (MULT(MULCNT) .EQ. 1) THEN
cx           JATOM=JATOM+1
cx           DO i=1,3
cx             SSFT(I,JATOM) = GAGORG(I) - CORD(I,ICENT)
cx           END DO
cx         ELSE
cx            DO 200 ISYMOP = 0, MAXOPR
cx               IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
cx                  JATOM = JATOM + 1
cx                 DO i=1,3
cx                   SSFT(I,JATOM) = GAGORG(I) - 
cx     &             PT(IBTAND(ISYMAX(I,1),ISYMOP))*CORD(I,ICENT)
cx                 END DO
cx               END IF
cx  200       CONTINUE
cx         END IF
cx  100 CONTINUE
cx      DO 300 INUCLEO = 1,NUCDEP
cx        DO 300 I =1,3
cx          SFTP(1,I,INUCLEO) = SFTP(1,I,INUCLEO)*SSFT(3,INUCLEO) 
cx     &                       -SFTM(1,I,INUCLEO)*SSFT(2,INUCLEO)
cx          SFTP(2,I,INUCLEO) = SFTP(2,I,INUCLEO)*SSFT(1,INUCLEO) 
cx     &                       -SFTM(2,I,INUCLEO)*SSFT(3,INUCLEO)
cx          SFTP(3,I,INUCLEO) = SFTP(3,I,INUCLEO)*SSFT(2,INUCLEO) 
cx     &                       -SFTM(3,I,INUCLEO)*SSFT(1,INUCLEO)
cx  300 CONTINUE
cx      RETURN
cx      END
c    ---------------------------------------
C
c    ---------------------------------------
C  /* Deck todoint */
cx      LOGICAL FUNCTION TODOINT(LABEL,IFILE)
C     07-06-2000 ALig
C
C     is a non generic function to check if some property (LABEL) are 
C     already stored in file (IFILE).
C     ... it search always in the AOPROPER file
C
C     It could be necessary to correct the PSO part for molecules with 
C     more then 33 atoms
C
cx#include <implicit.h>
cx#include <mxcent.h>
cx#include <trkoor.h>
cx#include <chrnos.h>
cLig DEBUG <>
cx#include <priunit.h>
cx      CHARACTER*8 LABEL,LABEL1
cx      LOGICAL FNDLAB
C
cx     TODOINT = .FALSE. 
cx      NOINT = 0
cx      CALL GPOPEN(IFILE,'AOPROPER','OLD',' ','UNFORMATTED',IDUMMY,
cx     &            .FALSE.)
cx      REWIND IFILE
cx      IF (LABEL.EQ.'DIPLEN  ') THEN
cx         IF (FNDLAB('XDIPLEN ',IFILE)) NOINT=NOINT+1
cx         IF (FNDLAB('YDIPLEN ',IFILE)) NOINT=NOINT+1
cx         IF (FNDLAB('ZDIPLEN ',IFILE)) NOINT=NOINT+1
cx         IF (NOINT.EQ.3) TODOINT=.TRUE.
cx      ELSE IF (LABEL.EQ.'DIPVEL  ') THEN
cx         IF (FNDLAB('XDIPVEL ',IFILE)) NOINT=NOINT+1
cx         IF (FNDLAB('YDIPVEL ',IFILE)) NOINT=NOINT+1
cx         IF (FNDLAB('ZDIPVEL ',IFILE)) NOINT=NOINT+1
cx         IF (NOINT.EQ.3) TODOINT=.TRUE.
cx      ELSE IF (LABEL.EQ.'ANGMOM  ') THEN
cx         IF (FNDLAB('XANGMOM ',IFILE)) NOINT=NOINT+1
cx         IF (FNDLAB('YANGMOM ',IFILE)) NOINT=NOINT+1
cx         IF (FNDLAB('ZANGMOM ',IFILE)) NOINT=NOINT+1
cx         IF (NOINT.EQ.3) TODOINT=.TRUE.
cx      ELSE IF (LABEL.EQ.'RPSO    ') THEN
cx         IF (FNDLAB('01RPSO X',IFILE))  TODOINT = .TRUE.
cx      ELSE IF (LABEL.EQ.'PSO     ') THEN
cx         DO I=1,NCOOR
cx            LABEL1='PSO '//CHRNOS(I/10)//CHRNOS(MOD(I,10))//'  '
cx            IF (FNDLAB(LABEL1,IFILE)) NOINT=NOINT+1
cx         END DO  
cx         IF (NOINT.EQ.NCOOR) TODOINT=.TRUE.
cx      ENDIF   
cx      CALL GPCLOSE(IFILE,'KEEP')
cx      RETURN  
cx      END     
