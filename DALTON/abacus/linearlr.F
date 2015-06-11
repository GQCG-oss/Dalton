!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C
      SUBROUTINE LINEARLR(KEY,WORK,LWORK)
C  
C...
C...  This subroutine was written by Juan Ignacio Melo using
C...  the subroutine ABACTOCD  as a model (2012)
C
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "trkoor.h"
c#include "sigma.h"
#include "maxorb.h"
#include "iratdef.h"
#include "priunit.h"
#include "cbilnr.h"
c#include "suscpt.h"
#include "infpri.h"
      LOGICAL CICLC, HFCLC, TRIPLE, EXECLC, FOUND
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL1,LABEL2,LISTA1(4*MXCOOR+9),LISTA2(4*MXCOOR+9)
      CHARACTER*6 LABEL22
      CHARACTER*4 KEY
      CHARACTER*3 char
      PARAMETER (D05=0.5D0,D025=0.25)
      LOGICAL TODOINT
C
#include "cbiexc.h"
#include "inflin.h"
#include "infvar.h"
#include "infdim.h"
#include "inforb.h"
#include "nuclei.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "maxmom.h"
#include "maxaqn.h"
#include "symmet.h"
#include "abainf.h"
#include "gnrinf.h"
c#include "infsop.h"

C
#include "lrescinf.h"
#include "chrxyz.h"
#include "chrnos.h"
#include "orgcom.h"
C
cxu  sacar lo de abajo
cc      IPRLNR = 5
      IPRLNR = JIMPRT
cxu
C
C
C  LINEAR LRESC SINGLET ROUTINE
C   KEY = 'DIAM'  DIAMAGNETIC <<DIA,MV>>  
C         'DIAD'  DIAMAGNETIC <<DIA,DW>> 
C         'PSOK'  PARAMAGNETIC with built-in integrals
C         'PSKI'  PARAMAGNETIC with angkin routine <<PSO.1.Kinenerg;PSO>>
C         'LKIN'  PARAMAGNETIC
C         'FCBS' 'SDBS'
C         'FCZK' 'FCZK'    
C
      CALL QENTER('LRSCLIN')
      CALL TIMER('START ',TIMEIN,TIMOUT)

      IF (JIMPRT .GE. 2) THEN
      WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
      WRITE(LUPRI,*)
         IF (KEY.EQ.'DIAM') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Diamagnetic Second Order Singlet : DIAM '
         IF (KEY.EQ.'DIAD') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Diamagnetic Second Order Singlet : DIAD '
         IF (KEY.EQ.'PSOK') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Paramagnetic Second Order Singlet : PSOK '
         IF (KEY.EQ.'PSOKI') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Paramagnetic Second Order Singlet : PSOKI'
         IF (KEY.EQ.'LKIN') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Paramagnetic Second Order Singlet : LKIN '
         IF (KEY.EQ.'FCZK') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Paramagnetic Second Order Triplet : FCZK'
      IF (KEY.EQ.'SDZK') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Paramagnetic Second Order Triplet : SDZK'
      IF (KEY.EQ.'FCBS') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Paramagnetic Second Order Triplet : FCBS'
      IF (KEY.EQ.'SDBS') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Paramagnetic Second Order Triplet : SDBS'
      WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
      END IF
C
       IPRRSP = JIMPRT-1

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
      IF (JIMPRT.GT.2) Then 
          write(lupri,'(A,3F12.8)') ' orgcom.h : GAGORG :', GAGORG
          write(lupri,'(A,3F12.8)') '            ORIGIN :', ORIGIN
          write(lupri,'(A,3F12.8)') '            CMXYZ  :', CMXYZ
          write(lupri,*)
          write(lupri,*) ' alocando 1 : ANTES MEMGET :'
          write(lupri,*) ' ----------------------------'
          write(lupri,*) '     KFREE = 1' 
          write(lupri,*) '     LFREE = LWORK :         ', LWORK 
          write(lupri,*) '                             '     
C      
          write(lupri,*) ' COMMON VARIABLES on LINEAR '
          write(lupri,*) ' ----------------------------'
          write(lupri,*) '     N2ASHX : ni idea ' , N2ASHX
          write(lupri,*) '     NASHT # Active Orbitals = 0 ? :', NASHT
          write(lupri,*) '     LISTA1y2 (4*MXCOOR+9) = ', 4*MXCOOR+9
          write(lupri,*) '     MXCOOR              = ', MXCOOR
          write(lupri,*) '     NCMOT = NORB * NORB = ', NCMOT
          write(lupri,*) '   '     
          write(lupri,*) ' memget....  '     
      ENDIF 
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KUDV  ,LUDV  ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KPVX  ,LPVX  ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KXINDX,LCINDX,WORK ,KFREE ,LFREE)
c                  TYPE, KBASE, LENGTH, WORK, KFREE, LFREE
c            dimensiona work(KCMO, KCMO+NCMOT)
C
      KWORK1 = KFREE
      LWORK1 = LFREE
      IF (JIMPRT.GT.5) Then 
          write(lupri,*) '   '
          write(lupri,*) ' AFTER MEMGET  '
          write(lupri,*) ' ----------------------------'
          write(lupri,*) '        KCMO, NCMOT     =  ', KCMO,NCMOT
          write(lupri,*) '        KUDV, LUDV      =  ', KUDV,LUDV
          write(lupri,*) '        KPVX, LPVX      =  ', KPVX,LPVX
          write(lupri,*) '        KPXINDX, LCINDX =  ', KXINDX,KXINDX
          write(lupri,*) '        KWORK1 = KFREE :   ', KFREE
          write(lupri,*) '        LWORK1 = LFREE :   ', LFREE
          write(lupri,*) '   '
      ENDIF

      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO),WORK(KFREE),LFREE)
C          RD_SIRIFC( KEY ,FOUND,   AMAT   ,  WRK      ,LWRK)
      IF (.NOT.FOUND) CALL QUIT('LRSCLIN error: CMO not on SIRIFC')
cx      write(lupri,*)' CMOS :  dan bien !!!'
cx      Do i =1, NORBT*NORBT
cx      write(LUPRI,*) ' orb ', I,WORK(KCMO+I)
cx      enddo
cx      CALL OUTPUT(WORK(KCMO),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,LUPRI)
cx      write(lupri,*) '   '

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
      IF (JIMPRT.GT.5) Then
          write(lupri,*) '   '
          write(lupri,*) ' about to call LNRVAR'
          write(lupri,*) ' ----------------------------'
          write(lupri,*) '           ISYM   : ', ISYM
          write(lupri,*) '   KWORK1=KFREE   : ', KWORK1
          write(lupri,*) '   LWORK1=LFREE   : ', LWORK1
          write(lupri,*) '   '
      ENDIF
      CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK1),LWORK1)
C
c tirar esto
      IPRCIX = -1
c
      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LWORK1,0)
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
      LWORK2 = LWORK1  - KWORK2
      IF (JIMPRT.GT.5) Then
          write(lupri,*) '        '
          write(lupri,*) ' sigo alocando : '
          write(lupri,*) ' ----------------------------'
          write(lupri,*) ' KIDSYM = KWORK1           : ' , KIDSYM
          write(lupri,*) ' KIDADR = KIDSYM + 9MXCENT : ' , KIDADR
          write(lupri,*) ' KWORK2 = KIDADR + 9MXCENT : ' , KWORK2
          write(lupri,*) ' LWORK2 = LWORK - KWORK2   : ' , LWORK2
      ENDIF
C
cxz       NLBTOT = 0
cxz       NLBSHIS = 4
C
cxz       NCOMP  = 0
cxz       NPATOM = 0
cLig  <> added the TODOINT to see if the property was already in the file
cxz       IF (TODOINT('DIPVEL  ',LUPROP)) THEN
cxz         CALL GETLAB('DIPVEL',6,LABINT,WORK(KIDSYM),LUPROP)
cxz       ELSE
cxz         CALL GET1IN(DUMMY,'DIPVEL ',NCOMP,WORK(KWORK2),LWORK2,
cxz      &              LABINT,WORK(KIDSYM),WORK(KIDADR),
cxz      &              IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
cxz      &              IPR1IN)
cxz       ENDIF
cxz       NLAB = 3

cxz       CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
cxz                        LAB1   ISYM1 -->   LAB2  ISYM2
cxz     NLAB   = local  numero de labels?
cxz     NLBTOT = local
cxz     LABINT = LABINT(4*MXCOOR+9) LISTA1Y2 character *8
cxz     WORK
cxz     LABAPP = common cbiexe, cbiqr, LABAPP(MAXPP) character *8
cxz              MAXAPP = parameter = 80
cxz     LABSYM =  common cbiexe, cbiqr,   LABSYM(MAXPP)

cxz       IF (MAGSUS) THEN
cxz          NCOMP  = 0
cxz          CALL GET1IN(DUMMY,'RANGMO ',NCOMP,WORK(KWORK2),LWORK2,
cxz      &               LABINT,WORK(KIDSYM),WORK(KIDADR),
cxz      &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
cxz      &               IPR1IN)
cxz          NLAB = 9
cxz          CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
cxz          NLBSHIS = 13
cxz       ENDIF
C  blanquear todas las labels antes
       NLAB1 = 0
       NLAB2 = 0 
c       LISTA1
c       LISTA2
c       LABSYM 
c       LABAPP es un common q se usa en rsp ??


C ==============================================================================
C
C  Staritng Labels stuff
C
C ==============================================================================

cb       write(lupri,*)' selected atom is :', LRATOM
cb       write(lupri,*)' LRATOM / 100' , LRATOM / 100
cb       write(lupri,*)' LRATOM / 10' , LRATOM / 10
 
       npos1 = 3*LRATOM-2
       aa = LRATOM / 100
       bb = LRATOM /10
cb       write(lupri,*)
cb       IF (aa.GE.1) THEN
cb               write(lupri,*)'LRATOM tiene 3 digits'
cb          ELSE IF (bb.GE.1) THEN
cb               write(lupri,*)'LRATOM tiene 2 digitos'
cb          ELSE IF (bb.LT.1) THEN
cb               write(lupri,*)'LRATOM tiene 1 digito'
cb       ENDIF
cb       write(lupri,*)
C
C===============================================================================
C  SINGLET CALCULATIONS : D1S and P1S 
C ==============================================================================
      IF((KEY.EQ."DIAM").OR.(KEY.EQ."DIAD")) THEN 
C        Look for LABELS : NSCO  
C        -----------------------
         LABEL1='NSCO'
         IJ = 1
         DO I=npos1,npos1+2
           DO J=1,3
             LISTA1(IJ)= CHRNOS(I/100)//CHRNOS(I/10)//CHRNOS(MOD(I,10))
     &       //'NSCO'//CHRXYZ(J)
             LABSYM(IJ)=1
             IJ = IJ + 1
           END DO
         END DO
         NLAB1 = 9
cx         NLAB1 = 3*NCOOR
C        Look for LABELS : DARWIN, MASSVELO
C        -----------------------------------
         IF (KEY.EQ."DIAD") LISTA2(1) ='DARWIN  '
         IF (KEY.EQ."DIAM") LISTA2(1) ='MASSVELO'
         NLAB2 = 1
C --------------------------------------------------------------
      ELSE IF (KEY.EQ."LKIN") THEN
C        Look for LABELS : PSO 
C        ----------------------       
         LABEL1='PSO'
         IJ = 1
         DO I=npos1,npos1+2
             LISTA1(IJ)= 'PSO '//CHRNOS(I/100)//CHRNOS(I/10)//
     &                    CHRNOS(MOD(I,10))
             LABSYM(IJ)=1
             IJ = IJ + 1
         END DO
         NLAB1 = 3
C        Look for LABELS : OZ-KE
C        ------------------------       
         LABEL2='OZKE   '
         DO J=1,3
           LISTA2(J)= CHRXYZ(J)//'OZKE   '
           LABSYM(J)=1 
         END DO 
         NLAB2 = 3
C -----------------------------------------------------------------
      ELSE IF (KEY.EQ."PSKI") THEN  !We use jmelo integrals. 


C AngKin makes L-PsoKin integrals from PSO and ANgMOM
c additionally also makes Lkin for free
C
           CALL ANGKIN(WORK(KCMO),WORK,LWORK2)
C
c        we repeat the labeling here, for simplicity.
C        Look for LABELS : PSO-KE : jmelo way
C        -------------------------------------       
         LABEL1='PSOKI'
         IJ = 1
         DO I=npos1,npos1+2
             LISTA1(IJ)= 'PSOKI'//CHRNOS(I/100)//CHRNOS(I/10)//
     &                    CHRNOS(MOD(I,10))
             LABSYM(IJ)=1
             IJ = IJ + 1
         END DO
         NLAB1 = 3
C        Look for LABELS : ANGMOM
C        -------------------------       
         LABEL2='ANGMOM'
         DO J=1,3
           LISTA2(J)= CHRXYZ(J)//'ANGMOM '
           LABSYM(J)=1
         END DO
         NLAB2 = 3
C
C -----------------------------------------------------------------
      ELSE IF (KEY.EQ."PSOK") THEN  ! we will use dalton built in integrals
C
C        Look for LABELS : PSO-KE : jmelo way
C        -------------------------------------       
         LABEL1='PSOKE '
         IJ = 1
         DO I=npos1,npos1+2
             LISTA1(IJ)= 'PSOKE'//CHRNOS(I/100)//CHRNOS(I/10)//
     &                    CHRNOS(MOD(I,10))
             LABSYM(IJ)=1
             IJ = IJ + 1
         END DO
         NLAB1 = 3
C        Look for LABELS : ANGMOM
C        -------------------------       
         LABEL2='ANGMOM'
         DO J=1,3
           LISTA2(J)= CHRXYZ(J)//'ANGMOM '
           LABSYM(J)=1
         END DO
         NLAB2 = 3
      END IF

C===============================================================================
C  TRIPLET CALCULATIONS : FcKin, SdKin and FcBso, SdBso
C ==============================================================================
C -----------------------------------------------------------------
       IF((KEY.EQ."FCZK").OR.(KEY.EQ."FCBS")) THEN
C
C        Look for LABELS : FERMI  
C        ------------------------
         LABEL22='FC '//NAMN(LRATOM)(1:3)
C jim : just in case integral label changes to 3 digits on atom list 
c         IF (aa.GE.1) THEN
c             write(char,'(i3)') LRATOM
c             LISTA2(1)= LABEL22//char
c          ELSE IF (bb.GE.1) 
          IF (bb.GE.1) THEN
             write(char,'(i2)') LRATOM
             LISTA2(1)= LABEL22//char
          ELSE IF (bb.LT.1) THEN
             write(char,'(i1)') LRATOM
             LISTA2(1)= LABEL22//'0'//char
         ENDIF
         LABSYM(1)=1
         NLAB2 = 1
C         
         IF(KEY.EQ."FCZK") THEN
C
C        Look for LABELS : KINENERG
C        -----------------------------------
            LISTA1(1) ='KINENERG'
            NLAB1 = 1
C
C        Look for LABELS : BSO
C        -----------------------------------
         ELSE IF(KEY.EQ."FCBS") THEN
            LISTA1(1) ='SOMF  XX'
            LISTA1(2) ='SOMF  YX'
            LISTA1(3) ='SOMF  ZX'
            LISTA1(4) ='SOMF  YY'
            LISTA1(5) ='SOMF  ZY'
            LISTA1(6) ='SOMF  ZZ'
            NLAB1 = 6
         ENDIF
      ENDIF
C --------------------------------------------------------------
      IF((KEY.EQ."SDZK").OR.(KEY.EQ."SDBS")) THEN
C
C        Look for LABELS : SD
C        ------------------------
         LABEL2='SD'
         IJ = 1
         DO I=npos1,npos1+2
           DO J=1,3
             LISTA2(IJ)= 'SD '//CHRNOS(J/100)//CHRNOS(J/10)//
     &                    CHRNOS(MOD(I,10))//' '//CHRXYZ(-J)
             LABSYM(IJ)=1
             IJ = IJ + 1
           END DO
         END DO
         NLAB2 = 9
C
C        Look for LABELS : dd/didj
C        -----------------------------------
         IF(KEY.EQ."SDZK") THEN
            LISTA1(1) ='dd/dxdx'
            LISTA1(2) ='dd/dxdy'
            LISTA1(3) ='dd/dxdz'
            LISTA1(4) ='dd/dydy'
            LISTA1(5) ='dd/dydz'
            LISTA1(6) ='dd/dzdz'
            NLAB1 = 6
            
         ELSE IF(KEY.EQ."SDBS") THEN
C        Look for LABELS : BSO
C        -----------------------------------
            LISTA1(1) ='SOMF  XX'
            LISTA1(2) ='SOMF  YX'
            LISTA1(3) ='SOMF  ZX'
            LISTA1(4) ='SOMF  YY'
            LISTA1(5) ='SOMF  ZY'
            LISTA1(6) ='SOMF  ZZ'
            NLAB1 = 6
         ENDIF
      ENDIF
C ---------------------------------------------------------------------
C
C  Print Section for LABELS on LISTAS 
C  ----------------------------------------------
C        
      IF (JIMPRT.GE.2) THEN 
        WRITE(LUPRI,*)' @LinearLR setting  LABEL1 :'
        DO i =1, nlab1
          WRITE(LUPRI,*)'   LABEL1 :', LISTA1(I)
        ENDDO
        WRITE(LUPRI,*)
C 
        WRITE(LUPRI,*)' @LinearLR setting  LABEL2 :'
        DO i =1, nlab2
          WRITE(LUPRI,*)'   LABEL2 :', LISTA2(I)
        ENDDO
      ENDIF
C
C
C ---------------------------------------------------------------------
C
C     Set variables for ABARSP and logicals
C
      CICLC  = .FALSE.      ! TRUE for CI calculations
      HFCLC  = NASHT .LE. 1 ! .T. RHF-closed shell or 1e in one active orbital
      TRIPLE = .FALSE.      ! .T. for triplet perturbation operators
      EXECLC = .FALSE.      ! false for linear response equations
      IF(KEY.EQ."DIAM".OR.KEY.EQ."DIAD") THEN 
         NABATY = 1      ! = 1 for real operators .. -1 for imm. op.
         ELSE IF(KEY.EQ."LKIN") THEN 
         NABATY = -1     ! = 1 for real operators .. -1 for imm. op. 
         ELSE IF ((KEY.EQ."PSOK").OR.(KEY.EQ."PSKI")) THEN 
         NABATY =  -1    ! = 1 for real operators .. -1 for imm. op. 
         ELSE IF((KEY.EQ."FCZK").OR.(KEY.EQ."SDZK")) THEN 
         NABATY = 1      ! = 1 for real operators .. -1 for imm. op. 
         ELSE IF((KEY.EQ."FCBS").OR.(KEY.EQ."SDBS")) THEN 
         NABATY = 1      ! = 1 for real operators .. -1 for imm. op. 
      END IF
C     IF   TRIPLET THEN NABATY = 1 !!!!!!!!!!!!!!!!!!
         NABAOP = 1      ! number of right hand sides. dejarlo asi . solo 1
C
C     Zero the property tensors
cdx      IF (MAGSUS) CALL DZERO(SUSDZD,9)

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
      KJ = 0
      DO 300 IDIP = 1,NLAB1
         LABEL1 = LISTA1(IDIP)
         ISYM = 1   ! ISYM deberia ser 1 ... 
c         ISYM=ISYMAX(IDIP,1)+1
c        set variables for response module
         IF(JIMPRT.GT.3) THEN 
            write(lupri,*) ' about to call LNRVAR' 
            write(lupri,*) ' ----------------------------'
         ENDIF
         CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK2),LWORK2)
C
cm         IF (NFRVAL.GT.0) THEN
cm           write(lupri,*)' # of freq val, NFRVAL:',NFRVAL 

C
C           3. Work Allocations:
C
            KGD1   = KWORK1
            KWRKG1 = KGD1
            LWRKG1 = LWORK - KWRKG1
            KSLV   = KGD1 + 2*NVARPT
            KLAST  = KSLV + 2*NVARPT
            IF (KLAST.GT.LWORK) CALL STOPIT('KLAST GT LWORK on LINEARLR
     $           ',' ',KLAST,LWORK)
            KWRK = KLAST
            LWRK = LWORK - KLAST + 1
cx            write(lupri,*) ' KLAST ',KLAST
cx            write(lupri,*) ' LWORK ',LWORK
cx            write(lupri,*) ' IF KLAST GT LWORK you will get an error '

C
C           Find right hand side for right operator and write to file
C           =========================================================
C
            KSYMOP = ISYM
            TRPLET = .FALSE.
            ANTSYM = 0                   !juan ANTSYM = -1 if triplet !!
c            IF(KEY.EQ.'SDBS') ANTSYM=-1
            IF((KEY.EQ."PSOK").OR.(KEY.EQ."PSKI")) ANTSYM=1  ! indistinto para lkin
C           ANTSYM : matrix symmetry of PRPMO matrix
C          (1: symmetric, -1: antisymmetric, 0: unknown)

            CALL GETGPV(LABEL1,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &           WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
            REWIND LUGDVE
            CALL WRITT(LUGDVE,2*NVARPT,WORK(KWRKG1))
            IF (JIMPRT.GE.3) THEN
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
            DO 200 IPL = 1, NLAB2
C
C              Find label and symmetry of the left side operator
C
               LABEL2 = LISTA2(IPL)
               KSYM   = 1
               KJ = KJ +1 
cb               write(lupri,*) '     else KSYM   = 1      ', KSYM 
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
     &             WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
C
               IF (JIMPRT.GT.3) THEN
                  WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL2
                  CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,2,1,LUPRI)
               ENDIF
C
C                 Form second order properties SNDPRP
C                 ===================================
C
               REWIND LUSOVE
               CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
               IF (JIMPRT.GT.3) THEN
                 WRITE (LUPRI,'(2A)') 'Solution Vector, label: ',LABEL1
                 CALL OUTPUT(WORK(KSLV),1,NVARPT,1,2,NVARPT,2,1,LUPRI)
               ENDIF
C
               SNDPRP = DDOT(2*NVARPT,WORK(KSLV),1,WORK(KGD1),1)
C
               IF (JIMPRT.GE.2) THEN
                  WRITE (LUPRI,'(1A,I2,5A,F20.12)')
     &         '#',KJ,' Second order property for <<',LABEL2,';',LABEL1,
     &            '>> = ',SNDPRP
               ENDIF
C
C              Write properties into the various property matrices
C              ===================================================
C
C         <<Fc,Kin>>
C       ------------------
               IF (KEY.EQ.'FCZK') THEN
                    LRFCZK(1,1)= SNDPRP
                    LRFCZK(2,2)= SNDPRP
                    LRFCZK(3,3)= SNDPRP
               ENDIF
C         <<Sd,Kin>> : (Good luck to find non diag elements)
C       ------------------
c          xx = (001x;xx)+(001y;xy)+(001z;xz) :  1 + 11 + 21 . ok
c          yy = (002x;xy)+(002y;yy)+(002z;zy) : 13 + 32 + 42 . ok
c          zz = (003x;xz)+(003y;zy)+(003z;zz) : 25 + 44 + 54 . ok
               IF (KEY.EQ.'SDZK') THEN
                  IF ((KJ.EQ.1).OR.(KJ.EQ.11).OR.(KJ.EQ.21)) THEN
                     LRSDZK(1,1)=LRSDZK(1,1)+SNDPRP
                  ENDIF
                  IF ((KJ.EQ.13).OR.(KJ.EQ.32).OR.(KJ.EQ.42)) THEN
                     LRSDZK(2,2)=LRSDZK(2,2)+SNDPRP
                  ENDIF
                  IF ((KJ.EQ.25).OR.(KJ.EQ.44).OR.(KJ.EQ.54)) THEN
                     LRSDZK(3,3)=LRSDZK(3,3)+SNDPRP
                  ENDIF
               ENDIF
C         <<Fc,Bso>>
C       ------------------
               IF(KEY.EQ.'FCBS') THEN
                 IF (KJ.EQ.1) LRFCBS(1,1)=SNDPRP
                 IF (KJ.EQ.4) LRFCBS(2,2)=SNDPRP
                 IF (KJ.EQ.6) LRFCBS(3,3)=SNDPRP
               ENDIF
C         <<Sd,Bso>> : (Good luck to find non diag elements)
C       ------------------
c          xx = (001x;xx)+(001y;xy)+(001z;xz) :  1 + 11 + 21 . ok
c          yy = (002x;xy)+(002y;yy)+(002z;zy) : 13 + 32 + 42 . ok
c          zz = (003x;xz)+(003y;zy)+(003z;zz) : 25 + 44 + 54 . ok
               IF (KEY.EQ.'SDBS') THEN
                  IF ((KJ.EQ.1).OR.(KJ.EQ.11).OR.(KJ.EQ.21)) THEN
                     LRSDBS(1,1)=LRSDBS(1,1)+SNDPRP
                  ENDIF
                  IF ((KJ.EQ.13).OR.(KJ.EQ.32).OR.(KJ.EQ.42)) THEN
                     LRSDBS(2,2)=LRSDBS(2,2)+SNDPRP
                  ENDIF
                  IF ((KJ.EQ.25).OR.(KJ.EQ.44).OR.(KJ.EQ.54)) THEN
                     LRSDBS(3,3)=LRSDBS(3,3)+SNDPRP
                  ENDIF
               ENDIF
C         <<L,Lkin>>
C       ------------------
               IF(KEY.EQ.'LKIN') THEN
                 IF (KJ.EQ.1) LRLKIN(1,1)=SNDPRP
                 IF (KJ.EQ.5) LRLKIN(2,2)=SNDPRP
                 IF (KJ.EQ.9) LRLKIN(3,3)=SNDPRP
               ENDIF
C         <<L,Psokin>>
C       ------------------
               IF(KEY.EQ.'PSOK') THEN
                 IF (KJ.EQ.1) LRPSOK(1,1)=SNDPRP
                 IF (KJ.EQ.5) LRPSOK(2,2)=SNDPRP
                 IF (KJ.EQ.9) LRPSOK(3,3)=SNDPRP
               ENDIF
C         <<L,Psokin>>  ! this is for debuigging !!!
C       ------------------
               IF(KEY.EQ.'PSKI') THEN
                 IF (KJ.EQ.1) LRPSKI(1,1)=SNDPRP
                 IF (KJ.EQ.5) LRPSKI(2,2)=SNDPRP
                 IF (KJ.EQ.9) LRPSKI(3,3)=SNDPRP
               ENDIF
C         <<Dia,Dw>>  
C       ------------------
               IF(KEY.EQ."DIAD") THEN
                 IF (KJ.EQ.1) LRDIAD(1,1)=SNDPRP
                 IF (KJ.EQ.5) LRDIAD(2,2)=SNDPRP
                 IF (KJ.EQ.9) LRDIAD(3,3)=SNDPRP
               ENDIF
C         <<Dia,Mv>>  
C       ------------------
               IF(KEY.EQ."DIAM") THEN
                 IF (KJ.EQ.1) LRDIAM(1,1)=SNDPRP
                 IF (KJ.EQ.5) LRDIAM(2,2)=SNDPRP
                 IF (KJ.EQ.9) LRDIAM(3,3)=SNDPRP
               ENDIF
  200 CONTINUE
cm         END IF  NFRVAL
  300 CONTINUE

      CALL GPCLOSE(LUSOVE,'DELETE')
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')
C
C
      write(lupri,*) 
      WRITE(LUPRI,*)
      CALL TIMER ('LRSCLIN',TIMEIN,TIMOUT)
C
      CALL QEXIT('LRSCLIN')
      RETURN
      END
C...
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



