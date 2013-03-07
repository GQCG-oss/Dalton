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
c    ---------------------------------------
c    ---------------------------------------
C
      SUBROUTINE AVELRSC(KEY,WORK,LWORK)
C  
C...
C...  This subroutine was written by Juan Ignacio Melo using
C...  the subroutine ABACTOCD  as a model (2012)
C
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <trkoor.h>
c#include <sigma.h>
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <cbilnr.h>
c#include <suscpt.h>
#include <infpri.h>
      LOGICAL FOUND
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL1,LISTA1(10)
      CHARACTER*6 LABEL2
      CHARACTER*4 KEY
      CHARACTER*3 char
      DIMENSION SIGMAMV(9),SIGMADW(9)
      PARAMETER (D05=0.5D0,D025=0.25)
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
c#include <infsop.h>
#include <ibtfun.h>
C
#include "lrescinf.h"
#include "chrxyz.h"
#include "chrnos.h"
#include "orgcom.h"
C
      IPRLNR = JIMPRT
cxu
C
C
C  LINEAR LRESC SINGLET ROUTINE
C   KEY = 'FCAV'  <Fc>
C         'DIAK'  <Dia.Kin>
C         'ANGP'  <L.Pso>
C
C
      CALL QENTER('AVELRSC')
      CALL TIMER('START ',TIMEIN,TIMOUT)
      WRITE(LUPRI,*)
c      IF (JIMPRT .GT. 5) THEN
         WRITE (LUPRI,'(A,/)') ' <<<<<<<<<< Output from LRESC-'//
     &       'AVERAGE >>>>>>>>>> '
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(2A)') '@AVERAGE to do ', KEY
c      END IF
C
      IPRRSP = JIMPRT
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
      IF (JIMPRT.GT.5) Then 
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
      WORK1  = LFREE
      IF (JIMPRT.GT.5) Then 
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
      ENDIF 

      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO),WORK(KFREE),LFREE)
C          RD_SIRIFC( KEY ,FOUND,   AMAT   ,  WRK      ,LWRK)
 
      IF (.NOT.FOUND) CALL QUIT('AVELRSC error: CMO not on SIRIFC')
      IF (JIMPRT.GT.5) THEN 
         write(lupri,*)' CMOS :  dan bien !!! aunq no parezca.' 
         CALL OUTPUT(WORK(KCMO),1,NORBT,1,5,NORBT,NORBT,1,LUPRI)
         write(lupri,*) '   '     
      ENDIF
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
          write(lupri,*) '    WORK1=LFREE   : ', WORK1
          write(lupri,*) '   '     
      ENDIF
C not needed for mean value
C      CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK1),LWORK1)
C
C     we keep this just in case
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
      LWORK2 = LWORK  - KWORK2
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
C
C  blanquear todas las labels antes 
C
      NLAB1 = 0
      NLAB2 = 0 
      DO I = 1, 10
       LISTA1(I)='        '
      ENDDO
C
C ==============================================================================
C
C  Starting Labels stuff
C
C ==============================================================================

c       write(lupri,*)
cb      write(lupri,*)' selected atom is :', LRATOM
cb      write(lupri,*)' LRATOM / 100' , LRATOM / 100
cb      write(lupri,*)' LRATOM / 10' , LRATOM / 10
 
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
cb
C===============================================================================
C  SINGLET CALCULATIONS : D1S and P1S 
C ==============================================================================
C
C      Look for LABELS : FERMI  
C     ---------------------------
      IF(KEY.EQ."FCAV") THEN 
        nlab=1
        LABEL2='FC '//NAMN(LRATOM)(1:3)
cx        IF (aa.GE.1) THEN
cx          write(char,'(i3)') LRATOM
cx          LISTA1(1) = LABEL2//char
cx        write(lupri,*)' label2 for FC :',LISTA1(1) 
cx          ELSE IF (bb.GE.1) 
        IF (bb.GE.1) THEN
           write(char,'(i2)') LRATOM
           LISTA1(1) = LABEL2//char
        ELSE IF (bb.LT.1) THEN
           write(char,'(i1)') LRATOM
           LISTA1(1) = LABEL2//'0'//char
        ENDIF
      ENDIF
C
C ---------------------------------------------------------------------
C
C      Look for LABELS : DIAKIN
C     ---------------------------
      IF(KEY.EQ."DIAK") THEN
        nlab=3 
        LABEL2='NSKE'
        IJ = 1
        DO I=npos1,npos1+2
           LISTA1(IJ)= CHRNOS(I/100)//CHRNOS(I/10)//
     &                    CHRNOS(MOD(I,10))//'NSKE'//CHRXYZ(IJ)
          IJ = IJ + 1
        END DO
      ENDIF
C
C
C ---------------------------------------------------------------------
C
C   AngPso is done with a unity matrix inside : AngPso = ANGMOM.1.PSO
C
      IF(KEY.EQ."ANGP") THEN
        nlab = 6
C      Look for LABELS : ANGPSO
C     ---------------------------
        LISTA1(1)='XANGMOM'
        LISTA1(2)='YANGMOM'
        LISTA1(3)='ZANGMOM'
C
C      Look for LABELS : PSO 
C     ----------------------       
        LABEL1='PSO'
        IJ = 4
        DO I=npos1,npos1+2
           LISTA1(IJ)= 'PSO '//CHRNOS(I/100)//CHRNOS(I/10)//
     &            CHRNOS(MOD(I,10))
           IJ = IJ + 1
        END DO
C
C ---------------------------------------------------------------------

        write(lupri,*)' NAOS  ' , NAOS(1)
        write(lupri,*)' NASHT ' , NASHT
        write(lupri,*)' NCMOT ' , NCMOT
        write(lupri,*)'   KCMO  ', KCMO
        write(lupri,*)'   NORBT ', NORBT
        write(lupri,*)'   NORB  ', NORB(1)

        write(lupri,*)'NISHT,NSSHT,NOCCT,NBAST,NRHFT,NVIRT :' ,
     &                 NISHT,NSSHT,NOCCT,NBAST,NRHFT,NVIRT
        DO j = 1 , NORBT
           write(lupri,*) ' ',j ,  WORK(KCMO+j)
        ENDDO
        write(lupri,*)' about to call angpso' 
        CALL ANGPSO(LISTA1,WORK(KCMO))
C     This is to avoid the next loop to call PRP1AVE
      NLAB=0
      ENDIF

C  Print LABELS
C -----------------------------------
      write(lupri,*) '@Average setting  LABEL1 :'
      DO i =1, nlab
         WRITE(LUPRI,*)'   LABEL :', LISTA1(I) 
      ENDDO
      WRITE(LUPRI,*)' ---------------------------------- '
C
C ---------------------------------------------------------------------
C
C     Zero the property tensors
cdx      IF (MAGSUS) CALL DZERO(SUSDZD,9)

C   
C        Loop over the right operators which are the
C        the dipole velocity operators 
C        ===========================================
C   
C     
      IDIP = 0 
      DO 300 IDIP = 1,NLAB
         LABEL1 = LISTA1(IDIP)
C
C           3. Work Allocations:
C
         KGD1   = KWORK1
         KWRKG1 = KGD1
         LWRKG1 = LWORK - KWRKG1
         KSLV   = KGD1 + 2*NVARPT
         KLAST  = KSLV + 2*NVARPT
cb         write(lupri,*) ' KLAST ',KLAST
cb         write(lupri,*) ' LWORK ',LWORK
cb         write(lupri,*) ' IF KLAST GT LWORK you will get an error '
         IF (KLAST.GT.LWORK) CALL STOPIT('KLAST GT LWORK on AVELRESC
     $        ',' ',KLAST,LWORK)
         KWRK = KLAST
         LWRK = LWORK - KLAST + 1
         write(lupri,*)' NVARPT :', NVARPT 
         CALL PRP1AVE(LABEL1,SNDPRP,WORK(KCMO),WORK(KUDV),
     $        WORK(KWRKG1),LWRKG1,JIMPRT)
c        CALL PRP1AVE(LBL1AV(IOP),AVEVAL,CMO,UDV,WRK,LWRK,IPRRSP)
cxz      WORK(KCMO),WORK(KUDV),WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
C
        IF (JIMPRT.GE.1) THEN
            WRITE (LUPRI,'(3A,F20.12)')
     &     ' Expectation Value for Operator : < ',LABEL1,
     &     ' > = ',SNDPRP
        ENDIF
C
C       =========================================================
C
cxzC              Nuclear Shieldings
cxzC              ------------------
cxz               IF (KEY.EQ.'D1S ') THEN
cxz                  IF (LABEL2.EQ.'MASSVELO') THEN 
cxz                      SIGMAMV(IDIP)=SNDPRP
cxz                      write(lupri,*)'massvelo' , IDIP
cxz                  ENDIF
cxz                  IF (LABEL2.EQ.'DARWIN  ') THEN 
cxz                     SIGMADW(IDIP)=SNDPRP
cxz                     write(lupri,*)'darwin'  , IDIP
cxz                  ENDIF
cxz               END IF      
cx                  IF (LABEL1.EQ.'XDIPVEL') THEN
cx                     IF (LABEL2(8:8).EQ.'Y') THEN
cx                        SIGMADZ(IPTAX(3,2),INUCLEO)=
cx     &                     SIGMADZ(IPTAX(3,2),INUCLEO)+SNDPRP
cx                     ELSE IF (LABEL2(8:8).EQ.'Z') THEN
cx                        SIGMADZ(IPTAX(2,2),INUCLEO)=
cx     &                     SIGMADZ(IPTAX(2,2),INUCLEO)-SNDPRP
cx                     ENDIF
cx                  ELSE IF (LABEL1.EQ.'ZDIPVEL') THEN
cx                     IF (LABEL2(8:8).EQ.'X') THEN
cx                        SIGMADZ(IPTAX(2,2),INUCLEO)=
cx     &                     SIGMADZ(IPTAX(2,2),INUCLEO)+SNDPRP
cx                     ELSE IF (LABEL2(8:8).EQ.'Y') THEN
cx                        SIGMADZ(IPTAX(1,2),INUCLEO)=
cx     &                     SIGMADZ(IPTAX(1,2),INUCLEO)-SNDPRP
cx                     ENDIF
cx                  ENDIF
cxz  200 CONTINUE
cm         END IF  NFRVAL
  300 CONTINUE
      write(lupri,*)' aca llegue despues de angpso' 


cxz      CALL GPCLOSE(LUREVE,'DELETE')
C   
C     for magnetizability conversion factor -0.25
C   
cb      IF (MAGSUS) THEN
cb        Do I=1,3
cb          Do J=1,3
cb            SUSDZD(I,J)=-D025*SUSDZD(I,J)
cb          ENDDO
cb        ENDDO
cb      ENDIF
cb      IF (SHIELD) THEN
C   
C     for shieldings conversion factor -0.5
C   
cb        DO J=1,3
cb          Do I=1,NCOOR
cb            SIGMADZ(J,I)=-D05*SIGMADZ(J,I)
cb            SIGMASFTP(J,I)=D05*SIGMASFTP(J,I)
cb            SIGMASFTM(J,I)=D05*SIGMASFTM(J,I)
cb          ENDDO
cb        ENDDO
cb      ENDIF
C
c      WRITE(LUPRI,*) ' llegue al final de AVERAGE.' 
      CALL TIMER ('AVELRSC',TIMEIN,TIMOUT)
C
      CALL QEXIT('AVELRSC')
      RETURN
      END
C...
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
