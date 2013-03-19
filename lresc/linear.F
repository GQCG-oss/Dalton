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
<<<<<<< HEAD
c    ---------------------------------------
C  /* Deck lrscdrv */
      SUBROUTINE LRSCDRV(WORK,LWORK)
C
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <lrescinf.h>
C

      write(lupri,*)'@ LRSCDRV : ' 
      IF (.NOT. LRESCALL) THEN 
         WRITE (LUPRI,'(A)')  '   all logic  '
      END IF

      IF (SIGMAP1S)  THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA1S for : ', WORD
          CALL LRSCLIN('P1S',WORK,LWORK)
      END IF 

      IF (SIGMAP1T) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA1T for : ', WORD
      END IF

      IF (SIGMAP3S) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA3S for : ', WORD
      END IF

      IF (SIGMAP3T) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA3T for : ', WORD
      END IF

      IF (SIGMAD0S) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA0S for : ', WORD
      END IF

      IF (SIGMAD1S) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA1S for : ', WORD
         CALL LRSCLIN('D1S',WORK,LWORK)
      END IF

      RETURN
      END
C...
c    ---------------------------------------
c    ---------------------------------------

      SUBROUTINE LRSCLIN(KEY,WORK,LWORK)
C based on CTOABA 
C...
C...  This subroutine was written by Juan Ignacio Melo using
C...  the subroutine ABACTOCD  as a model (2012)
=======
      SUBROUTINE LRSCLIN(WORK,LWORK)
C based on CTOABA 
C...
C...  This subroutine was written by Andrea Ligabue using
C...  the subroutine LNRABA as a model (1999)
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
C...
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <trkoor.h>
<<<<<<< HEAD
c#include <sigma.h>
=======
#include <sigma.h>
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <cbilnr.h>
<<<<<<< HEAD
c#include <suscpt.h>
#include <infpri.h>
      LOGICAL CICLC, HFCLC, TRIPLE, EXECLC, FOUND
      DIMENSION WORK(LWORK)
      CHARACTER*3  KEY
      CHARACTER*8 LABEL1,LABEL2,LISTA1(4*MXCOOR+9),LISTA2(4*MXCOOR+9)
=======
#include <suscpt.h>
#include <infpri.h>
      LOGICAL CICLC, HFCLC, TRIPLE, EXECLC, FOUND
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL1, LABEL2, LABINT(4*MXCOOR+9)
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
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
<<<<<<< HEAD
#include "chrxyz.h"
#include "chrnos.h"
#include "orgcom.h"
C
cxu  sacar lo de abajo
      IPRLNR =2
=======
cxu  sacar lo de abajo
      IPRLNR = 4
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
cxu
      CALL QENTER('LRSCLIN')
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPRLNR .GT. 0) WRITE (LUPRI,'(A,/)')
<<<<<<< HEAD
     *    ' <<<<<<<<<< Output from CTOABA jim >>>>>>>>>> '
      WRITE(LUPRI,'(A)') 'Called to do ', KEY
=======
     *    ' <<<<<<<<<< Output from linear-lresc >>>>>>>>>> '
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
C
      IPRRSP = IPRLNR
C
C     Get reference state
C     ===================
C
C     1. Work Allocations:
C
<<<<<<< HEAD
      LUDV   = N2ASHX
      LPVX   = 0
      KFREE  = 1
      LFREE  = LWORK
      write(lupri,*) ' from orgcom.h : GAGORG '
      write(lupri,*) GAGORG
      write(lupri,*) ' from orgcom.h : ORIGIN '
      write(lupri,*) ORIGIN
      write(lupri,*) ' from orgcom.h : CMXYZ  '
      write(lupri,*) CMXYZ
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

=======
cdx      IF (ABASOP) THEN
cdx        LUDV   = NORBT * NORBT
cdx        LPVX   = LPVMAT
cdx      ELSE
      write(lupri,*)'jim : N2ASHX ;' , N2ASHX
      LUDV   = N2ASHX
      LPVX   = 0
cdx      ENDIF
      KFREE  = 1
      LFREE  = LWORK
C      
      write(lupri,*) 'NASHT INTEGER common inforb :', NASHT
      write(lupri,*) '      Total Number of Active orbitals, '
      write(lupri,*) '      deberia ser ZERO para este caso '

      write(lupri,*)'   LABINT (4*MXCOOR+9) = ', 4*MXCOOR+9
      write(lupri,*)'   MXCOOR              = ', MXCOOR
      write(lupri,*)'   NCMOT = NORB * NORB = ', NCMOT
      write(lupri,*)'   KCMO  no se que sea = ', KCMO
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KUDV  ,LUDV  ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KPVX  ,LPVX  ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KXINDX,LCINDX,WORK ,KFREE ,LFREE)
<<<<<<< HEAD
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
=======
C
      KWORK1 = KFREE
      LWORK1 = LFREE
      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO),WORK(KFREE),LFREE)
ccc        RD_SIRIFC( KEY ,FOUND,   AMAT   ,  WRK      ,LWRK)
      write(lupri,*)' jaja    los cmo dan bien !!!' 
c.      CALL OUTPUT(WORK(KCMO),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,LUPRI)
      IF (.NOT.FOUND) CALL QUIT('LRSCLIN error: CMO not on SIRIFC')
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c

cx ACA es para alguna capa activa
cx      IF (NASHT .GT. 0) THEN
cx         CALL RD_SIRIFC('DV',FOUND,WORK(KWORK1),WORK(KWORK1),LWORK1)
cx         write(lupri,*)'jim  DV found on RD_SIFC '
cx         IF (.NOT.FOUND)
<<<<<<< HEAD
cx     &      CALL QUIT('ROUTINE error: DV not found on SIRIFC')
=======
cx     &      CALL QUIT('DIAMV error: DV not found on SIRIFC')
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
cx         CALL DSPTSI(NASHT,WORK(KWORK1),WORK(KUDV))
cx      END IF
C
      ISYM = 1
<<<<<<< HEAD
      write(lupri,*) '   ' 
      write(lupri,*) ' about to call LNRVAR' 
      write(lupri,*) ' ----------------------------'
      write(lupri,*) '           ISYM   : ', ISYM 
      write(lupri,*) '   KWORK1=KFREE   : ', KWORK1
      write(lupri,*) '    WORK1=LFREE   : ', WORK1
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
=======
      write(lupri,*)' about to call LNRVAR' 
      write(lupri,*)'       ISYM ', ISYM 
      write(lupri,*)'       IPRLNR ', IPRLNR
      write(lupri,*)'       KWORK1 ', KWORK1
      write(lupri,*)'       LWORK1 ', LWORK1
      CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK1),LWORK1)
C
cc      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LWORK1,0)
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
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
<<<<<<< HEAD
cdx            CALL OUTPUT(WORK(KUDV),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,
cdx     &                  LUPRI)
=======
            write(lupri,*)' jeje WORK(KUDV) ni idea q sea' 
            CALL OUTPUT(WORK(KUDV),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,
     &                  LUPRI)
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
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
<<<<<<< HEAD
      write(lupri,*) '        '
      write(lupri,*) ' sigo alocando : '
      write(lupri,*) ' ----------------------------'
      write(lupri,*) ' KIDSYM = KWORK1           : ' , KIDSYM
      write(lupri,*) ' KIDADR = KIDSYM + 9MXCENT : ' , KIDADR
      write(lupri,*) ' KWORK2 = KIDADR + 9MXCENT : ' , KWORK2
      write(lupri,*) ' LWORK2 = LWORK - KWORK2   : ' , LWORK2
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

C
C  Look for NSCO LABELS 
C
      write(lupri,*)
      LABEL1='NSCO'
      write(lupri,*) '@Linear setting  LABEL1 : NSCO '
      IJ = 1
      DO I=1,NCOOR
        DO J=1,3
          LISTA1(IJ)= CHRNOS(I/10)//CHRNOS(I/10)//CHRNOS(MOD(I,10))//
     &     'NSCO'//CHRXYZ(J)
          LABSYM(IJ)=1
          write(lupri,*) '   LABEL1 : ', LISTA1(IJ), LABSYM(IJ)
          IJ = IJ + 1
        END DO
      END DO
      NLAB1 = 3*NCOOR
C
C  Look for DARWIN and MASSVELO LABELS 
C
      write(lupri,*)
      LISTA2(1) ='DARWIN  '
      LISTA2(2) ='MASSVELO'
      NLAB2 = 2
      write(lupri,*) '@Linear setting  LABEL2 :'
      write(lupri,*) '   LABEL2 : ', LISTA2(1)
      write(lupri,*) '   LABEL2 : ', LISTA2(2)

C
C     Set variables for ABARSP and logicals
C
      CICLC  = .FALSE.   ! TRUE for CI calculations
      HFCLC  = NASHT .LE. 1 ! .T. RHF-closed shell or 1e in one active orbital
      TRIPLE = .FALSE. ! .T. for triplet perturbation operators
      EXECLC = .FALSE. ! false for linear response equations
      NABATY = 1 ! = 1 for real operators .. -1 for imm. op.
      NABAOP = 1  ! number of right hand sides. dejarlo asi . solo 1
C
C     Zero the property tensors
C   por ahora solo un numero sale de aca
cdx      IF (MAGSUS) CALL DZERO(SUSDZD,9)

=======
C
      NLBTOT = 0
      NLBSHIS = 4
C
      NCOMP  = 0
      NPATOM = 0
cLig  <> added the TODOINT to see if the property was already in the file
c busco DAREWI
      
cdbg       write(lupri,*)' about to call TODOINT, q no anda...pa esto' 
cdbg trabajar con esto de todoint, por ahora supongo que las int estan 
cdbg en luprop

cf      IF (TODOINT('DARWIN ',LUPROP)) THEN
cf      write(lupri,*) ' to do INT dio TRUE'
        CALL GETLAB('ANGMOM  ',8,LABINT,WORK(KIDSYM),LUPROP)
        write(lupri,*)' after getlab, ANGMOM '  
        write(lupri,*)'       WORK(KIDSYM)', WORK(KIDSYM)
        write(lupri,*)'       LABINT(1)', LABINT(1)
        write(lupri,*)'       LABINT(2)', LABINT(2)  
        write(lupri,*)'       LABINT(3)', LABINT(3)  
cf      ELSE
        
cd        CALL GET1IN(DUMMY,'DARWIN ',NCOMP,WORK(KWORK2),LWORK2,
cd     &              LABINT,WORK(KIDSYM),WORK(KIDADR),
cd     &              IDUMMY,.TRUE.,NPATOM,.TRUE.,IPR1IN)
        
cf      ENDIF
      NLAB = 3
      CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      write(lupri,*)'after labcop ANGMOM,LABSYM NI idea q tiene'
C
C
cf      IF (TODOINT('NSTCGO ',LUPROP)) THEN
        CALL GETLAB('PSO',3,LABINT,WORK(KIDSYM),LUPROP)
        write(lupri,*)' after getlab, PSO  '  
        write(lupri,*)'       WORK(KIDSYM)', WORK(KIDSYM)
        write(lupri,*)'       LABINT(1)', LABINT(1)
        write(lupri,*)'       LABINT(2)', LABINT(2)  
        write(lupri,*)'       LABINT(3)', LABINT(3)  
cf       ELSE
cf         CALL GET1IN(DUMMY,'NSTCGO ',NCOMP,WORK(KWORK2),LWORK2,
cf     &              LABINT,WORK(KIDSYM),WORK(KIDADR),
cf     &              IDUMMY,.TRUE.,NPATOM,.TRUE.,IPR1IN)
cf       ENDIF
       NLAB = 3*NCOOR
       CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
       write(lupri,*)'after labcop PSO,LABSYM NI idea q tiene'



cdx      IF (MAGSUS) THEN
cdx         CALL GET1IN(DUMMY,'RANGMO ',NCOMP,WORK(KWORK2),LWORK2,
cdx     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
cdx     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,IPR1IN)
cdx         NLAB = 9
cdx         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
cdx         NLBSHIS = 13
cdx      ENDIF
cLig  <> added the TODOINT to see if the property was already in the file

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
cdx         NLBSHIE = NLBTOT
cdx      IF (TODOINT('PSO     ',LUPROP)) THEN
cdx           write(lupri,*)'todoint PSO dio FALSE  '
cdx           CALL GETLAB('PSO',3,LABINT,WORK(KIDSYM),LUPROP)
cdx         ELSE
cdx           write(lupri,*)'debtro dle ELSE GEt1in  PSO  '
cdx           CALL GET1IN(DUMMY,'PSO    ',NCOMP,WORK(KWORK2),LWORK2,
cdx     &                 LABINT,WORK(KIDSYM),WORK(KIDADR),
cdx     &                 IDUMMY,.TRUE.,NPATOM,.TRUE.,IPR1IN)
cdx           write(lupri,*)'volvi de get1in  '
cdx         ENDIF
cdx         NLAB = NCOOR
cdx         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
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

>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
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
<<<<<<< HEAD
      NLAB1 = 9
      DO 300 IDIP = 1,NLAB1
         LABEL1 = LISTA1(IDIP)
C
         ISYM=ISYMAX(IDIP,1)+1
         write(lupri,*)' isimax :', ISYM   ! ISYM =1  deberia ser 1        
C     
c        set variables for response module
         write(lupri,*) ' about to call LNRVAR' 
         write(lupri,*) ' ----------------------------'
         CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK2),LWORK2)
C
cm         IF (NFRVAL.GT.0) THEN
           write(lupri,*)' # of freq val, NFRVAL:',NFRVAL 
=======
cxv      DO 300 IDIP = 1,3
         IDIP = 1
cv         LABEL1 = 'DARWIN  ' 
         LABEL1 = '001NSCOX' 
C
cx         ISYM=ISYMAX(IDIP,1)+1
C
         ISYM =1          
         write(lupri,*)' ISYM :',ISYM
c  set variables for response module
         write(lupri,*)' about to call LNRVAR' 
         IPRLNR = 4
         write(lupri,*)' IPRLNR = ', IPRLNR 
         CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK2),LWORK2)
C
         IF (NFRVAL.GT.0) THEN
         write(lupri,*)' antes NFRVAL:',NFRVAL 
         write(lupri,*)'       common de cbilnr.h ' 
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
C
C           3. Work Allocations:
C
            KGD1   = KWORK1
            KWRKG1 = KGD1
            LWRKG1 = LWORK - KWRKG1
            KSLV   = KGD1 + 2*NVARPT
            KLAST  = KSLV + 2*NVARPT
<<<<<<< HEAD
            write(lupri,*) ' KLAST ',KLAST
            write(lupri,*) ' LWORK ',LWORK
            write(lupri,*) ' KLAST GT LWORK ! error '
            IF (KLAST.GT.LWORK) CALL STOPIT('LINEARjim',' ',KLAST,LWORK)
=======
            IF (KLAST.GT.LWORK) CALL STOPIT('CTOABjim',' ',KLAST,LWORK)
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
            KWRK = KLAST
            LWRK = LWORK - KLAST + 1
            write(lupri,*)' NVARPT :', NVARPT 

C
C           Find right hand side for right operator and write to file
C           =========================================================
C
            KSYMOP = ISYM
            TRPLET = .FALSE.
<<<<<<< HEAD
            IPRSP = -2 
            ANTSYM = 0
C   
            write(lupri,*) ' about to call GETPV'  
C jim sym to be read on luprop
C ANTSYM : matrix symmetry of PRPMO matrix
C          (1: symmetric, -1: antisymmetric, 0: unknown)

            CALL GETGPV(LABEL1,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &           WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
=======
C   
            write(lupri,*) ' about to call GETPV'  
            write(lupri,*) ' ANTSYM', ANTSYM  

            CALL GETGPV(LABEL1,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &           WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
cx                 WORD = LABEL1
cx                 FC = DUMMY
cx                 FV = DUMMY
cx                 CMO = WORK(KCMO)
cx                 UDV = WORK(KUDV)
cx                 PV = WORK(KPVX)    --- 4 SOPPA
cx                 XINDX = WORK(KXINDX)
cx                 ANTSYM_GP = ANTSYM
cx                 WRK = WORK(KWRKG1)
cx                 LWRK = LWRKG1
cx    DIMENSION FC(*),CMO(*),UDV(*),PV(*),XINDX(*),WRK(LWRK)
      write(lupri,*) 'KUDV ' , WORK(KUDV)

cx      CALL PRPCTL(WORD,WRK(KGP),CMO,UDV,PV,XINDX,ANTSYM_GP,
cx     &               WRK(KWRK1),LWRK1)

>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
            REWIND LUGDVE
            CALL WRITT(LUGDVE,2*NVARPT,WORK(KWRKG1))
cvx            IF (IPRLNR.GT.-1) THEN
               WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL1
               CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,2,1,LUPRI)
cvx            ENDIF
<<<<<<< HEAD
=======
         write(lupri,*)' antes del lio.?' 
         write(lupri,*)' jim hasta aca. ' 
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
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
<<<<<<< HEAD
C 
            NLAB2 = 2
            DO 200 IPL = 1, NLAB2
C
C              Find label and symmetry of the left side operator
C
               LABEL2 = LISTA2(IPL)
               write(lupri,*)' LABAPP IPL ' , LABEL2
cx             KSYM   = LABSYM(IPL)
               write(lupri,*) '  KSYM   = LABSYM(IPL) ', KSYM 
               KSYM   = 1
               write(lupri,*) '  else KSYM   = 1', KSYM 
=======
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

>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
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
<<<<<<< HEAD
                     WRITE (LUPRI,'(5A,F15.12)')
     &               ' Second order property for ',LABEL2,' , ',LABEL1,
=======
                     WRITE (LUPRI,'(4A,F15.8)')
     &               ' Second order property for ',LABEL2,LABEL1,
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
     &               ' = ',SNDPRP
                  ENDIF
C
C                 Write properties into the various property matrices
C                 ===================================================
C
<<<<<<< HEAD
C                 Nuclear Shieldings
C                 ------------------
cx                  IF ((IPRLBL.GE.NLBSHIS).AND.
cx     &               (IPRLBL.LE.NLBSHIE))  THEN
=======
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
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
cx                     INUCLEO = (IPRLBL-NLBSHIS)/3+1
cx                     IF (LABEL1.EQ.'XDIPVEL') THEN
cx                        IF (LABEL2(8:8).EQ.'Y') THEN
cx                           SIGMADZ(IPTAX(3,2),INUCLEO)=
cx     &                        SIGMADZ(IPTAX(3,2),INUCLEO)+SNDPRP
cx                        ELSE IF (LABEL2(8:8).EQ.'Z') THEN
cx                           SIGMADZ(IPTAX(2,2),INUCLEO)=
cx     &                        SIGMADZ(IPTAX(2,2),INUCLEO)-SNDPRP
cx                        ENDIF
<<<<<<< HEAD
=======
cx                     ELSE IF (LABEL1.EQ.'YDIPVEL') THEN
cx                        IF (LABEL2(8:8).EQ.'Z') THEN
cx                           SIGMADZ(IPTAX(1,2),INUCLEO)=
cx     &                        SIGMADZ(IPTAX(1,2),INUCLEO)+SNDPRP
cx                        ELSE IF (LABEL2(8:8).EQ.'X') THEN
cx                           SIGMADZ(IPTAX(3,2),INUCLEO)=
cx     &                        SIGMADZ(IPTAX(3,2),INUCLEO)-SNDPRP
cx                        ENDIF
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
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
<<<<<<< HEAD
cx               ENDIF
  200 CONTINUE
cm         END IF  NFRVAL
  300 CONTINUE
=======
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
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
      CALL GPCLOSE(LUSOVE,'DELETE')
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')
C   
C     for magnetizability conversion factor -0.25
C   
      IF (MAGSUS) THEN
<<<<<<< HEAD
cb        Do I=1,3
cb          Do J=1,3
cb            SUSDZD(I,J)=-D025*SUSDZD(I,J)
cb          ENDDO
cb        ENDDO
=======
        Do I=1,3
          Do J=1,3
            SUSDZD(I,J)=-D025*SUSDZD(I,J)
cDEBUG	
c	    write(LUPRI,*) 'Dentro ABACTODC'
c	    write(LUPRI,*) 'SUSDZD= ',SUSDZD
cDEBUG
          ENDDO
        ENDDO
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
      ENDIF
      IF (SHIELD) THEN
C   
C     for shieldings conversion factor -0.5
C   
<<<<<<< HEAD
cb        DO J=1,3
cb          Do I=1,NCOOR
cb            SIGMADZ(J,I)=-D05*SIGMADZ(J,I)
cb            SIGMASFTP(J,I)=D05*SIGMASFTP(J,I)
cb            SIGMASFTM(J,I)=D05*SIGMASFTM(J,I)
cb          ENDDO
cb        ENDDO
=======
        DO J=1,3
          Do I=1,NCOOR
            SIGMADZ(J,I)=-D05*SIGMADZ(J,I)
            SIGMASFTP(J,I)=D05*SIGMASFTP(J,I)
            SIGMASFTM(J,I)=D05*SIGMASFTM(J,I)
          ENDDO
        ENDDO
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
      ENDIF
C
      CALL TIMER ('LRSCLIN',TIMEIN,TIMOUT)
C
      CALL QEXIT('LRSCLIN')
      RETURN
      END
C...
c    ---------------------------------------
<<<<<<< HEAD
=======
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
>>>>>>> 46d1f194ad6db75a493ae3805b16ad5c830dd66c
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
