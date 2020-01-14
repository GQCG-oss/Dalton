!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
 
C  FILE : abacus/abaloc.F

C  /* Deck locini */
      SUBROUTINE LOCINI
C
C     Initialize LOCINF
C
#include "implicit.h"
#include "maxorb.h"
C
#include "locinf.h"
C
      CALL QENTER('LOCINI')
C
      FOSOCC = .FALSE.
      FBOVIR = .FALSE.
      FBSETV = .FALSE.
      FBSTVO = .FALSE.
      MULLOC = .FALSE.
      FBOCIN = .FALSE.
      FBOOCC = .FALSE.
      LABOCC = .FALSE.
      LABVIR = .FALSE.
C
      CALL QEXIT('LOCINI')
C
      RETURN
      END
C***********************************************************************
C  /* Deck locinp */
       SUBROUTINE LOCINP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
#include "maxorb.h"
C
       PARAMETER (NTABLE = 10, D0 = 0.0D0, NLAB = 20)
       PARAMETER (IZERO = 0)
       INTEGER NAUXV, NAUXO
C
       LOGICAL NEWDEF
       CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include "locinf.h"
#include "inforb.h"
C
CSPAS:31.05.11: MULLOC is not yet working correctly
C      DATA TABLE /'.FOSBOY','.FBOVIR','.FBSETV','.FBSTVO','.MULLOC',
       DATA TABLE /'.FOSBOY','.FBOVIR','.FBSETV','.FBSTVO','.XXXXXX',
CKeinSPASmehr
     &             '.FBOCIN','.FBOOCC','.LABOCC','.LABVIR','.XXXXXX'/
C
       CALL QENTER('LOCINP')
C
       NEWDEF = (WORD .EQ. '*LOCALI')
       ICHANG = 0
       NO2LOC = 0
       NV2LOC = 0
       DO I = 1, NLAB
         TABOCL(I) = '-**--**-'
         TABVIL(I) = '-**--**-'
         NTVI2L(I) = IZERO
         NTOC2L(I) = IZERO
       END DO
C
       IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
            READ (LUCMD,'(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GOTO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                    GOTO (1,2,3,4,5,6,7,8,9,10), I
                  END IF
 200           CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GOTO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "', WORD,
     &               '" not recognized in LOCALI.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Input keywords in LOCALI')
 1             CONTINUE
                  FOSOCC = .TRUE.
               GOTO 100
 2             CONTINUE
CPFP Jan/28/05
                  IF (FOSOCC .OR. FBOCIN .OR. FBOOCC) THEN
C
C                   FOSOCC = .TRUE.
                    FBOVIR = .TRUE.
C
                    IF(NO2LOC .EQ. IZERO) THEN
                      NOCAUX = NOCCT
                    ELSE
                      NOCAUX = NO2LOC
                    END IF
C
                    AUXV  = NVIRT / NOCAUX
                    NAUXV = INT(AUXV)
                    NCHSV = NAUXV*NOCAUX
                    IF (NCHSV .EQ. NVIRT ) THEN
                       NVSET = NAUXV
                    ELSE
                       NVSET = NAUXV + 1
                    END IF
C
                  ELSE
C
                     WRITE (LUPRI,'(/,2X,A,A,/)') 
     &              ' .LABOCC keyword needs specification of .FOSOCC,',
     &              ' .FBOOCC, .MULLOC  or .FBOCIN '
C
                  END IF
Cend-PFP Jan/28/05
C
               GOTO 100
 3             CONTINUE
CPFP Jan/28/05
                  IF (FOSOCC .OR. FBOCIN .OR. FBOOCC) THEN
C
C                   FOSOCC = .TRUE.
                    FBOVIR = .FALSE.
                    FBSETV = .TRUE.
                    READ (LUCMD,*) NFBSET
                    NAUXV = INT(NFBSET)
                    NCHSV = NAUXV*NOCCT
                    NVSET = NAUXV
C
                    IF(NO2LOC .EQ. IZERO) THEN
                      NOCAUX = NOCCT
                    ELSE
                      NOCAUX = NO2LOC
                    END IF
C
                    IF (NCHSV .GT. NVIRT) THEN
                      AUXV  = NVIRT / NOCAUX
                      NAUXV = INT(AUXV)
                      NCHSV = NAUXV*NOCAUX
                      IF (NCHSV .EQ. NVIRT ) THEN
                        NVSET = NAUXV
                      ELSE
                        NVSET = NAUXV + 1
                      END IF
CPFP Jan/28/05
C                     FOSOCC = .TRUE.
Cend-PFP Jan/28/05
                      FBOVIR = .TRUE.
                      FBSETV = .FALSE.                
C
                      WRITE (LUPRI,'(/,A,I4,A,I4,/,A,I4,/)') 
     &                  ' Number of virtual orbitals chosen is',NVSET,
     &                  'times', NOCCT,
     &                  ', which is larger than the total',NVIRT
                       WRITE (LUPRI,'(/,A,/)')
     &                  ' Therefore all virtual orbitals are localized'
C
                    END IF
C
                  ELSE
C
                     WRITE (LUPRI,'(/,2X,A,A,/)') 
     &              ' .LABOCC keyword needs specification of .FOSOCC,',
     &              ' .FBOOCC, .MULLOC  or .FBOCIN '
C
                  END IF
Cend-PFP Jan/28/05
C
               GOTO 100
C
 4             CONTINUE
CPFP Jan/28/05
                  IF ( FOSOCC .OR. FBOCIN .OR. FBOOCC ) THEN
C
C                   FOSOCC = .TRUE.
                    FBOVIR = .FALSE.
                    FBSETV = .FALSE.
                    FBSTVO = .TRUE.
                    READ (LUCMD,*) NFBSET,NV2LOC
C
                    IF(NO2LOC .EQ. IZERO) THEN
                      NOCAUX = NOCCT
                    ELSE
                      NOCAUX = NO2LOC
                    END IF
C
                    IF (NV2LOC .GT. NOCAUX) NV2LOC = NOCAUX
                    READ (LUCMD,*) (NTVI2L(I), I = 1,NV2LOC)
                    NAUXV = INT(NFBSET)
                    NCHSV = NAUXV*NV2LOC
                    NVSET = NAUXV
C
                    IF (NCHSV .GT. NVIRT) THEN
                       AUXV  = NVIRT / NV2LOC
                       NAUXV = INT(AUXV)
                       NCHSV = NAUXV*NV2LOC
                       IF (NCHSV .EQ. NVIRT ) THEN
                         NVSET = NAUXV
                       ELSE
                         NVSET = NAUXV + 1
                       END IF
C
                       WRITE (LUPRI,'(/,2x,A,I4,A,I4,/,A,I4,/,A,/)') 
     &                  ' Number of virtual orbitals chosen is',NVSET,
     &                  'times', NV2LOC,
     &                  ', which is larger than the total',NVIRT,
     &                  'Therefore all virtual orbitals are localized'
C
                    END IF
C
                  ELSE
C
                     WRITE (LUPRI,'(/,2X,A,A,/)') 
     &              ' .LASTV keyword needs specification of .FOSOCC,',
     &              ' .FBOOCC, .MULLOC  or .FBOCIN '
C
                  END IF
Cend-PFP Jan/28/05
C
               GOTO 100
C
 5             CONTINUE
CSPAS:31.05.11: MULLOC is not yet working correctly
C                 MULLOC = .TRUE.
CKeinSPASmehr
               GOTO 100
 6             CONTINUE
                  FBOCIN = .TRUE.
                  READ (LUCMD,*) NAUXO
                  NO2LOC = INT(NAUXO)
                  IF (NAUXO .GT. NOCCT) NO2LOC  = NOCCT
                  READ (LUCMD,*) (NTOC2L(I),I=1,NO2LOC)
C
               GOTO 100
CPFP  Be Crefull here, because this option is for occupied orbitals to be not localized
 7             CONTINUE
                  FBOOCC = .TRUE.
                  READ (LUCMD,*) NAUXO
CPFP
C                  NO2LOC = INT(NAUXO)
                  NO2LOC = INT(NOCCT-NAUXO)
                  IF (NAUXO .GT. NOCCT) NO2LOC  = NOCCT
C                  READ (LUCMD,*) (NTOC2L(I),I=1,NO2LOC)
                  READ (LUCMD,*) (NTOC2L(I),I=1,NAUXO)
Cend-PFP
C
               GOTO 100
 8             CONTINUE
                  LABOCC = .FALSE.
                  READ (LUCMD,*) NOCLAB
                  IF (FOSOCC .OR. FBOOCC .OR. MULLOC .OR. FBOCIN) THEN
C
                    LABOCC = .TRUE.
                    READ (LUCMD,*) (TABOCL(I),I=1,NOCLAB)
C
                  ELSE
C
                     WRITE (LUPRI,'(/,2X,A,A,/,2X,A,/)') 
     &              ' .LABOCC keyword needs specification of .FOSOCC,',
     &              ' .FBOOCC, .MULLOC  or .FBOCIN ',
     &              ' and only allows 20 labels of 8 characters each'
C
                  END IF
C
               GOTO 100
 9             CONTINUE
                  LABVIR = .FALSE.
                  READ (LUCMD,*) NVILAB
                  IF (FBSTVO .OR. FBOVIR .OR. FBSETV) THEN
C
                    LABVIR = .TRUE.
                    READ (LUCMD,*) (TABVIL(J),J=1,NVILAB)
C
                   ELSE
C
                     WRITE (LUPRI,'(/,2X,A,A,/,2X,A,/)') 
     &              ' .LABVIR keyword needs specification of .FBOVIR,',
     &              ' .FBSETV or .FBSTVO ',
     &              ' and only allows 20 labels of 8 characters'
C
                  END IF
C
               GOTO 100
 10            CONTINUE
               GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GOTO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &               '" not recognized in LOCALI'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT ('Illegal keyword in LOCALI')
            END IF
       END IF
 300   CONTINUE
       IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for LOCALI:',0)
C
         IF (FOSOCC) WRITE (LUPRI,'(A,A)')
     &      ' FOSOCC: Foster Boys localization of all',
     &      ' occupied orbitals'
C
         IF (MULLOC) WRITE (LUPRI,'(A,A)')
     &   ' MULLOCM: ulliken population localization of',
     &   ' all occupied orbitals'
C
         IF (FBOVIR) WRITE (LUPRI,'(/,A,A,I4,A)')
     &      ' FBOVIR: Localization of all virtual orbital',
     &      ' is calculated', NVSET, ' sets '
C
         IF (FBSETV) WRITE (LUPRI,'(/,A,I4,A)')
     &      ' FBSETV: Localization of  ',NVSET,
     &      '  virtual orbital(s) per occupied orbital '
C
         IF (FBSTVO) THEN
            WRITE (LUPRI,'(/,A,I4,A,/,I4,A)')
     &      ' FBSTVO: Localization of  ',NVSET,
     &      '  set(s) of virtual orbital(s) respect to',NV2LOC,
     &      '  occupied orbital, which are :'
            WRITE(LUPRI,*) (NTVI2L(I), I = 1, NV2LOC)
         END IF
C
CPFP Nov/21/04
         IF (FBOCIN) THEN
            WRITE (LUPRI,'(/,2X,A,I5,A,/)')
     &      ' FBOCIN: Localize all occ. and chose ',NAUXO,
     &      '  occupied orbital(s) [to be delocalized], which are :'
            WRITE(LUPRI,*)  (NTOC2L(I),I=1,NAUXO)
         END IF
C
         IF (FBOOCC) THEN
            WRITE (LUPRI,'(/,2X,A,I5,A,/)')
     &      ' FBOOCC: Localization of ',NOCCT-NAUXO,
     &      '  occupied orbital(s). Non-localization for :'
            WRITE(LUPRI,*)  (NTOC2L(I),I=1,NAUXO)
         END IF
C
         IF (LABOCC) THEN
            WRITE (LUPRI,'(/,2X,A,A,/)')
     &      ' LABOCC: Labels are added to localized occupied',
     &      ' orbitals, they are:'
            WRITE(LUPRI,*)  NOCLAB
            WRITE(LUPRI,*)  (TABOCL(J), J=1,NOCLAB)
         END IF
C
         IF (LABVIR) THEN
            WRITE (LUPRI,'(/,2X,A,A,/)')
     &      ' LABVIR: Labels are added to localized virtual',
     &      ' orbitals, they are:'
            WRITE(LUPRI,*)  NVILAB
            WRITE(LUPRI,*)  (TABVIL(J),J=1,NVILAB)
         END IF
Cend-PFP Nov/21/04
C
       END IF
C
       CALL QEXIT('LOCINP')
C
       RETURN
       END
C***********************************************************************
C  /* Deck getloc */
      SUBROUTINE GETLOC(IPRSOS,NLOCC,NLVIR,CLOCO,CLOCV,PRPMO,
     &                  WORK,LWORK)
C
C      Stephan P. A. Sauer 18/11-1998
C
C      This routine is 
C
#include "implicit.h"
C
C 
C IRAT used from include file iratdef.h
C LUPRI used from include file priunit.h
C MXCENT is needed in COMMON /
C inforb for NBAST
#include "iratdef.h"
#include "priunit.h"
#include "maxash.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
C NBAST, NNBASX, N2BASX, NORBT, NSYM, NOCCT, NVIRT,
C NRHF(8), NVIR(8), NBAS(8), NORB(8) used from COMMON /INFORB/
C NWOPPT used from COMMON /INFLIN/
C LUPROP used from COMMON /INFTAP/
#include "inforb.h"
#include "inflin.h"
#include "inftap.h"
#include "locinf.h"
#include "infpri.h"
#include "molde.h"
C
      PARAMETER (D0 = 0.0D+00, D1 = 1.0D+00, D2 = 2.0D+00)
      PARAMETER (THRLOC = 1.0D-12)
      PARAMETER (MXCGTO = 8)
C
      DIMENSION   CLOCO(NOCCT,NLOCC)
      DIMENSION   CLOCV(NVIRT,NLVIR)
      DIMENSION   PRPMO(NORBT,NORBT,3)
      DIMENSION   WORK(LWORK)
C
      LOGICAL     FNDLB2
      DIMENSION   LABSYM(3)
      CHARACTER*8 LABEL, LABAPP(3), RTNLBL(2)
C
      DIMENSION   OCCUP(NORBT)
C
C-----------------------------------------------------------
C     Initialize transformation matrix to localized orbitals
C-----------------------------------------------------------
C
      CALL QENTER('GETLOC')
C
      CALL DZERO(CLOCO,NOCCT*NLOCC)
      CALL DZERO(CLOCV,NVIRT*NLVIR)
C
C---------------------------------
C     This is only for close shell
C---------------------------------
C
      DO INDEX = 1, NBAST
         IF ( INDEX .LE. NOCCT ) THEN 
            OCCUP(INDEX) = D2
         ELSE
            OCCUP(INDEX) = D0
         END IF
      END DO
C
C---------------------------------
C     This is only for close shell
C---------------------------------
C
      DO IOCCT = 1, NOCCT
         DO JLOCC = 1, NLOCC
            IF (IOCCT .EQ. JLOCC) CLOCO(IOCCT,JLOCC) = D1
         END DO
      END DO
      IF (IPRSOS .GT. 5) THEN
         WRITE(LUPRI,'(/A)')' CLOCO matrix before localization: '
         CALL OUTPUT(CLOCO,1,NOCCT,1,NLOCC,NOCCT,NLOCC,1,LUPRI)
      END IF               
C
      DO IVIR = 1, NVIRT
         DO JLVIR = 1, NLVIR
            IF (IVIR .EQ. JLVIR) CLOCV(IVIR,JLVIR) = D1
         END DO
      END DO
C
      IF (IPRSOS .GT. 5) THEN
         WRITE(LUPRI,'(/A)')' CLOCV matrix before localization: '
         CALL OUTPUT(CLOCV,1,NVIRT,1,NLVIR,NVIRT,NLVIR,1,LUPRI)
      END IF              
C
C-----------------------------
C     Allocation of work space
C-----------------------------
C 
      KCMO   = 1
      KEND1  = KCMO   + NCMOT
      LWORK1 = LWORK  - KEND1
C     
      IF (LWORK1 .LT. 0) CALL STOPIT('GETLOC.1',' ',KEND1,LWORK)
C
C---------------------------------------------------
C     Open the file containing overlap matrix and...
C---------------------------------------------------
C
      call flshfo(lupri)
C     
C-------------------------
C     Read MO coefficients
C-------------------------
C
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ',' ',
     &                               IDUMMY,.FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
C
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READT (LUSIFC,NCMOT,WORK(KCMO))
C     
C***********************************************************
C     Generate transformation matrices to localized orbitals
C***********************************************************
C
      IF (ABALOC) THEN
C
C---------------------------------------------------------
C       MULLIKEN LOCALIZATION METHOD (According to ec. 31)
C---------------------------------------------------------
C
#ifdef NOT_FINISHED_YET
CSPAS:31.05.11: MULLOC is not yet working correctly
         IF(MULLOC) THEN
C
            MATCMO = 0
            DO  ISYM = 1,NSYM
               MATCMO = MATCMO + NBAS(ISYM)*NOCC(ISYM)
            END DO   
C
C---------------------
C           Allocation
C---------------------
C
            KCMOC  = KEND1
            KIAMN  = KCMOC+MATCMO
            KIPR   = KIAMN+NBAST
            KOVLP  = KIPR+NBAST
            KJTRAN = KOVLP+NNBAST
            KITRAN = KJTRAN+NBAST
            KCTRAN = KITRAN+NBAST*MXCGTO/IRAT
            KMOG   = KCTRAN+NBAST*MXCGTO
C
C----------------------------------
C           kss1 is for the overlap
C----------------------------------
C
            KSS1   = KMOG+NBAST*NOCCT
            KNUMB  = KSS1+NBAST*(NBAST+1)/2
            KLCMO  = KNUMB+NBAST
            KPAST  = KLCMO+NBAST*NOCCT
            KAIJ   = KPAST+NOCCT*NOCCT*NUCDEP
            KBIJ   = KAIJ+NOCCT*NOCCT
C
C-------------------------------------------------------------
C           kss2 is for the localization transformation matrix
C-------------------------------------------------------------
C
            KSS2   = KBIJ+NOCCT*NOCCT
            KEND6  = KSS2+NOCCT*NOCCT
            LWORK6 = LWORK  - KEND6
            IF (LWORK6 .LT. 0) CALL STOPIT('GETLOC.6',' ',KEND6,LWORK)
C
C---------------------------------------------------------
C           End of allocation and copy MO's to WORK(KCMOC)
C---------------------------------------------------------
C
            JCMO  =  KCMOC
C
            DO 47 ISYM=1,NSYM
               NOCCI=NOCC(ISYM)
               NBASI=NBAS(ISYM)
               IF(NOCCI.EQ.0) GOTO 47
            CALL DCOPY(NOCCI*NBASI,WORK(KCMO+ICMO(ISYM)),1,WORK(JCMO),1)
               JCMO = JCMO + NOCCI*NBASI
   47       CONTINUE
C
C----------------------------------------------------------
C           Necessary information for poptra transformation
C----------------------------------------------------------
C
            CALL READSYM(IPRSOS,WORK(KJTRAN),WORK(KITRAN),WORK(KCTRAN),
     &                   WORK(KIAMN),WORK(KIPR))
C
C--------------------------------------------------------
C           Here we make the transformation from symmetry 
C           orbitals to cgto's orbitals (POPTRA)
C--------------------------------------------------------
C
            CALL POPTRA (NNBAST,MXCGTO,MATCMO,WORK(KCMOC),WORK(KOVLP),
     &                   WORK(KJTRAN),WORK(KITRAN),WORK(KCTRAN),
     &                   WORK(KMOG),WORK(KSS1))
C
C--------------------------------------------------- 
C           Here we calculated the matrix population 
C           PAST betwen two molecular orbitals
C---------------------------------------------------
C
            CALL POPMAT (IPRSOS,WORK(KPAST),WORK(KIAMN),WORK(KIPR),
     &                   WORK(KNUMB),WORK(KMOG),WORK(KSS1))
C
C-----------------------------------------------
C           Finally we proceed with localization 
C           method and we obtain SS -> CLOCO
C-----------------------------------------------
C
            CALL MULLLOC(IPRSOS,CLOCO,THRLOC,WORK(KPAST),WORK(KSS2),
     &                   WORK(KBIJ),WORK(KAIJ),NLOCC)
C
C---------------------------------------
C           Print the localized orbitals
C---------------------------------------
C            
            CALL DZERO(WORK(KLCMO),NORBT*NOCCT)
C     
            CALL DGEMM('N','N',NORBT,NLOCC,NOCCT,D1,WORK(KCMOC),NORBT,
     &                 CLOCO,NOCCT,D0,WORK(KLCMO),NORBT)
C
            CALL HEADER('Localized Mulliken population Orbitals',-1)
            CALL PLOCORB(WORK(KLCMO),.TRUE.,NLOCC,LUW4)
C
         END IF
CKeinSPASmehr
#endif
C
CPFP Nov/21/04
C        IF ( FOSOCC .OR. FBOVIR .OR. FBSETV ) THEN
         IF ( FOSOCC .OR. FBOOCC .OR. FBOCIN 
     &               .OR. FBOVIR .OR. FBSETV ) THEN
Cend-PFP
C
C============================================
C           Calculate dipole length integrals
C============================================
C
C------------------------------------
C           Allocation of work space.
C------------------------------------
C     
            KDLAB  = KEND1
            KIDSYM = KDLAB  + 3
            KIDADR = KIDSYM + (3+1)/IRAT
            KPRPAO = KIDADR + (9*MXCENT+1)/IRAT
            KTMP   = KPRPAO + N2BASX
            KEND2  = KTMP   + NNBASX
            LWORK2 = LWORK  - KEND2
C     
            IF (LWORK2 .LT. 0) CALL STOPIT('GETLOC.2',' ',KEND2,LWORK)
C     
            NLBTOT = 0
            NCOMP  = 0
            NPATOM = 0
C
C-----------------------------------------------
C           Calculate AO dipole length integrals
C-----------------------------------------------
C
            IPRINT = IPRSOS - 20
            CALL GET1IN(DUMMY,'DIPLEN ',NCOMP,WORK(KEND2),LWORK2,
     &                  WORK(KDLAB),WORK(KIDSYM),WORK(KIDADR),
     &                  IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &                  DUMMY,IPRINT)
C
            NLAB = 3
            CALL LABCOP(NLAB,NLBTOT,WORK(KDLAB),WORK(KIDSYM),LABAPP,
     &                  LABSYM)
C
cspaS       LUPROP = 0
C
            CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
C
C---------------------------------------------
C           Read in AO dipole length integrals
C---------------------------------------------
C
            DO IPRLBL = 1, NLBTOT
               LABEL = LABAPP(IPRLBL)
               KSYM = LABSYM(IPRLBL)
C     
               IF (LABEL(1:1) .EQ. 'X') IDIP = 1
               IF (LABEL(1:1) .EQ. 'Y') IDIP = 2
               IF (LABEL(1:1) .EQ. 'Z') IDIP = 3
C
               REWIND LUPROP
C
               IF (FNDLB2(LABEL,RTNLBL,LUPROP,LUERR)) THEN
                  IF (RTNLBL(2) .EQ. 'SYMMETRI') THEN
                     ANTSYM = D1

                     CALL READT(LUPROP,NNBASX,WORK(KTMP))
C
                     CALL DSPTSI(NBAST,WORK(KTMP),WORK(KPRPAO))
                  ELSE
                     CALL QUIT('Error: No antisymmetry label on LUPROP')
                  END IF
               ELSE
                  WRITE (LUPRI,'(//3A)') ' PROPERTY "',LABEL,
     &                                   '" NOT FOUND ON LUPROP.'
                  CALL QUIT('PROPERTY NOT FOUND ON LUPROP.')
               END IF
C     
               IF (IPRSOS .GT. 10) THEN 
                  WRITE(LUPRI,'(/3A)') ' PROPERTY MATRIX ',LABEL,
     &                              ' IN AO BASIS AS READ FROM LUPROP '
                  CALL OUTPAK(WORK(KEND2),NBAST,1,LUPRI)
                  WRITE(LUPRI,'(/3A)') ' PROPERTY MATRIX ',LABEL,
     &                       ' IN AO BASIS AFTER (ANTI)SYMMETRIZATION '
                  CALL OUTPUT(WORK(KPRPAO),1,NBAST,1,NBAST,NBAST,NBAST,
     &                        1,LUPRI)
               END IF
C
C--------------------------------------------------------
C           Transform dipole length integrals to MO basis
C--------------------------------------------------------
C
               CALL DZERO(PRPMO(1,1,IDIP),NORBT*NORBT)
C
               DO ISYM = 1, NSYM
C
                  JSYM = MULD2H(ISYM,KSYM)
C
                  IF ((ISYM .GE. JSYM) .AND. 
     &              (NORB(ISYM) .GT. 0) .AND. (NORB(JSYM) .GT. 0)) THEN
C     
                     IF (LWORK2 .LT. NBAS(ISYM)) 
     &                   CALL STOPIT('GETLOC.2.1',' ',KEND2+LWORK2,
     &                               KEND2+NBAS(ISYM))
C     
                     CALL UTHV(WORK(KCMO+ICMO(ISYM)),WORK(KPRPAO),
     &                      WORK(KCMO+ICMO(JSYM)),ISYM,JSYM,NBAS(ISYM),
     &                         NBAS(JSYM),PRPMO(1,1,IDIP),WORK(KEND2))
C
                     IF (IPRSOS. GT. 5) THEN 
                        WRITE(LUPR I,'(/,A,I5,A,I5)') ' ISYM= ',ISYM,
     &                                                ' JSYM= ',JSYM
                        WRITE(LUPRI,'(/4A)') 
     &                      ' PROPERTY: ',LABEL,' IN MO. BASIS',
     &                      ' BEFORE (ANTI)SYMMETRIZATION '
                        CALL OUTPUT(PRPMO(1,1,IDIP),1,NORBT,1,NORBT,
     &                              NORBT,NORBT,1,LUPRI)
                     END IF
C
                     IF (IPRSOS .GT. 10) THEN 
                        WRITE(LUPRI,'(/A,I5,A)')
     &                      ' MO. COEFFICIENTS FOR SYMMETRY',ISYM
                        CALL OUTPUT(WORK(KCMO+ICMO(ISYM)),1,NBAS(ISYM),
     &                            1,NORB(ISYM),NBAS(ISYM),NORB(ISYM),1,
     &                              LUPRI)
                        IF (ISYM .NE. JSYM) THEN
                           WRITE(LUPRI,'(/A,I5,A)')
     &                         ' MO. COEFFICIENTS FOR SYMMETRY',JSYM
                           CALL OUTPUT(WORK(KCMO+ICMO(JSYM)),1,
     &                                 NBAS(JSYM),1,NORB(JSYM),
     &                                 NBAS(JSYM),NORB(JSYM),1,LUPRI)
                        END IF
                     END IF
C     
                  END IF
C     
               END DO
C     
               IF (KSYM .GT. 1) THEN
CSPAS:27/3-06:TRANSA is changed to TRANSX
C                 CALL TRANSA(PRPMO(1,1,IDIP),PRPMO(1,1,IDIP),NORBT,
C    &                        NORBT,ANTSYM)
                  CALL TRANSX(PRPMO(1,1,IDIP),PRPMO(1,1,IDIP),NORBT,
     &                        NORBT,ANTSYM,IPRSOS)
CKeinSPASmehr
               END IF
C
               IF (IPRSOS .GE. 5) THEN
               WRITE(LUPRI,'(/4A)')' PROPERTY: ',LABEL,' IN MO. BASIS'
                  CALL OUTPUT(PRPMO(1,1,IDIP),1,NORBT,1,NORBT,
     &                        NORBT,NORBT,1,LUPRI)
               END IF
      
            END DO
C
         CALL GPCLOSE(LUPROP,'KEEP')
C
        END IF
C    
C============================================================
C     Generate transformation matrix to localized orbitals
C============================================================
C     
CPFP Nov/21/04
C        IF ( FOSOCC ) THEN
         IF ( FOSOCC .OR. FBOCIN .OR. FBOOCC ) THEN
Cend-PFP
C     
C-----------------------------------
C     Allocation of work space   
C-----------------------------------
C
            KRM    = KEND1
            KSS    = KRM   + NOCCT * NOCCT * 3
            KCC    = KSS   + NOCCT * NOCCT
            KDD    = KCC   + NOCCT * NOCCT
            KLCMO  = KDD   + NOCCT * NOCCT
            KEND3  = KLCMO + NORBT * NOCCT
            LWORK3 = LWORK - KEND3
C     
            IF (LWORK3 .LT. 0) CALL STOPIT('GETLOC.3',' ',KEND3,LWORK)
C     
C------------------------------------------------------------
C     Foster Boys localization of the occupied orbitals
C------------------------------------------------------------
C     
            CALL FOSBOY(IPRSOS,PRPMO,NLOCC,CLOCO,THRLOC,
     &          WORK(KRM),WORK(KSS),WORK(KCC),WORK(KDD))

C
C---------------------------------------------
C     Print localized molecular orbitals    
C---------------------------------------------
C 
         call flshfo(lupri)   
            CALL DZERO(WORK(KLCMO),NORBT*NOCCT)
C
            CALL DGEMM('N','N',NORBT,NLOCC,NOCCT,D1,WORK(KCMO),NORBT,
     &          CLOCO,NOCCT,D0,WORK(KLCMO),NORBT)
C     
            CALL HEADER('Localized Occupied Molecular Orbitals',-1)
            CALL PLOCORB(WORK(KLCMO),.TRUE.,NLOCC,LUPRI)
C
         call flshfo(lupri)
         END IF
C    
C==================================================================
C     Generate transformation matrix to localized virtuals orbitals
C===================================================================
C   
         IF (FBOVIR .OR. FBSETV .OR. FBSTVO) THEN
C     
C-----------------------------
C     Allocation of work space 
C-----------------------------
C
            KRM    = KEND1
            KTT    = KRM   + NVIRT*NOCCT*3
            KLCMV  = KTT   + NVIRT*NVIRT
            KEND4  = KLCMV + NORBT*NLVIR
            LWORK4 = LWORK - KEND4
C     
            IF (LWORK4 .LT. 0) CALL STOPIT('GETLOC.4',' ',KEND4,LWORK)
C
C------------------------------------------------------
C     Intensity localization of the virtual orbitals
C------------------------------------------------------
C     
            call flshfo(lupri)
C
            CALL LOCALV(IPRSOS,PRPMO,NLOCC,CLOCO,NLVIR,CLOCV,
     &                  WORK(KRM),WORK(KTT),WORK(KEND4),LWORK4)
C     
C---------------------------------------------
C     Print localized molecular orbitals    
C---------------------------------------------
C     
            CALL DZERO(WORK(KLCMV),NORBT*NLVIR)
C     
            IF (IPRSOS .GT. 5) THEN
               CALL OUTPUT(WORK(KCMO),1,NORBT,1,NORBT,
     &                     NORBT,NORBT,1,LUPRI)
C     
               CALL OUTPUT(WORK(KCMO+NORBT*NOCCT),1,NORBT,1,NVIRT,
     &                     NORBT,NVIRT,1,LUPRI)
            END IF
C
            CALL DGEMM('N','N',NORBT,NLVIR,NVIRT,D1,
     &                 WORK(KCMO+NORBT*NOCCT),NORBT,
     &                 CLOCV,NVIRT,D0,WORK(KLCMV),NORBT)
C
            CALL HEADER('Localized Virtual Molecular Orbitals',-1)
            CALL PLOCORB(WORK(KLCMV),.FALSE.,NLVIR,LUPRI)
C     
         END IF
C
#ifdef NOT_FINISHED_YET
CSPAS:2.6.2011: this breaks the parallel code
         CALL NEWMOLDEN(IPRSOS,WORK(KLCMO),OCCUP)
#endif
C
      END IF
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      CALL QEXIT('GETLOC')
C
      RETURN
      END

C**********************************************************************
C     /* Deck fosboy */
      SUBROUTINE FOSBOY(IPRSOS,PRPMO,NLOCC,CLOCO,THRLOC,RM,SS,CC,DD)
C     
C      Stephan P. A. Sauer 20/2-1999
C      essentially stolen from the RPAC program
C
C      This routine does a FOSTER-BOYS LOCALIZATION
C      implemented after the equations from 
C      T.D.Bouman,B.Voigt,Aa.E.Hansen, JACS 101,550(1979)
C      
C
C     PRPMO  dipole length integrals in MO basis
C     RM     work array for dipole length integrals over occupied MOs
C     CLOCO  rectangular transformation matrix to be constructed
C            with the first index over all molecular orbitals
C     SS     square transformation matrix to be constructed
C     THRLOC Threshold for localization
C     CC     WORK ARRAY, LATER CONTAINING CENTROIDS
C     DD     WORK ARRAY
C
C     DR     INITIAL SUM OF (R**2)
C     DS     FINAL SUM
C     NR     ITERATION COUNTER
C
#include "implicit.h"
C
C LUPRI used from include file priunit.h
C MAXORB is needed in COMMON /INFIND/
C MAXASH is needed in COMMON /INFIND/
#include "maxorb.h"
#include "maxash.h"
C
C NSYM, NOCCT, NORBT, NRHF(8), NVIR(8) used from COMMON /INFORB/
C ISX(MAXORB) used from COMMON /INFIND/
#include "priunit.h"
#include "inforb.h"
#include "infind.h"
CPFP
#include "locinf.h"
Cend-PFP
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
      PARAMETER (DP25 = 0.25D0, THRSLD = 1.0D-15)
C
      DIMENSION   PRPMO(NORBT,NORBT,3)
      DIMENSION   RM(NOCCT,NOCCT,3)
      DIMENSION   CLOCO(NOCCT,NLOCC)
      DIMENSION   SS(NOCCT,NOCCT)
      DIMENSION   CC(NOCCT,NOCCT)
      DIMENSION   DD(NOCCT,NOCCT)
C
      LOGICAL     CONV
C
      CALL QENTER('FOSBOY')
C
C----------------------------------------
C     initialize some numerical constants
C----------------------------------------
C
      PID4 = ATAN(D1)
      PID8 = PID4/D2
      THR2 = THRLOC**2 
C
CPFP Nov/21/04
      IF ( FOSOCC .OR. FBOCIN .OR. FBOOCC) THEN
C
C------------------------------------------------------
C     copy the dipole integrals over occupied MOs to RM
C------------------------------------------------------
C
      CALL DZERO(RM,NOCCT*NOCCT*3)
C
      DO IDIP = 1, 3
C
         DO IOCCT = 1, NOCCT
            IORBT = ISX(IOCCT)
C     
            DO JOCCT = 1, NOCCT
               JORBT = ISX(JOCCT)
               RM(IOCCT,JOCCT,IDIP) = PRPMO(IORBT,JORBT,IDIP)
            END DO
         END DO
C
         IF (IPRSOS .GE. 5) THEN
            IF (IDIP .EQ. 1) 
     &         WRITE(LUPRI,'(/A)')' occupied RM matrix : X'
            IF (IDIP .EQ. 2) 
     &         WRITE(LUPRI,'(/A)')' occupied RM matrix : Y'
            IF (IDIP .EQ. 3) 
     &         WRITE(LUPRI,'(/A)')' occupied RM matrix : Z'
            CALL OUTPUT(RM(1,1,IDIP),1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,
     &                  LUPRI)
         END IF     
C
CPFP  For FBOOCC only
C
         IF (FBOOCC) THEN
C
           DO IO2L = 1, NO2LOC
             DO IOC = 1, NOCCT
               DO JOC = 1, NOCCT
                 IF ((IOC.EQ.NTOC2L(IO2L)).OR.(JOC.EQ.NTOC2L(IO2L))) 
     &             RM(IOC,JOC,IDIP) = D0
               END DO
             END DO
           END DO
C
          IF (IPRSOS .GE. 5) THEN
            IF (IDIP .EQ. 1) 
     &         WRITE(LUPRI,'(/A)')' occupied RM.2 matrix : X'
            IF (IDIP .EQ. 2) 
     &         WRITE(LUPRI,'(/A)')' occupied RM.2 matrix : Y'
            IF (IDIP .EQ. 3) 
     &         WRITE(LUPRI,'(/A)')' occupied RM.2 matrix : Z'
            CALL OUTPUT(RM(1,1,IDIP),1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,
     &                  LUPRI)
          END IF
C
        END IF
C
Cend-PFP
C
      END DO
C
C=============================================
C     initialize variables for first iteration
C=============================================
C
      NR = 0
C
      DS   = D0
      DMAX = D0
C
      CALL DZERO(SS,NOCCT*NOCCT)
      CALL DZERO(DD,NOCCT*NOCCT)
      CALL DZERO(CC,NOCCT*NOCCT)
C
      DO J = 1, NOCCT
         DD(J,J) = RM(J,J,1)**2 + RM(J,J,2)**2 + RM(J,J,3)**2
         DS = DS + DD(J,J)
         IF (J .LT. NOCCT)  THEN
            K = J + 1
            DO I = K, NOCCT
               CALL DAB (I,J,RM,CC,DD)
               DMAX = DMAX + DD(I,J)
            END DO
         END IF
      END DO
C
      IF (IPRSOS .GE. 10) THEN
         WRITE(LUPRI,'(/A)')' DD matrix :'
         CALL OUTPUT(DD,1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,LUPRI)
      END IF     
C
C--------------------------------------
C     DR = INITIAL VALUE OF SUM OF R**2
C--------------------------------------
C
      DR = DS
C
C---------------------------------
C     INITIALIZE SS TO UNIT MATRIX
C---------------------------------
C
      CALL DUNIT (SS,NOCCT)
C
C================================================
C     SUCCESSIVE 2*2 TRANSFORMATIONS USING EQ. 41
C================================================
C
      NLOC = 1
      FLNO = NOCCT - NLOC + 1
C
  300 DMAX = DMAX / FLNO
      IF (DMAX .LT. THRLOC)  DMAX = THRLOC
C

  400 CONV = .TRUE.
C
           DO I = NLOC + 1, NOCCT
C
              DO J = NLOC, I - 1
C
C--------------------
C     J .LT. I Always
C--------------------
C
                 IF (DD(I,J) .GT. DMAX) THEN
C
                    CONV = .FALSE.
                    NR = NR + 1
                    C1 = CC(J,I)
                    C2 = CC(I,J)
C
                    IF (ABS(C1) .GE. THR2) THEN
                       THETA = -DP25 * ATAN(C2/C1)
C
C---------------------------------------
C     Find value of theta giving maximun
C---------------------------------------
C
                       IF (C1 .LT. D0) THEN
                          THETA = THETA - SIGN(PID4,C2)
                          IF (C2 .EQ. D0) THETA = PID4
                       END IF
                    ELSE
C
C--------------------------
C     Special case - C1 = 0
C--------------------------
C
                       THETA = -SIGN(PID8,C2)
                    END IF
C
                    DT = ABS(THETA)
C
C----------------------
C     Conserve Symmetry
C----------------------
C
                    IF (ABS(DT-PID4).LE.THRLOC) THETA = SIGN(PID4,THETA)
                    IF (ABS(DT-PID8).LE.THRLOC) THETA = SIGN(PID8,THETA)
C
C-------------------
C     2 X 2 Rotation
C-------------------
C
                    CALL LOCROT(I,J,THETA,NLOC,RM,SS)
C
C---------------------------
C     UPDATE DD(I,J) AND SUM
C---------------------------
C
                    RIJ = DD(I,J)
                    CALL DAB(I,J,RM,CC,DD)
                    DS = DS +  RIJ - DD(I,J)
C
C-------------------------------------------
C     UPDATE ROW AND COLUMN I,J IN CC AND DD
C-------------------------------------------
C
                    DO K = NLOC, NOCCT
                       IF (J .NE. K .AND. I .NE. K) THEN
                          CALL DAB (MAX0(J,K),MIN0(J,K),RM,CC,DD)
                          CALL DAB (MAX0(I,K),MIN0(I,K),RM,CC,DD)
                       END IF
                    END DO
C
                 END IF
C
              END DO
           END DO
C     
      IF (.NOT. CONV) GO TO 400
C
      IF (DMAX .GT. THRLOC)  GO TO 300
C
C========================
C     CONVERGENCE REACHED
C========================
C
      END IF
C
C-----------------------------------------
C     Copy the transformation matrix in SS 
C     to the right positions in CLOCO
C-----------------------------------------
C
      DO IOC = 1, NOCCT
        DO JOC = 1, NLOCC
          CLOCO(IOC,JOC) = SS(IOC,JOC)
        END DO
      END DO
C
CPFP For FBOCIN only 
C
      IF ( FBOCIN ) THEN
        DO IO2L = 1, NO2LOC
C
          DO IOC = 1, NOCCT
            DO JOC = 1, NLOCC
              IF ((IOC.EQ.NTOC2L(IO2L)).OR.(JOC.EQ.NTOC2L(IO2L))) THEN
                 CLOCO(IOC,JOC) = D0
                 IF (IOC.EQ.JOC) CLOCO(IOC,JOC) = D1
              END IF
            END DO
          END DO
C
        END DO
      END IF
C
Cend-PFP
C 
      IF (IPRSOS .GT. 10) THEN
         WRITE(LUPRI,'(/A)')' SS matrix : '
         CALL OUTPUT(SS,1,NOCCT,1,NOCCT,NOCCT,NOCCT,1,LUPRI)
      END IF
C
      IF (IPRSOS .GT. 5) THEN
         WRITE(LUPRI,'(/A)')' CLOCO matrix : '
         CALL OUTPUT(CLOCO,1,NOCCT,1,NLOCC,NOCCT,NLOCC,1,LUPRI)
      END IF
C
C-------------------------
C     COPY CENTROIDS TO CC
C-------------------------
C
C      DO J = 1, 3
C         DO I = 1, NOCCT
C            CC(I,J) = RM(I,I,J)
C         END DO
C      END DO
C
      DSS = DS - DR
      WRITE (LUPRI,'(//,A,/,A,/,A,G14.6,A,G14.6,/,A,I5,A,G9.2,//)') 
     &      ' FOSTER-BOYS LOCALIZATION:',
     &      ' SUM OF (R-VECTOR)**2 FOR THE ORBITAL DISTRIBUTIONS',
     &      ' INCREASED BY',DSS,'  TO',DS,
     &      ' NUMBER OF (2*2)-OPTIMIZATIONS:',NR,'  ACCURACY:',THRLOC
C
      CALL QEXIT('FOSBOY')
C
      RETURN
      END

C**********************************************************************
C  /* Deck dab */
      SUBROUTINE DAB(I,J,RM,CC,DD)
C
C      Stephan P. A. Sauer 20/2-1999
C      essentially stolen from the RPAC program
C
C      This routine calculates the 
C      NUMERATOR & DENOMINATOR OF EQ. 41 OF 
C      T.D.Bouman,B.Voigt,Aa.E.Hansen, JACS 101,550(1979)
C      
C
C     RM     DIPOLE LENGTH INTEGRALS
C     CC     CENTROIDS
C     DD     WORK ARRAY
C
#include "implicit.h"
#include "maxorb.h"
#include "maxash.h"
C
      PARAMETER (D4 = 4.0D+00)
C
      DIMENSION   RM(NOCCT,NOCCT,3)
      DIMENSION   CC(NOCCT,NOCCT)
      DIMENSION   DD(NOCCT,NOCCT)
C
C NOCCT used from COMMON /INFORB/
#include "inforb.h"
C
      CALL QENTER('DAB')
C
      AX = RM(I,I,1) - RM(J,J,1)
      AY = RM(I,I,2) - RM(J,J,2)
      AZ = RM(I,I,3) - RM(J,J,3)
      BX = RM(J,I,1)
      BY = RM(J,I,2)
      BZ = RM(J,I,3)
      C1 = (AX*AX + AY*AY + AZ*AZ)/D4 - BX*BX - BY*BY - BZ*BZ
      C2 = AX*BX + AY*BY + AZ*BZ
      CC(J,I) = C1
      CC(I,J) = C2
C
C--------------------------------------------------------
C     CHANGE IN (R(I) - R(J))**2.  R(I), R(J) = CENTROIDS
C--------------------------------------------------------
C
      DD(I,J) = SQRT(C1*C1 + C2*C2) - C1
C
      CALL QEXIT('DAB')
C
      RETURN
      END

C***********************************************************************
C  /* Deck locrot */
      SUBROUTINE LOCROT(I,J,THETA,NLOC,RM,SS)
C
C      Stephan P. A. Sauer 20/2-1999
C      essentially stolen from the RPAC program
C
C      This routine calculates the 
C      PERFORMS UNITARY TRANSFORMATIONS OF ORBITALS I, J
C
C      THETA  TRANSFORMATION ANGLE
C      RM     DIPOLE LENGTH INTEGRALS TO BE UPDATED
C      SS     TRANSFORMATION MATRIX TO BE UPDATED
C
#include "implicit.h"
#include "maxash.h"
#include "maxorb.h"
C
      DIMENSION   RM(NOCCT,NOCCT,3)
      DIMENSION   SS(NOCCT,NOCCT)
C
C JTINAC, JTSEC, IOBTYP, ISW(MAXORB) used from COMMON /INFIND/
C JWOP used from COMMON /INFVAR/
C NOCCT used from COMMON /INFORB/
#include "infind.h"
#include "infvar.h"
#include "inforb.h"
C
      CALL QENTER('LOCROT')
C
      CA = COS(THETA)
      SA = SIN(THETA)
C
      DO KK = NLOC, NOCCT
C
         TI = SS(KK,I)
         TJ = SS(KK,J)
         SS(KK,J) = SA*TI + CA*TJ
         SS(KK,I) = CA*TI - SA*TJ
C
         DO L = 1, 3
C
            TI = RM(KK,I,L)
            TJ = RM(KK,J,L)
            RM(KK,J,L) = SA*TI + CA*TJ
            RM(KK,I,L) = CA*TI - SA*TJ
C
         END DO
C
      END DO
C     
      DO KK = NLOC, NOCCT
C
         II = (KK-I) * (KK-J)
C
         DO L = 1, 3
C
            IF (II.NE.0) THEN
C
               RM(I,KK,L) = RM(KK,I,L)
               RM(J,KK,L) = RM(KK,J,L)
C
            ELSE
C
               TI = RM(I,KK,L)
               TJ = RM(J,KK,L)
               RM(J,KK,L) = SA*TI + CA*TJ
               RM(I,KK,L) = CA*TI - SA*TJ
C
            END IF
C
         END DO
C
      END DO
C
      CALL QEXIT('LOCROT')
C
      RETURN
      END

C***********************************************************************
C  /* Deck locai */
      SUBROUTINE LOCAI(NLOCC,NLVIR,NVARPT,CLOCO,CLOCV,GD,GDLOC)
C
C      Stephan P. A. Sauer 29/10-1997
C
C      This routine transforms a gradient <a|O|i> to localized orbitals.
C
C
#include "implicit.h"
C
C MAXASH used from include file maxash.h
C MAXORB, MAXOCC used from include file maxorb.h
#include "maxash.h"
#include "maxorb.h"
C
C JTINAC, JTSEC, IOBTYP, ISW(MAXORB) used from COMMON /INFIND/
C JWOP used from COMMON /INFVAR/
C NOCCT, NVIRT used from COMMON /INFORB/
#include "infind.h"
#include "infvar.h"
#include "inforb.h"
C
      DIMENSION CLOCO(NOCCT,NLOCC)
      DIMENSION CLOCV(NVIRT,NLOCC)
      DIMENSION GD(NVARPT)
      DIMENSION GDLOC(NLVIR,NLOCC)
C
C---------------------------------------------------------
C     Transform to localized orbitals
C---------------------------------------------------------
C
      CALL QENTER('LOCAI')
C
      CALL DZERO (GDLOC,NLOCC*NLVIR)
C
      DO IVAR = 1, NVARPT
         M = JWOP(1,IVAR)
         N = JWOP(2,IVAR)
         MTYP = IOBTYP(M)
         NTYP = IOBTYP(N)
         IF (MTYP .EQ. JTINAC) I  = ISW(M)
         IF (NTYP .EQ. JTSEC)  IA = ISW(N) - NOCCT
         IF (NTYP .EQ. JTINAC) I  = ISW(N)
         IF (MTYP .EQ. JTSEC)  IA = ISW(M) - NOCCT
C     
         DO ILOCC = 1, NLOCC
            DO ILVIR = 1, NLVIR
                GDLOC(ILVIR,ILOCC) = GDLOC(ILVIR,ILOCC)
     &                   + GD(IVAR) * CLOCO(I,ILOCC) * CLOCV(IA,ILVIR)
            END DO
         END DO
      END DO
C
      CALL QEXIT('LOCAI')
C
      RETURN
      END             
      
C***********************************************************************
C  /* Deck localv */
      SUBROUTINE LOCALV(IPRSOS,PRPMO,NLOCC,CLOCO,NLVIR,CLOCV,
     &                  RM,TT,WORK,LWORK)
C
C      Stephan P. A. Sauer 23/2-1999
C      essentially stolen from the RPAC program
C
C      This routine calculates the localized virtual orbitals. 
C      The following criterion is used:
C                the dipole strenghs is maximized between the
C                canonical virtual orbitals and the already 
C                localized valence occupied orbitals.
C
C      PRPMO  ARRAY CONTAINING DIPOLE LENGTH INTEGRALS
C      CLOCO  TRANSFORMATION MATRIX FOR OCCUPIED CANONICAL TO LOCALIZED
C             MO'S.
C      TT   SAME FOR VIRTUAL MO'S.
C      VEC  WORK ARRAY, NVIRT*NVIRT
C      TM   DITTO
C
#include "implicit.h"
C
C MAXASH used from include file maxash.h
C MAXORB, MAXOCC used from include file maxorb.h
#include "maxash.h"
#include "maxorb.h"
C
C JWOP used from COMMON /INFVAR/
C NOCCT, NVIRT used from COMMON /INFORB/
C NWOPPT used from COMMON /INFLIN/
#include "priunit.h"
#include "infvar.h"
#include "inforb.h"
#include "inflin.h"
#include "locinf.h"
C
      PARAMETER ( S0 = 0, D0 = 0.0D+00, D1 = 1.0D+00 )
C
      DIMENSION   PRPMO(NORBT,NORBT,3)
      DIMENSION   RM(NVIRT,NOCCT,3)
      DIMENSION   TT(NVIRT,NVIRT)
      DIMENSION   CLOCO(NOCCT,NLOCC)
      DIMENSION   CLOCV(NVIRT,NLVIR)
      DIMENSION   WORK(LWORK)
C
      CALL QENTER('LOCALV')
C
C=========================
C     initialize variables
C=========================
C
      CALL DUNIT(TT,NVIRT)
C
C-----------------------
C     Allocate workspace
C-----------------------
C
      NVIRT2 = NVIRT*NVIRT
      NOCCT2 = NOCCT*NOCCT
C
      KRMW   = 1
      KDA    = KRMW   + NWOPPT
      KENGVA = KDA    + NVIRT2
      KENGVE = KENGVA + NVIRT
      KTMP1  = KENGVE + NVIRT2
      KTMP2  = KTMP1  + NVIRT
      KOVRLP = KTMP2  + NVIRT
      KWORK1 = KOVRLP + NVIRT2
      LWORK1 = LWORK  - KWORK1
C
      IF (LWORK1 .LT. 0) CALL STOPIT('LOCALV.1',' ',KWORK1,LWORK)
C
C------------------------------------------------------
C     copy the dipole transition moment integrals to RM
C------------------------------------------------------
C
      DO IDIP = 1, 3
C
         CALL DZERO(WORK(KRMW),NWOPPT)
C         
         DO IVAR = 1, NWOPPT
C
            M = JWOP(1,IVAR)
            N = JWOP(2,IVAR)
C
            WORK(KRMW-1+IVAR) = PRPMO(M,N,IDIP)
C
         END DO
C
        IF (IPRSOS .GE. 50) THEN
            IF (IDIP .EQ. 1) 
     &         WRITE(LUPRI,'(/2A)')' LOCALV before localization :',
     &                             ' transition RM matrix : X'
            IF (IDIP .EQ. 2) 
     &         WRITE(LUPRI,'(/2A)')' LOCALV before localization :',
     &                             ' transition RM matrix : Y'
            IF (IDIP .EQ. 3) 
     &         WRITE(LUPRI,'(/2A)')' LOCALV before localization :',
     &                             ' transition RM matrix : Z'
            CALL OUTPUT(WORK(KRMW),1,NWOPPT,1,1,NWOPPT,1,1,LUPRI)
         END IF     
C
C------------------------------------------------------
C     Form the matrices  IVAR localized, IDIP canonical
C------------------------------------------------------
C
         CALL LOCAI(NLOCC,NVIRT,NWOPPT,CLOCO,TT,WORK(KRMW),
     &              RM(1,1,IDIP))
C
         IF (IPRSOS .GE. 50) THEN
            IF (IDIP .EQ. 1) 
     &         WRITE(LUPRI,'(/2A)')' LOCALV after localization :',
     &                             ' transition RM matrix : X'
            IF (IDIP .EQ. 2) 
     &         WRITE(LUPRI,'(/2A)')' LOCALV after localization :',
     &                             ' transition RM matrix : Y'
            IF (IDIP .EQ. 3) 
     &         WRITE(LUPRI,'(/2A)')' LOCALV after localization :',
     &                             ' transition RM matrix : Z'
            CALL OUTPUT(RM(1,1,IDIP),1,NVIRT,1,NLOCC,NVIRT,NLOCC,1,
     &                  LUPRI)
         END IF
C
      END DO
C
C----------------------------------------------------------          
C     Form matrix  DA(M,M') = (M/R/LAM)*(LAM/R/M'), and 
C     diagonalize for each occupied orbital. Choose largest
C     eigenvalue and corresponding eigenvector.
C----------------------------------------------------------
C
      LSTOCC = MOD(NLVIR,NLOCC)
      IF (LSTOCC .EQ. 0) LSTOCC = NLOCC
      NCNT = S0
C
      DO  N = 1, NLOCC
         DO  J = 1, NVIRT
            IOFF = KDA - 1 + (J-1)*NVIRT
C
            DO  I = 1, NVIRT  
               WORK(IOFF+I) = RM(I,N,1)*RM(J,N,1)
     &                      + RM(I,N,2)*RM(J,N,2)
     &                      + RM(I,N,3)*RM(J,N,3)
            END DO
         END DO
C
         CALL RS(NVIRT,NVIRT,WORK(KDA),WORK(KENGVA),1,WORK(KENGVE),
     &           WORK(KTMP1),WORK(KTMP2),IERR)
C
         IF (FBSTVO) THEN
C
            NVSTX = S0
            DO M = 1, NOV
               NM = NTVI2L(M)
               IF (N .EQ. NM) THEN 
                  NVSTX = NVSET
                  NCNT = NCNT + 1
                  GO TO 101
               END IF
            END DO
C
 101        CONTINUE
C
            DO  L = 1, NVSTX
               KOFF = KENGVE - 1 + (NVIRT-L)*NVIRT
               DO  K = 1, NVIRT
                  CLOCV(K,(L-1)*NV2LOC+NCNT) = WORK(KOFF+K)
               END DO
            END DO
C
         ELSE
C
            NVSTX = NVSET
            IF (N .GT. LSTOCC) NVSTX =  NVSET - 1
C
            DO  L = 1, NVSTX
               KOFF = KENGVE - 1 + (NVIRT-L)*NVIRT
               DO  K = 1, NVIRT
                  CLOCV(K,(L-1)*NLOCC+N) = WORK(KOFF+K)
               END DO
            END DO
C
         END IF
C
      END DO
C
      IF (IPRSOS .GE. 10) THEN
         WRITE(LUPRI,'(/2A)')' CLOCV before L.Ortho.and overlap'
         CALL OUTPUT(CLOCV,1,NVIRT,1,NLVIR,NVIRT,NLVIR,1,LUPRI)
      ENDIF
C
C----------------------------------------------
C     Orthonormalization of each set separately
C----------------------------------------------
C
      IF (FBSTVO) THEN
C
         DO  I = 1, NVSET
            NOCC1 = NVIRT - (I-1)*NV2LOC
C
            IF (NOCC1 .GE. NLOCC) THEN
               NLOC1 = NV2LOC
            ELSE
               NLOC1 = NOCC1
            END IF
C
           CALL LOC_LOWDIN(IPRSOS,CLOCV(1,(I-1)*NV2LOC + 1),NVIRT,NLOC1,
     &                  LFLAG,WORK(KWORK1),LWORK1)
         END DO
C
      ELSE
C
         DO  I = 1, NVSET
            NOCC1 = NVIRT - (I-1)*NLOCC
C
            IF (NOCC1 .GE. NLOCC) THEN
               NLOC1 = NLOCC
            ELSE
               NLOC1 = NOCC1
            END IF
            CALL LOC_LOWDIN(IPRSOS,CLOCV(1,(I-1)*NLOCC + 1),NVIRT,NLOC1,
     &                  LFLAG,WORK(KWORK1),LWORK1)
         END DO
C
      END IF
C
C---------------------------------------------------------------
C     Orthonormalization of the later sets respect to the erlier
C---------------------------------------------------------------
C
      IF (FBSTVO) THEN
C
         IF (NVSET .GT. 1) THEN
            NVS1 = NVSET - 1
            DO  I = 1, NVS1
               NREST = ABS(NVIRT - (I+1)*NV2LOC)
               IF (NREST .GE. NV2LOC ) THEN
C
                  DO  J = 1, NV2LOC
                     NVEC = I*NV2LOC + J - 1
                     CALL ORTHON(IPRSOS,NVIRT,NVEC,CLOCV(1,1),
     &                           CLOCV(1,NVEC+1),IFLAG)
                   END DO
C
               ELSE
C
               NLOC1 = NV2LOC - NREST
                  DO  J = 1, NLOC1
                     NVEC = I*NV2LOC + J - 1
C
                     CALL ORTHON(IPRSOS,NVIRT,NVEC,CLOCV(1,1),
     &                           CLOCV(1,NVEC+1),IFLAG)
C
                  END DO
C
               END IF
C
            END DO
         END IF
C
      ELSE
C
         IF (NVSET .GT. 1) THEN
            NVS1 = NVSET - 1
            DO  I = 1, NVS1
               NREST = ABS(NVIRT - (I+1)*NLOCC)
               IF (NREST .GE. NLOCC ) THEN
C
                  DO  J = 1, NLOCC
                     NVEC = I*NLOCC + J - 1
                     CALL ORTHON(IPRSOS,NVIRT,NVEC,CLOCV(1,1),
     &                           CLOCV(1,NVEC+1),IFLAG)
                   END DO
C
               ELSE
C
               NLOC1 = NLOCC - NREST
                  DO  J = 1, NLOC1
                     NVEC = I*NLOCC + J - 1
C
                     CALL ORTHON(IPRSOS,NVIRT,NVEC,CLOCV(1,1),
     &                           CLOCV(1,NVEC+1),IFLAG)
C
                  END DO
C
               END IF
C
            END DO
         END IF
C
      END IF
C
C---------------------------------------
C     Calculation of the overlap matrix
C---------------------------------------
C
      DO L = 1, NVIRT
         IOFF = KOVRLP - 1 + (L-1)*NVIRT
         DO I = 1, NVIRT
            WORK(IOFF+I) = D0
            DO J = 1, NVIRT
               WORK(IOFF+I) = WORK(IOFF+I) + CLOCV(J,I)*CLOCV(J,L)
            END DO
         END DO
      END DO
C
      IF (IPRSOS .GT. 50) THEN
         WRITE(LUPRI,'(/A)')' OVERLAP matrix : in LOCALV '
         CALL OUTPUT(WORK(KOVRLP),1,NVIRT,1,NLVIR,NVIRT,NLVIR,1,LUPRI)
      END IF
C
C--------------------------------
C     Print matrixes TT and CLOCV
C--------------------------------
C
      IF (IPRSOS .GT. 10) THEN
         WRITE(LUPRI,'(/A)')' TT matrix : '
         CALL OUTPUT(TT,1,NVIRT,1,NLVIR,NVIRT,NLVIR,1,LUPRI)
      END IF
      IF (IPRSOS .GT. 5) THEN
         WRITE(LUPRI,'(/A)') ' CLOCV matrix : after orthonorm. '
         CALL OUTPUT(CLOCV,1,NVIRT,1,NLVIR,NVIRT,NLVIR,1,LUPRI)
      END IF
C
      CALL QEXIT('LOCALV')
C
      RETURN
      END

C*********************************************************************
C  /* Deck LOC_lowdin */  'STOLEN BY ... ' from RPAC
      SUBROUTINE LOC_LOWDIN(IPRSOS,OMAT,NOMAT,NSMAT,IFLAG,WORK,LWORK)
C
C      Stephan P. A. Sauer 18/12-1999
C
C     LOC_LOWDIN ORTHOGONALIZATION OF NVEC VECTORS STORED IN OMAT.
C
C     OMAT : Matrix containing Orbitals to be orthogonalized.
C     NOMAT: N. of rows.
C     NSMAT: N. of colomns.
C-------------------------------------------------------------------
C
#include "implicit.h"
C
C MAXASH used from include file maxash.h
C MAXORB, MAXOCC used from include file maxorb.h
#include "maxash.h"
#include "maxorb.h"
C
C JWOP used from COMMON /INFVAR/
C NOCCT, NVIRT used from COMMON /INFORB/
C NWOPPT used from COMMON /INFLIN/
#include "priunit.h"
#include "infvar.h"
#include "inforb.h"
#include "inflin.h"
C
      PARAMETER ( D0 = 0.0D+00, D1 = 1.0D+00 )
      PARAMETER ( DCTE = 1.0D-02 )
C
      DIMENSION   OMAT(NOMAT,NSMAT)
      DIMENSION   WORK(LWORK)
C
      CALL QENTER('LOC_LOWDIN')
C
C-----------------------
C     Allocate workspace
C-----------------------
C
      NSMAT2 = NSMAT*NSMAT
      NOMAT2 = NOMAT*NOMAT
C
      KEIGVE  = 1
      KEIGVA  = KEIGVE + NSMAT2
      KTMP1   = KEIGVA + NSMAT
      KTMP2   = KTMP1  + NSMAT
      KSS     = KTMP2  + NSMAT
      KTMP3   = KSS    + NSMAT2
      KOVRLP  = KTMP3  + NOMAT*NSMAT
      KWORK1  = KOVRLP + NSMAT2
      LWORK1  = LWORK  - KWORK1
C
      IF (LWORK1 .LT. 0) CALL STOPIT('LOC_LOWDIN.1',' ',LWORK1,LWORK)
C
      IFLAG = 0
C
C---------------------------------
C      Calculation of the overlaps
C---------------------------------
C
      DO L = 1, NSMAT
         IOFF = KOVRLP - 1 + (L-1)*NSMAT
         DO I = 1, NSMAT 
            WORK(IOFF+I) = D0
            DO J = 1, NOMAT
               WORK(IOFF+I) = WORK(IOFF+I) + OMAT(J,I)*OMAT(J,L)
            END DO
         END DO
      END DO
C
      IF (IPRSOS .GT. 10) THEN
         WRITE(LUPRI,'(/A)')' WORK(KOVRLP) matrix : overlaps matrix'
         CALL OUTPUT(WORK(KOVRLP),1,NSMAT,1,NSMAT,NSMAT,NSMAT,1,LUPRI)
      END IF
C
      CALL RS(NSMAT,NSMAT,WORK(KOVRLP),WORK(KEIGVA),1,WORK(KEIGVE),
     &         WORK(KTMP1),WORK(KTMP2),IERR)
C
           DO  I = 1, NSMAT
              IOFF = KEIGVA - 1 + I
              IF (WORK(IOFF).LE.DCTE) THEN
C
C     Stops Orthogonalization end gets out
C
                IFLAG = 1
C
                CALL QEXIT('LOC_LOWDIN')
C
                RETURN
              END IF
              WORK(IOFF) = D1/DSQRT(WORK(IOFF))
           END DO
C
      IF (IPRSOS .GT. 50) THEN
         WRITE(LUPRI,'(/A)')
     &         ' WORK(KEIGVA) matrix : after it does a=1/sqrt(a)'
         WRITE(LUPRI, '(F12.8)') (WORK(KEIGVA-1+i),i=1,NSMAT)
C
         WRITE(LUPRI,'(/A)') 
     &        ' WORK(KEIGVE) matrix : before multiplication '
         CALL OUTPUT(WORK(KEIGVE),1,NSMAT,1,NSMAT,NSMAT,NSMAT,1,LUPRI)
      END IF
C
           DO  I = 1, NSMAT
              ISS = KSS - 1 + (I-1)*NSMAT
C
              DO  J = 1, NSMAT
                 WORK(ISS+J) = D0
C
                 DO  K = 1, NSMAT
                    KOFF = KEIGVE - 1 + (K-1)*NSMAT
                    LOFF = KEIGVA - 1 + K
C
                    WORK(ISS+J) = WORK(ISS+J)
     &                          + WORK(KOFF+J)*WORK(LOFF)*WORK(KOFF+I)
C
                 END DO
              END DO
           END DO
C
      IF (IPRSOS .GT. 50) THEN
         WRITE(LUPRI,'(/A)') ' WORK(KSS) matrix : after multiplication'
         CALL OUTPUT(WORK(KSS),1,NSMAT,1,NSMAT,NSMAT,NSMAT,1,LUPRI)
         WRITE(LUPRI,'(/A)') ' OMAT matrix : before DGEMM'
         CALL OUTPUT(OMAT,1,NOMAT,1,NSMAT,NOMAT,NSMAT,1,LUPRI)
      END IF
C
C-----------------------------------
C     Generetes a new set of vectors
C-----------------------------------
C
      CALL DGEMM('N','N',NOMAT,NSMAT,NSMAT,D1,OMAT,NOMAT,WORK(KSS),
     &           NSMAT,D0,WORK(KTMP3),NOMAT)
C
      IF (IPRSOS .GT. 10) THEN
         WRITE(LUPRI,'(/A)') ' WORK(KTMP3) matrix: after DGEMM '
         CALL OUTPUT(WORK(KTMP3),1,NOMAT,1,NSMAT,NOMAT,NSMAT,1,LUPRI)
      END IF
C
       DO  I = 1, NSMAT
          IOFF = KTMP3 - 1 + (I-1)*NOMAT
          DO  J = 1, NOMAT
             OMAT(J,I) =  WORK(IOFF+J)
          END DO
       END DO
C
      CALL QEXIT('LOC_LOWDIN')
C
      RETURN
      END

C***********************************************************************
C  /* Deck plocorb */
      SUBROUTINE PLOCORB(CMO,PROCC,KSET,IOUT)
C
C      Stephan P. A. Sauer 18/12-1999
C
C
C Purpose:
C  Print localized molecular orbital coefficients on unit IOUT.
C
C Input:
C  CMO: MO orbital coefficients (symmetry blocked)
C  PROCC: If 'true' print only occupied orbitals
C         If 'false' print only virtual orbitals 
C  KSET: Number of rows of either OMO or VMO (NLOCC or NLVIR).
C  IOUT: output file unit
C
#include "implicit.h"
      DIMENSION CMO(*)
      LOGICAL   PROCC
      CHARACTER *8 ORBKND
C
C Used from common blocks:
C   INFINP : CENT,TYPE,SUPSYM,?
C   INFORB : NSYM,...
C   INFIND : ISSMO(),
C
#include "maxorb.h"
#include "maxash.h"
#include "infinp.h"
#include "inforb.h"
#include "infind.h"
C
      CALL QENTER('PLOCORB')
C
C---------------------------------------------
C     FORMAT statements for molecular orbitals
C---------------------------------------------
C
  100 FORMAT(/5X,I2,'  Localized ',A8,' Molecular Orbitals ')
  200 FORMAT(/' Orbital  ',5X,7I9)
  300 FORMAT(1X,I3,2X,A4,2X,A4,2X,7F9.4)
C
C-----------------------------
C     Print molecular orbitals
C-----------------------------
C
      ISTBAS = 0
      DO  ISYM = 1,NSYM
         IF (PROCC) THEN
            NENDI = KSET
            ORBKND = 'Occupied'
         ELSE
            NENDI = KSET
            ORBKND = 'Virtual '
         END IF
         NBASI = NBAS(ISYM)
      IF (NENDI.EQ.0) GO TO 400
         WRITE(IOUT,100) KSET,ORBKND
C
         ICMOI  = ICMO(ISYM)
         ISTORB = IORB(ISYM)
         IEND   = 0
  500    IST    = IEND + 1
         ISTMO  = IEND*NBASI + ICMOI
         IEND   = IEND + 7
         IF(IEND.GT.NENDI) IEND = NENDI
         IEMO   = NBASI*(IEND - 1) + ICMOI
         WRITE(IOUT,200) (I,I=IST,IEND)
C
            DO  I = 1, NBASI
               JSMO=ISTMO+I
               JEMO=IEMO+I
               WRITE(IOUT,300) I,CENT(I+ISTBAS),TYPE(I+ISTBAS),
     &              (CMO(J),J=JSMO,JEMO,NBASI)
            END DO
         IF (IEND.NE.NENDI) GO TO 500
C
 400  CONTINUE
        ISTBAS = ISTBAS + NBASI
      END DO
C
      CALL QEXIT('PLOCORB')
C
      RETURN
      END

C***********************************************************************
C  /* Deck orthon */

      SUBROUTINE ORTHON(IPRSOS,N,NVEC,CX,CY,IFLAG)
C
C      Stephan P. A. Sauer 18/12-1999
C
C
C     The subroutine orthogonalizes one vector of length N 
C     stored in the array CY(N) to NVEC vectors stored in 
C     CX(N,NVEC) sequentially
C
C     The vector in CX-array must be normalized.
C     The vector in CY-array is normalized.
C
#include "implicit.h"
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "infind.h"
C
      DIMENSION CX(N,NVEC)
      DIMENSION CY(N)
C
      PARAMETER ( D1 = 1.0D0, D10M6 = 1.0D-6 )
C
      CALL QENTER('ORTHON')
C
      IFLAG = 0
      ITAB = 1
C
      IF (IPRSOS .GT. 10) THEN
         WRITE(LUPRI,'(/A)')' CLOCV(N,NVEC): in ORTHO before orthonor.'
         CALL OUTPUT(CX,1,N,1,NVEC,N,NVEC,1,LUPRI)
         WRITE(LUPRI,'(/A)')' CLOCV(N,1): in ORTHO before orthonor.'
         CALL OUTPUT(CY,1,N,1,1,N,1,1,LUPRI)
      END IF
C
           DO  NV = 1, NVEC
C
              ATMP =  DDOT(N,CX(1,NV),1,CY,1)
                DO  K = 1, N
                   CY(K) = CY(K) - ATMP*CX(K,NV)
                END DO
C
           END DO 
C
      ATMP =  DDOT(N,CY,1,CY,1)
      IF (ATMP .LT. D10M6) IFLAG = 1
      ATMP = D1/SQRT(ATMP)
           DO  K = 1, N
           CY(K) = CY(K)*ATMP
           END DO
C
      IF (IPRSOS .GT. 10) THEN
         WRITE(LUPRI,'(/A)')' CLOCV(N,NVEC+1): in ORTHO after Orthonor.'
         CALL OUTPUT(CX,1,N,1,NVEC+1,N,NVEC+1,1,LUPRI)
      END IF
C
      CALL QEXIT('ORTHON')
C
      RETURN
      END

C***********************************************************************
C  /* Deck mullloc */
      SUBROUTINE MULLLOC(IPRSOS,CLOCO,THRLOC,PAST,SS,BIJ,AIJ,NLOCC)
C
C      Localize molecular orbitals with  Mulliken population functional
C 
#include "implicit.h"
C
C LUPRI used from include file priunit.h
C MAXORB is needed in COMMON /INFIND/
C MAXASH is needed in COMMON /INFIND/
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
C
C NSYM, NOCCT, NORBT, NRHF(8), NVIR(8) used from COMMON /INFORB/
C ISX(MAXORB) used from COMMON /INFIND/
#include "priunit.h"
#include "inforb.h"
#include "infind.h" 
#include "nuclei.h"
C
      PARAMETER (D0 = 0.0D+00, D1 = 1.0D+00, D2 = 2.0D+00)
      PARAMETER (DP25 = 0.25D+00, DP5 = 0.5D0)
C 
      DIMENSION   CLOCO(NOCCT,NLOCC)
      DIMENSION   SS(NOCCT,NOCCT)
      DIMENSION   PAST(NOCCT,NOCCT,NUCDEP)
      DIMENSION   AIJ(NOCCT,NOCCT),BIJ(NOCCT,NOCCT)
C  
      LOGICAL CONV
C
      CALL QENTER('MULLLOC')
C
C-----------------------------------------
C     Some definitions and initializations
C-----------------------------------------
C
      THETA = D0
      THR2  = THRLOC**2
      PID4  = ATAN(D1)
      PID8  = DP5*PID4
C
      CALL DZERO(BIJ,NOCCT*NOCCT)
      CALL DZERO(AIJ,NOCCT*NOCCT)
      CALL DZERO(SS,NOCCT*NOCCT)
      CALL DUNIT(SS,NOCCT)
C
C--------------------------------
C     Generation of  AIJ BIJ BMAX
C--------------------------------
C
      CALL GENERA_IJ(PAST,AIJ,BIJ)
      call GENERA_BMAX(BIJ,AIJ,BMAX)  
C
  300 BMAX = BMAX/NOCCT
      IF (BMAX .LT. THRLOC) BMAX = THRLOC
  400 CONV = .TRUE.
C
      DO I=2,NOCCT
         DO J = 1, I-1
            ARS = AIJ(I,J)
            BRS = BIJ(I,J)
            IF (ABS(ARS) .GE. THR2) THEN 
               THETA = DP25*ATAN(BRS/ARS)        
               IF(ARS .LT. D0) THEN
                  THETA = THETA + SIGN(PID4,BRS) 
                  IF (BRS .EQ. D0) THETA = PID4    
               END IF
            ELSE
C
C-----------------------
C     Now the case Ars=0 
C-----------------------
C
               THETA = SIGN(PID8,BRS)
            END IF
C
            DT = ABS(THETA)
C            
            IF (ABS(DT-PID4) .LE. THRLOC) THETA = SIGN(PID4,THETA) 
            IF (ABS(DT-PID8) .LE. THRLOC) THETA = SIGN(PID8,THETA)
C
            CALL LOC_ROT(I,J,THETA,NLOCC,SS,BIJ,AIJ,PAST)
C
         END DO
      END DO
C
      IF( .NOT. CONV )      GOTO 400
C      
      IF ( BMAX .GT. THRLOC ) GOTO 300
C
C---------------------------------------------------
C     Convergence reached
C     Copying the transformation matrix to the CLOCO
C---------------------------------------------------
C      
      DO IOCCT = 1, NOCCT
         DO JOCCT = 1, NOCCT
            CLOCO(IOCCT,JOCCT) = SS(IOCCT,JOCCT)
         END DO
      END DO
C
      CALL QEXIT('MULLLOC')
C
      RETURN
      END

C***********************************************************************
C  /* Deck genera_ij */
      SUBROUTINE GENERA_IJ(PAST,AIJ,BIJ)
C
C----------------------------------------------------------------
C     This subroutine generates sum over al the atoms calculating
C     Aij Bij acording with formulas 29a 29b
C----------------------------------------------------------------
C      
#include "implicit.h"
C     
C     LUPRI used from include file priunit.h
C     MAXORB is needed in COMMON /INFIND/
C     MAXASH is needed in COMMON /INFIND/
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
C     NSYM, NOCCT, NORBT, NRHF(8), NVIR(8) used from COMMON /INFORB/
C     ISX(MAXORB) used from COMMON /INFIND/
#include "priunit.h"
#include "inforb.h"
#include "infind.h"
#include "nuclei.h"
C
      PARAMETER (D0 = 0.0D+00, D1 = 1.0D+00, D2 = 2.0D+00)
      PARAMETER (DP25 = 0.25D+00)
C       
      DIMENSION   PAST(NOCCT,NOCCT,NUCDEP)
      DIMENSION   AIJ(NOCCT,NOCCT),BIJ(NOCCT,NOCCT)
C
      CALL QENTER('GENERA_IJ')
C
      CALL DZERO(BIJ,NOCCT*NOCCT)
      CALL DZERO(AIJ,NOCCT*NOCCT)
C
      DO I = 1, NOCCT
         DO J = 1, I
            DO IATOM = 1, NUCDEP
C              
               AIJ(I,J) = AIJ(I,J) + DP25*(PAST(I,I,IATOM) - 
     &             PAST(J,J,IATOM))**2 - PAST(I,J,IATOM)**2
               BIJ(I,J) = BIJ(I,J) + 
     &             PAST(I,J,IATOM)*(PAST(I,I,IATOM) -
     &             PAST(j,j,iatom))
C
            END DO
C
            IF (I .NE. J) THEN
               BIJ(J,I) = -BIJ(I,J)
               AIJ(J,I) = -AIJ(I,J)
            END IF
C            
         END DO
      END DO
C
      CALL QENTER('GENERA_IJ')
C
      RETURN
      END

C***********************************************************************
C  /* Deck genera_bmax */
      SUBROUTINE GENERA_BMAX(BIJ,AIJ,BMAX)   
C
C
C     This subroutine calculates DMAX
C    (some kind of threshold for convergence)
C
C
#include "implicit.h"
C
C LUPRI used from include file priunit.h
C MAXORB is needed in COMMON /INFIND/
C MAXASH is needed in COMMON /INFIND/
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
C
C NSYM, NOCCT, NORBT, NRHF(8), NVIR(8) used from COMMON /INFORB/
C ISX(MAXORB) used from COMMON /INFIND/
#include "priunit.h"
#include "inforb.h"
#include "infind.h" 
#include "nuclei.h"
C   
      DIMENSION   AIJ(NOCCT,NOCCT),BIJ(NOCCT,NOCCT)
C
      CALL QENTER('GENERA_BMAX')
C
      DO I = 1, NOCCT
         DO J = 1, I-1
            BMAX = BMAX + SQRT(AIJ(I,J)*AIJ(I,J) + 
     &             BIJ(I,J)*BIJ(I,J)) - AIJ(I,J)
         END DO
      END DO
C
      CALL QEXIT('GENERA_BMAX')
C
      RETURN
      END

C***********************************************************************
C  /* Deck genera_bmax */
      SUBROUTINE LOC_ROT(I,J,THETA,NLOC,SS,BIJ,AIJ,PAST)
C
C
C     This subroutine does the actualization of the orbitals 
C
C
#include "implicit.h"
C
C LUPRI used from include file priunit.h
C MAXORB is needed in COMMON /INFIND/
C MAXASH is needed in COMMON /INFIND/
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
C
C NSYM, NOCCT, NORBT, NRHF(8), NVIR(8) used from COMMON /INFORB/
C ISX(MAXORB) used from COMMON /INFIND/
#include "priunit.h"
#include "inforb.h"
#include "infind.h" 
#include "nuclei.h"
C
      PARAMETER (D0 = 0.0D+00, D1 = 1.0D+00, D2 = 2.0D+00)
      PARAMETER (DP25 = 0.25D+00)
C
      DIMENSION   SS(NOCCT,NOCCT)
      DIMENSION   PAST(NOCCT,NOCCT,NUCDEP)
      DIMENSION   AIJ(NOCCT,NOCCT),BIJ(NOCCT,NOCCT) 
C
      CALL QENTER('LOC_ROT')
C
      CO = COS(THETA)
      SI = SIN(THETA)
C
C--------------------------------
C     Updating the transformation
C--------------------------------
C
      DO KK = 1, NOCCT
         TI = SS(KK,I)
         TJ = SS(KK,J)
         SS(KK,I) = CO*TI + SI*TJ
         SS(KK,J) = CO*TJ - SI*TI
      END DO
C
C----------------------------
C     Updates matrix elements
C----------------------------
C
      DO IA = 1, NUCDEP
         DO JJ = 1, NOCCT
            TI = PAST(JJ,I,IA)
            TJ = PAST(JJ,J,IA)
            PAST(JJ,I,IA) = CO*TI + SI*TJ
            PAST(JJ,J,IA) = CO*TJ - SI*TI
         END DO   
C
         DO JJ = 1, NOCCT
            TI = PAST(I,JJ,IA)
            TJ = PAST(J,JJ,IA)
            PAST(I,JJ,IA) = CO*TI + SI*TJ
            PAST(J,JJ,IA) = CO*TJ - SI*TI
         END DO   
      END DO
C
      CALL GENERA_IJ(PAST,AIJ,BIJ)
C
      CALL QEXIT('LOC_ROT')
C
      RETURN
      END

C***********************************************************************
C  /* Deck newmolden */
#ifdef NOT_FINISHED_YET
        SUBROUTINE NEWMOLDEN(IPRSOS,CLMO,OCCUP)
C
C
C       This subroutine copies the molden.inp ( with nonlocalized mo's )
C       to molden1.inp ( with localized orbitals )
C
C
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "inftap.h"
#include "maxorb.h"
#include "molde.h"
       PARAMETER( N1 = 4, N2 = 7)
       DIMENSION CLMO(*),OCCUP(*)       
CPFP
C      INTEGER*4 UNIT2, MATCH1, MATCH2, NRF1, NRF2, NOCCT, NBAST
       INTEGER   UNIT2, MATCH1, MATCH2, NRF1, NRF2
Cend-PFP
       CHARACTER LEE(77), WORD1(n1), WORD2(n2)
C
      CALL QENTER('NEWMOLDEN')
C
       word1(1)='['
       word1(2)='M'
       word1(3)='O'
       word1(4)=']'
C
       word2(1)='['
       word2(2)='A'
       word2(3)='t'
       word2(4)='o'
       word2(5)='m'
       word2(6)='s'
       word2(7)=']'
C
       call flshfo(lupri) 
C
       LU_MOLDENLOC = -1
       CALL GPOPEN(LU_MOLDENLOC,'moldenloc.inp','UNKNOWN',' ',
     &             'FORMATTED',IDUMMY,.FALSE.)
C
       call flshfo(lupri) 
C
       REWIND(LUMOLDEN)     
C
 1     READ(LUMOLDEN,'(77A)',err=10,end=10)lee
C
       IF ( lee(1) .EQ. '[') THEN
          match1 = 0
          CALL EQUALLY(lee,word1,n1,match1)
          IF ( match1 .EQ. 4 ) THEN
C
C--------------------------------
C    here we should call to MOMOS
C--------------------------------
C
             DONEIV = .FALSE.
             CALL lmomos(LU_MOLDENLOC,1,CLMO,OCCUP)
 2           READ (LUMOLDEN,'(77A)',err=10,end=10)lee
             IF ( lee(1) .eq. '[' ) THEN
                match2 = 0
                CALL EQUALLY(lee,word2,n2,match2)
                IF ( match2 .EQ. 7 ) THEN
                   WRITE(LU_MOLDENLOC,'(77A)')lee
                   GOTO 1
                ELSE
                   GOTO 2
                END IF
             ELSE
                GOTO 2
             END IF
          ELSE
             WRITE(LU_MOLDENLOC,'(77A)')lee
             GOTO 1
          END IF
       ELSE
          WRITE(LU_MOLDENLOC,'(77A)')lee
          GOTO 1
       END IF
C
  10   CONTINUE
C
       WRITE(LUPRI,'(/,2X,A,/,2X,A,/,2X,A,/,2X,A,/)')
     &        '----------------------------------------------------',
     &        '---------- Molden (loc) input generated ------------',
     &        '--- Use [-v] opt. when running dalton to see it ----',
     &        '----------------------------------------------------'
C
       CALL GPCLOSE(LU_MOLDENLOC,'KEEP')
C
      CALL QEXIT('NEWMOLDEN')
C
       RETURN
       END
#endif
C***********************************************************************
C  /* Deck equally */
       SUBROUTINE EQUALLY(lee,word1,n1,match1)
C
C      Compare to characters
C
       CHARACTER LEE,WORD1
       INTEGER*4 N1,MATCH1
       DIMENSION LEE(N1),WORD1(N1)
C
      CALL QENTER('EQUALLY')
C
       DO INDEX1 = 1, N1
          IF ( LEE(INDEX1) .EQ. WORD1(INDEX1) ) THEN
             MATCH1 = MATCH1 + 1
          END IF
       END DO
C
      CALL QEXIT('EQUALLY')
C
       RETURN 
       END 


C***********************************************************************
C  /* Deck lmomos */
      SUBROUTINE LMOMOS(LU_MOLDENLOC,ITASK,ORVAL,OCCUP)
C
C     This is just momos in for localized orbtitals
C
C ORVAL  =  contains MO coefficients
C itask  =  1 print everything to file
C           2 save orbital energies in ORVAL
C
#include "implicit.h"
#include "priunit.h"
ckr#include "bassel.h"
#include "cbieri.h"
#include "inforb.h"
#include "inftap.h"
#include "maxorb.h"
C
#include "maxaqn.h"
#include "mxcent.h"
#include "aosotr.h"
#include "chrsgn.h"
C
#include "nuclei.h"
#include "symmet.h"
#include "molde.h"
C    
      LOGICAL WRTELEM,WRTZERO
      DIMENSION ORVAL(*), OCCUP(*)
C
      CALL QENTER('LMOMOS')
C
      DONEIV = .FALSE.
C
      IF (ITASK .EQ. 1) THEN
C
         IF (.NOT. DONEIV) THEN 
C
          WRITE(LU_MOLDENLOC,'(/A)') '[MO]'
C          
            ISYMCLAS = 1
            ISYMORB = 0
            IADD = 0
            ICMMO = 1
            DO 1 I=1,NBAST
C
               IF (ISYMORB .EQ. NAOS(ISYMCLAS)) THEN
                  ISYMORB = 0
                  IADD = IADD + NAOS(ISYMCLAS)
                  ISYMCLAS = ISYMCLAS + 1
               END IF
C
               WRITE(LU_MOLDENLOC,'(A,1X,F9.4)') 'Ene=',OREN(I)
               WRITE(LU_MOLDENLOC,'(A,1X)') 'Spin= Alpha'
               WRITE(LU_MOLDENLOC,'(A,1X,F7.4)') 'Occup=',OCCUP(I)
C
               DO 4 M=1,NBAST
C
                  WRTELEM = .FALSE.
                  WRTZERO = .TRUE.
C
                  DO 2 K=1,NAOS(ISYMCLAS)
C
                     DO 3 J=1,NUCDEG(IPCEN(K+IADD))
                        IF (M.EQ.ITRAN(K+IADD,J))THEN
C
                           IF (J.EQ.NUCDEG(IPCEN(K+IADD)))THEN
                              WRTELEM = .TRUE.
                           ELSE 
                              WRTZERO = .FALSE.
                           END IF
C
                           IF (CHRSGN(NINT(CTRAN(K+IADD,J))).EQ.'+')THEN
                              WRITE(LU_MOLDENLOC,'(1X,I2,1X,F10.6)')
     &                             M,ORVAL(ICMMO)
                           ELSE
                             WRITE(LU_MOLDENLOC,'(1X,I2,1X,F10.6)')
     &                             M,-1.0D0*ORVAL(ICMMO)
                           END IF
C
                        END IF
 3                   CONTINUE
 2                CONTINUE
C
                  IF(WRTELEM) THEN
                     ICMMO = ICMMO+1
                  ELSE IF(WRTZERO) THEN
                     WRITE(LU_MOLDENLOC,'(1X,I2,1X,A)')M,'  0.000000'
                  END IF
C
 4             CONTINUE
C
               ISYMORB=ISYMORB+1
 1          CONTINUE
C
            CALL LMOATOMS('ATOM',LU_MOLDENLOC) 
C
         END IF
         DONEIV = .TRUE.
C
      END IF
C
      IF (ITASK .EQ. 2) THEN
C
         DO I = 1, NBAST
            OREN(I)=ORVAL(I)
         END DO
C
      END IF
C
      CALL QEXIT('LMOMOS')
C
      RETURN 
      END

C***********************************************************************
C  /* Deck lmoatoms */
      SUBROUTINE LMOATOMS(WORD,LU_MOLDENLOC)
C
#include "implicit.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"

      CHARACTER*4 NAME
      CHARACTER*4 WORD
CPFP with pgi77
#include "molde.h"
Cend-PFP
#include "priunit.h"
#include "nuclei.h"
#include "symmet.h"
#include "pgroup.h"
#include "cbirea.h"
#include "inftap.h"
#include "chrxyz.h"
#include "chrsgn.h"
#include "chrnos.h"

#include "codata.h"
C from codata: XTANG
CPFP with g77
C#include "molde.h"
Cend-PFP
C
C      CHARACTER*2 ASYMB
C
C      DATA (ASYMB(I),I = 1,103)
C     1/'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
C     2 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
C     3 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
C     4 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
C     5 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
C     6 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
C     7 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
C     8 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
C     9 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',

C     O 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
C     1 'Md', 'No', 'Lr' /
CPFP
C       IF (.NOT. DONEIU .OR. WORD .EQ. 'XYZ' .OR. WORD .EQ. 'FREQ') THEN
C
C          NCOOR = 3*NUCDEP
C          IF (WORD .EQ. 'XYZ ') WRITE(LU_MOLDENLOC,'(I3/)') NUCIND 
C          IF (WORD .EQ. 'ATOM') WRITE(LU_MOLDENLOC,'(A)') '[Atoms] AU' 
C          IF(WORD .EQ. 'FREQ') WRITE(LU_MOLDENLOC,'(A)') ' [FR-COORD]' 
C
C          ICRX = 1
c          ICRY = 2
C          ICRZ = 3
C          DO 100 ICENT = 1, NUCIND
C             MULCNT = ISTBNU(ICENT)
C             NAME   = NAMEX(3*ICENT)(1:4)
C             IF (MULT(MULCNT) .EQ. 1) THEN
C                ICHARGE = NINT(CHARGE(ICENT))
C                   
C                IF (WORD .EQ. 'ATOM')
C     &               WRITE (LU_MOLDENLOC,'(A,1X,I5,1X,I5,3(1X,F20.10))')
C     &               NAME,ICENT,ICHARGE,(CORD(K,ICENT),K=1,3)
C
C                IF(WORD .EQ. 'FREQ' .OR. WORD .EQ. 'XYZ ' )
C     &               WRITE (LU_MOLDENLOC,'(A,3(5X,F15.10))')
C     &               NAME,(CORD(K,ICENT),K=1,3)
C
C                ICRX = ICRX + 3
C                ICRY = ICRY + 3
C                ICRZ = ICRZ + 3
C             ELSE
C                JATOM = 0
C                DO 200 ISYMOP = 0, MAXOPR
C                   IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
C                      JATOM = JATOM + 1
C                      ICHARGE = NINT(CHARGE(ICENT))
C
C                      IF (WORD .EQ. 'ATOM')
C     &                  WRITE (LU_MOLDENLOC,'(A,1X,I5,1X,I5,3(1X,F20.10))')
C     &                         NAME,ICENT+JATOM-1,ICHARGE,
C     &           ((PT(IAND(ISYMAX(K,1),ISYMOP))*CORD(K,ICENT)),K=1,3)
C                      IF(WORD .EQ. 'FREQ')
C     &                   WRITE (LU_MOLDENLOC,'(A,3(1X,F20.10))') NAME,
C     &           ((PT(IAND(ISYMAX(K,1),ISYMOP))*CORD(K,ICENT)),K=1,3)
C                      ICRX = ICRX + 3
C                      ICRY = ICRY + 3
C                      ICRZ = ICRZ + 3
C                   END IF
C 200            CONTINUE
C             END IF        
C 100      CONTINUE
C
C       END IF
C       DONEIU = .TRUE.
C       RETURN
C       END

       IF (.NOT. DONEIU .OR. WORD .EQ. 'XYZ' .OR. WORD .EQ. 'FREQ') THEN

          NCOOR = 3*NUCDEP
          IF (WORD .EQ. 'XYZ ') WRITE(LU_MOLDENLOC,'(I5/)') NUCDEP 
          IF (WORD .EQ. 'ATOM') WRITE(LU_MOLDENLOC,'(A)') '[Atoms] AU' 
          IF (WORD .EQ. 'FREQ') WRITE(LU_MOLDENLOC,'(/A)')'[FR-COORD]' 
       
          IATOM = 0
          DO 100 ICENT = 1, NUCIND
             MULCNT = ISTBNU(ICENT)
             NAME   = '      '
             J = 0
             DO I = 1,4
                IF (NAMN(ICENT)(I:I) .NE. ' ') THEN
                   J = J + 1
                   NAME(J:J) = NAMN(ICENT)(I:I)
                END IF
             END DO
             IF (MULT(MULCNT) .EQ. 1) THEN
                IF (WORD .EQ. 'ATOM') THEN
                   ICHARGE = NINT(CHARGE(ICENT))
C FIXME ?          hjaaj Oct 2003: should we remove point charges ??
                   IATOM = IATOM + 1
                   WRITE (LU_MOLDENLOC,'(A,2I6,3F21.10)')
     &                    NAME,IATOM,ICHARGE,(CORD(K,ICENT),K=1,3)
                ELSE IF (WORD .EQ. 'FREQ') THEN
                   WRITE (LU_MOLDENLOC,'(A,3F21.10)')
     &                    NAME,(CORD(K,ICENT),K=1,3)
                ELSE IF (WORD .EQ. 'XYZ ') THEN
                   WRITE (LU_MOLDENLOC,'(A,3F21.10)')
     &                    NAME,(XTANG*CORD(K,ICENT),K=1,3)
                END IF
             ELSE
                JATOM = 0
                J = J + 1
                NAME(J:J) = '_'
                J = J + 1
                DO 200 ISYMOP = 0, MAXOPR
                   IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
                      JATOM = JATOM + 1
                      NAME(J:J) = CHRNOS(JATOM)
                      CRX = PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
                      CRY = PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
                      CRZ = PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
                      IF (WORD .EQ. 'ATOM') THEN
                         ICHARGE = NINT(CHARGE(ICENT))
                         IATOM = IATOM + 1
                         WRITE (LU_MOLDENLOC,'(A,2I6,3F21.10)')
     &                        NAME,IATOM,ICHARGE,CRX,CRY,CRZ
                      ELSE IF (WORD .EQ. 'FREQ') THEN
                         WRITE (LU_MOLDENLOC,'(A,3F21.10)')
     &                        NAME,CRX,CRY,CRZ
                      ELSE IF (WORD .EQ. 'XYZ ') THEN
                         WRITE (LU_MOLDENLOC,'(A,3F21.10)')
     &                        NAME,XTANG*CRX,XTANG*CRY,XTANG*CRZ
                      END IF
                   END IF
 200            CONTINUE
             END IF        
 100      CONTINUE

       END IF
       DONEIU = .TRUE.
       CALL FLSHFO(LU_MOLDENLOC)
       RETURN
       END
