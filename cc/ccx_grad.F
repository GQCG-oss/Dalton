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
*=====================================================================*
      SUBROUTINE CCX_GRAD(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: calculate relaxed excited state expectation values
*             of one- and two- electron operators the dummy way,
*             i.e., setting up the O1 vectors and some other vectors
*             in the MO basis and calculating LE A^(1) RE.
*
*        
*     Written by Christof Haettig, July 2002
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccexgr.h>
#include <ccgr.h>
#include "maxorb.h"
#include "cclists.h"
#include "mxcent.h"
#include "nuclei.h"
#include "energy.h"
#include "taymol.h"
#include "ccroper.h"
#include "ccr1rsp.h"
#include "cco1rsp.h"
#include "ccsdinp.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .TRUE.)
 
      INTEGER LWORK
      INTEGER IXETRAN(MXDIM_XEVEC,NAXGRO*NXGRST)
      INTEGER IOTRAN(NXGRST)
      INTEGER IEDOTS(NAXGRO*NXGRST), IODOTS(NAXGRO,NXGRST)

#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL ECONS(NAXGRO*NXGRST), OCONS(NAXGRO,NXGRST)
      REAL RDUM, PROPT, PROPE, PROPN, GRDNRM
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION ECONS(NAXGRO*NXGRST), OCONS(NAXGRO,NXGRST)
      DOUBLE PRECISION RDUM, PROPT, PROPE, PROPN, GRDNRM
#endif 

      CHARACTER CDUM*1
      INTEGER NXETRAN, IOP, IOPER, IRELAX, IOPTRES, IORDER, IDUM, IX,
     &        IDXO1, ISYST, IEX, ISCOOR

      INTEGER IR1TAMP, IRHSR1, ILSTSYM, IDXSYM
#if defined (SYS_CRAY) 
      REAL CC_NUCCON, DDOT
#else 
      DOUBLE PRECISION CC_NUCCON, DDOT
#endif 

      CALL QENTER('CCX_GRAD')

      IF (LOCDBG) WRITE(LUPRI,*) 'entered CCX_GRAD...'


      DO IX = 1, NXGRST
        IOTRAN(IX) = IX
        DO IOP = 1, NAXGRO
          IODOTS(IOP,IX) = 0
          OCONS(IOP,IX)  = 0.0D0
        END DO
      END DO
   
      NXETRAN = 0
      DO IOP = 1, NAXGRO
        IOPER = IAXGRO(IOP)
        IF (LPDBSOP(IOPER) .AND. ISYOPR(IOPER).EQ.1) THEN
         IRELAX = IR1TAMP(LBLOPR(IOPER),.TRUE.,0.0D0,ISYOPR(IOPER))
         IDXO1  = IRHSR1(LBLOPR(IOPER),.TRUE.,0.0D0,ISYOPR(IOPER))
cch
          write(lupri,*) 'ccx_grad> IOP,IOPER,LBLOPR(IOPER):',
     &                              IOP,IOPER,LBLOPR(IOPER)
          write(lupri,*) 'ccx_grad> IRELAX:',IRELAX
          write(lupri,*) 'ccx_grad> IDXO1 :',IDXO1
cch
         DO IX = 1, NXGRST
cch
          write(lupri,*) 'ccx_grad> IX,IXGRST(IX):',IX,IXGRST(IX)
cch
          NXETRAN = NXETRAN + 1
          IXETRAN(1,NXETRAN) = IOPER
          IXETRAN(2,NXETRAN) = IXGRST(IX)
          IXETRAN(3,NXETRAN) = -1
          IXETRAN(4,NXETRAN) =  0
          IXETRAN(5,NXETRAN) = IRELAX
          IXETRAN(6,NXETRAN) = -1
          IXETRAN(7,NXETRAN) = -1
          IXETRAN(8,NXETRAN) = -1
          IEDOTS(NXETRAN)    = IXGRST(IX)
          IOTRAN(IX)         = IX
          IODOTS(IOP,IX)     = IDXO1
         END DO
        END IF
      END DO

      ! check if anything at all to do...
      IF (NXETRAN.EQ.0) THEN
        CALL QEXIT('CCX_GRAD')
        RETURN
      END IF

      CALL AROUND(' Relaxed and/or two-electron perturbations ')
      WRITE(LUPRI,'(/10x,A)')
     & 'Operator   Total       Electronic  Nuclear'

      IF (.NOT.(CCS.OR.CC2.OR.CCSD)) THEN
        WRITE(LUPRI,*) 'Requested excited states properties not     '
        WRITE(LUPRI,*) 'available for the current wavefunction model'
      END IF

      IF (LOCDBG) WRITE(LUPRI,*) 'call cc_xieta...'

      IF (XGROPT) THEN
        ! set all gradient components to zero
        CALL DZERO(GRADKE,MXCOOR)
        CALL DZERO(GRADNA,MXCOOR)
        CALL DZERO(GRADEE,MXCOOR)
        CALL DZERO(GRADNN,MXCOOR)
        CALL DZERO(GRADFS,MXCOOR)
      END IF

      CALL DZERO(ECONS,NXETRAN)
      CALL DZERO(OCONS,NAXGRO*NXGRST)

      IOPTRES = 5
      IORDER  = 1
      CALL CC_XIETA(IXETRAN,NXETRAN,IOPTRES, IORDER, 'LE ',
     &              CDUM, IDUM, RDUM, 'RE ', IEDOTS, ECONS,
     &              1, WORK,LWORK) 

      
      IF (LOCDBG) WRITE(LUPRI,*) 'returned from cc_xieta...'
      IF (LOCDBG) WRITE(LUPRI,*) 'call cc_dotdrv...'

      CALL CC_DOTDRV('E0 ','O1 ',NXGRST,NAXGRO,IOTRAN,IODOTS,OCONS,
     &               WORK,LWORK)
      
      IF (LOCDBG) WRITE(LUPRI,*) 'returned from cc_dotdrv...'


      ! calculate nuclear repulsion contribution to gradient
      CALL NUCREP(WORK,WORK(MXCOOR*MXCOOR+1),WORK(2*MXCOOR*MXCOOR+1))

      NXETRAN = 0
      DO IOP = 1, NAXGRO
        IOPER = IAXGRO(IOP)
        IF (LPDBSOP(IOPER)) THEN
         IF (ISYOPR(IOPER).EQ.1) THEN
           IRELAX = IR1TAMP(LBLOPR(IOPER),.TRUE.,0.0D0,ISYOPR(IOPER))
           IDXO1  = IRHSR1(LBLOPR(IOPER),.TRUE.,0.0D0,ISYOPR(IOPER))
           PROPN  = CC_NUCCON(LBLOPR(IOPER),ISYOPR(IOPER))
         ELSE
           PROPN  = 0.0D0
         END IF
         DO IX = 1, NXGRST
           NXETRAN = NXETRAN + 1

           IF (ISYOPR(IOPER).EQ.1) THEN
             PROPE = ECONS(NXETRAN)+OCONS(IOP,IX)+AVEO1(IDXO1)
             PROPT = PROPN + PROPE
           ELSE
             PROPE = 0.0D0
             PROPT = PROPN + PROPE
           END IF

           ISYST = ILSTSYM('LE ',IXGRST(IX))
           IEX   = IDXSYM('LE ',ISYST,IXGRST(IX))

           IF (XGROPT .AND. ISYST.EQ.IXSTSY .AND. IEX.EQ.IXSTAT .AND.
     &         LBLOPR(IOPER)(1:5).EQ.'1DHAM') THEN
             READ(LBLOPR(IOPER)(6:8),'(I3)') ISCOOR
             GRADEE(ISCOOR) = PROPE
           END IF

           IF (ISYOPR(IOPER).EQ.1) THEN
            WRITE(LUPRI,'(10x,A9,2x,F10.6,2x,F10.6,2x,F10.6)')
     &        LBLOPR(IOPER)//':',PROPT,PROPE,PROPN
            IF (LOCDBG) THEN
             WRITE(LUPRI,'(20X,A9,1X,F10.6)') LBLOPR(IOPER)//':', PROPT
             WRITE(LUPRI,'(20X,A9,1X,F10.6)')'E A^(1) E:',ECONS(NXETRAN)
             WRITE(LUPRI,'(20X,A9,1X,F10.6)')'E0 O1    :',OCONS(IOP,IX)
             WRITE(LUPRI,'(20X,A9,1X,F10.6)')'AVERAGE  :',AVEO1(IDXO1)
             WRITE(LUPRI,'(20X,A9,1X,F10.6)')'NUCCON   :',PROPN
            END IF
           ELSE
             WRITE(LUPRI,'(10X,A9,1X,A)') LBLOPR(IOPER)//':', 
     &         ' zero by symmetry '
           END IF

         END DO
        END IF
      END DO

      WRITE(LUPRI,*) ' '

      IF (XGROPT) THEN
        IF ( IPRINT.GT.1) THEN
          CALL CC_SETDORPS('1DHAM   ',.TRUE.,IPRINT)
          ! calculate nuclear repulsion contribution to gradient
          CALL NUCREP(WORK,WORK(MXCOOR*MXCOOR+1),
     &                     WORK(2*MXCOOR*MXCOOR+1))
          CALL HEADER('Electronic gradient',-1)
          CALL PRIGRD(GRADEE,WORK,WORK(MXCOOR*MXCOOR+1))
        END IF

        CALL ZERGRD
        CALL ADDGRD(GRADNN)
        CALL ADDGRD(GRADEE)

        CALL HEADER('Molecular gradient',-1)
        CALL PRIGRD(GRDMOL,WORK,WORK(MXCOOR*MXCOOR+1))

        GRDNRM = DDOT(3*NUCDEP,GRDMOL,1,GRDMOL,1)
        WRITE (LUPRI,'(/19X,A,1P,E10.2)')
     *     'Molecular gradient norm:', GRDNRM

      END IF

      
      CALL QEXIT('CCX_GRAD')

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCX_GRAD                             *
*---------------------------------------------------------------------*
