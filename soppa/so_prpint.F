#include <single.h>
C
C  /* Deck so_prpint */
      SUBROUTINE SO_PRPINT(LTYPE,NLBTOT,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Request that needed property integrals are calculated
C              and written to file LUPROP. Property labels and
C              symmetry are stored in LABAPP and LABSYM.
C
#include <implicit.h>
#include <priunit.h>
C
#include <mxcent.h>
#include <cbiexc.h>
#include <cbilnr.h>
#include <inftap.h>
CSPAS:10/11-2003: is merged with maxorb.h
#include <maxorb.h>
C#include <mxorb.h>
CKeinSPASmehr
#include <maxaqn.h>
#include <symmet.h>
#include <gnrinf.h>
CSPAS:29/03-2006: we need IPRSOP from soppinf.h
#include <soppinf.h>
CKeinSPASmehr
C
      DIMENSION   WORK(LWORK)
      CHARACTER*8 LABEL
      CHARACTER*6 LTYPE
      LOGICAL     LDIPLEN, LLONMAG, LHBDO, LDIPVEL, LANGMOM
      LOGICAL     LSECMOM, LTHETA
CSPAS:23/5-11: second and third moment sum rules
      LOGICAL     LTHIRDM
CKeinSPASmehr
      LOGICAL     EXTST   , OPENED
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_PRPINT')
C
C-------------------------------------------
C     Set print level for integral routines.
C-------------------------------------------
C
      IPRINT = IPRSOP - 10
      IF (IPRINT .LT. 0) IPRINT = 0
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      KDLAB   = 1
      KIDSYM  = KDLAB  + 9*MXCENT
      KIDADR  = KIDSYM + 9*MXCENT
      KEND1   = KIDADR + 9*MXCENT
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_PRPINT',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_PRPINT',' ',KEND1,LWORK)
C
C================================================================
C     Determine the types of property integrals which are needed.
C================================================================
C
      LDIPLEN  = .FALSE.
      LLONMAG  = .FALSE.
      LHBDO    = .FALSE.
      LDIPVEL  = .FALSE.
      LANGMOM  = .FALSE.
      LSECMOM  = .FALSE.
CSPAS:23/5-11: second and third moment sum rules
      LTHIRDM  = .FALSE.
CKeinSPASmehr
      LTHETA   = .FALSE.
C
C-----------------------------------------------
C     For the calculation of excitation energies
C-----------------------------------------------
C
      IF (LTYPE .EQ. 'EXCITA') THEN
C
C-----------------------------------------------------------------------
C        Dipole strength and conventional rotatory strength is set true.
C-----------------------------------------------------------------------
C
         LSECMOM = .TRUE.
C
         DIPSTR  = .TRUE.
         ROTVEL  = .TRUE.
C
C
CSPAS:23/5-11: second and third moment sum rules
         IF (SUMRUL) THEN           ! sum rules, dip, quad, oct.
            LDIPLEN = .TRUE.
            LDIPVEL = .TRUE.
            LSECMOM = .TRUE.
            LTHIRDM = .TRUE.
         END IF
CKeinSPASmehr
         IF ( DIPSTR )              LDIPLEN = .TRUE.
         IF ( DIPSTR )              LDIPVEL = .TRUE.
         IF ( ROTSTR )              LDIPLEN = .TRUE.
         IF ( ROTSTR )              LLONMAG = .TRUE.
cKeld      IF ( ROTSTR .AND. NODIFC ) LHBDO   = .TRUE.
         IF ( ROTVEL )              LDIPVEL = .TRUE.
         IF ( ROTVEL )              LANGMOM = .TRUE.
C
C
C------------------------------------------------------
C     For the calculation of linear response properties
C------------------------------------------------------
C
      ELSE IF (LTYPE .EQ. 'LINEAR') THEN
C
C
         IF ( ALFA .OR. ROAA .OR. ROAG )    LDIPLEN = .TRUE.
         IF ( ROAA )                        LTHETA  = .TRUE.
         IF ( ROAG )                        LLONMAG = .TRUE.
C
      ENDIF
C
C----------------------------------------------------------
C     Initialize the number of types of property integrals.
C----------------------------------------------------------
C
      NLBTOT = 0
C
C------------------------------------------
C     Close AOPROPER if it is already open.
C-------------------------------------------
C
      CALL GPINQ('AOPROPER','EXIST',EXTST)
      IF (EXTST) THEN
         CALL GPINQ('AOPROPER','OPENE',OPENED)
         IF (OPENED) THEN
            INQUIRE (FILE='AOPROPER',NUMBER=LUPROP)
            WRITE(LUPRI,'(2A,I3)') ' SO_PRPINT: file AOPROPER is ',
     &                             'already opened with unit number ',
     &                             LUPROP
            CALL GPCLOSE(LUPROP,'KEEP')
            WRITE(LUPRI,'(A)') ' SO_PRPINT: file AOPROPER was closed'
         ENDIF
      ENDIF
C
C--------------------------------------
C     Determine dipole length integrals.
C---------------------------------------
C
      IF ( LDIPLEN ) THEN
C
         NCOMP  = 0
         NPATOM = 0
C
         CALL GET1IN(DUMMY,'DIPLEN ',NCOMP,WORK(KEND1),LWORK1,
     &               WORK(KDLAB),WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &               DUMMY,IPRINT)
C
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,WORK(KDLAB),WORK(KIDSYM),LABAPP,LABSYM)
C
      END IF
C
C----------------------------------------------------------
C     Registrate that LONMAG integrals should be available.
C----------------------------------------------------------
C
      IF ( LLONMAG ) THEN
C
         CALL LABCOP(1,NLBTOT,'XLONMAG ',ISYMAX(1,2),LABAPP,LABSYM)
         CALL LABCOP(1,NLBTOT,'YLONMAG ',ISYMAX(2,2),LABAPP,LABSYM)
         CALL LABCOP(1,NLBTOT,'ZLONMAG ',ISYMAX(3,2),LABAPP,LABSYM)
C
         NCOMP  = 0
         NPATOM = 0
C
         CALL GET1IN(DUMMY,'ANGMOM ',NCOMP,WORK(KEND1),LWORK1,
     &               WORK(KDLAB),WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &               DUMMY,IPRINT)
C
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,WORK(KDLAB),WORK(KIDSYM),LABAPP,LABSYM)
C
      END IF
C
C-------------------------------------------------------
C     Determine half differentiated B overlap integrals.
C-------------------------------------------------------
C
      IF ( LHBDO ) THEN
C
         NCOMP  = 0
         NPATOM = 0
C
         CALL GET1IN(DUMMY,'HBDO   ',NCOMP,WORK(KEND1),LWORK1,
     &               WORK(KDLAB),WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &               DUMMY,IPRINT)
C
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,WORK(KDLAB),WORK(KIDSYM),LABAPP,LABSYM)
C
      END IF
C
C------------------------------------------
C     Determinte dipole velocity intergals.
C------------------------------------------
C
      IF ( LDIPVEL ) THEN
C
         NCOMP  = 0
         NPATOM = 0
C
         CALL GET1IN(DUMMY,'DIPVEL ',NCOMP,WORK(KEND1),LWORK1,
     &               WORK(KDLAB),WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &               DUMMY,IPRINT)
C
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,WORK(KDLAB),WORK(KIDSYM),LABAPP,LABSYM)
C
      END IF
C
C-------------------------------------------
C     Determinte angular momentum intergals.
C-------------------------------------------
C
      IF ( LANGMOM ) THEN
C
         NCOMP  = 0
         NPATOM = 0
C
         CALL GET1IN(DUMMY,'ANGMOM ',NCOMP,WORK(KEND1),LWORK1,
     &               WORK(KDLAB),WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &               DUMMY,IPRINT)
C
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,WORK(KDLAB),WORK(KIDSYM),LABAPP,LABSYM)
C
      END IF
C
C-------------------------------------------
C     Determine quadrupole moment integrals.
C-------------------------------------------
C
      IF ( LSECMOM ) THEN
C
         NCOMP  = 0
         NPATOM = 0
C
         CALL GET1IN(DUMMY,'SECMOM ',NCOMP,WORK(KEND1),LWORK1,
     &               WORK(KDLAB),WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &               DUMMY,IPRINT)
C
         NLAB = 6
         CALL LABCOP(NLAB,NLBTOT,WORK(KDLAB),WORK(KIDSYM),LABAPP,LABSYM)
C 
      END IF
C
C-----------------------------------------------------
C     Determine traceless quadrupole moment integrals.
C-----------------------------------------------------
C
      IF ( LTHETA ) THEN
C
         NCOMP  = 0
         NPATOM = 0
C
         CALL GET1IN(DUMMY,'THETA  ',NCOMP,WORK(KEND1),LWORK1,
     &               WORK(KDLAB),WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &               DUMMY,IPRINT)
C
         NLAB = 6
         CALL LABCOP(NLAB,NLBTOT,WORK(KDLAB),WORK(KIDSYM),LABAPP,LABSYM)
C 
      END IF
CSPAS:23/5-11: second and third moment sum rules
C
C-----------------------------------------------
C     Determine third electric moment integrals.
C-----------------------------------------------
C
      IF ( LTHIRDM ) THEN
C
         NCOMP  = 0
         NPATOM = 0
C
         CALL GET1IN(DUMMY,'THRMOM ',NCOMP,WORK(KEND1),LWORK1,
     &               WORK(KDLAB),WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &               DUMMY,IPRINT)
C
         NLAB = 10
         CALL LABCOP(NLAB,NLBTOT,WORK(KDLAB),WORK(KIDSYM),LABAPP,LABSYM)
C
      END IF
CKeinSPASmehr
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_PRPINT')
C
      RETURN
      END
