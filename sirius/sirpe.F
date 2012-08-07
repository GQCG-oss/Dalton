      SUBROUTINE pe_gradient(CREF,CMO,INDXCI,DV,G,EQMMM,WRK,LFREE)

      use polarizable_embedding, only: pe_master
C
C
C     Written by Erik Donovan Hedegård (edh) based on PCMGRAD 
C
C     Purpose:  calculate MCSCF energy and gradient contribution
C               from a PE medium
C
C     Output:
C     G          MCSCF gradient with solvation contribution added
C     ESOLT      total solvation energy
C
C Used from common blocks:
C   INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
C   INFORB: NNASHX, NNBASX, NNORBX, etc.
C   INFIND: IROW(*)
C   INFTAP: LUSOL,  LUIT2
C   INFPRI: IPRSOL
C
#include "implicit.h"
#include "priunit.h" 
#include "pi.h"     
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "qmmm.h"
#include "nuclei.h"
#include "orgcom.h"
#include "infvar.h"
#include "inforb.h"
#include "infind.h"
#include "inftap.h"
#include "infpri.h"
C edh 09/11 2011
#include "gnrinf.h"

      DIMENSION CREF(*), CMO(*), INDXCI(*)
      DIMENSION DV(*),   G(*),   WRK(LFREE)
      PARAMETER ( D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0,
     &            DCVAL = D2, FPI = 4.0D0 * PI )
      LOGICAL LOCDEB,FNDLAB,FIRST
      CHARACTER*8 STAR8,SOLVDI,EODATA
      DATA        FIRST/.TRUE./, STAR8/'********'/
      DATA        SOLVDI/'SOLVDIAG'/, EODATA/'EODATA  '/

C
C     Statement functions;
C     define automatic arrays (dynamic core allocation)
C
C
      CALL QENTER('PEGRD') 
C
C     Core allocation
C
      LOCDEB = .FALSE.

      KDENC  = 1
      KDENV  = KDENC  + N2BASX
      KDENT  = KDENV  + N2BASX
      KDENTF = KDENT  + N2BASX
C     -------------------------------
      KFPE   = KDENTF + NNBASX
      KUCMO  = KFPE   + NNBASX
      KFPEMO = KUCMO  + NORBT*NBAST
      KFPEAC = KFPEMO + NNORBX
C     ------------------------------
      KGRDPE = KFPEAC + NNASHX        
      KDIAPE = KGRDPE + NVARH       
C     ------------------------------  
      KWRK1  = KDIAPE + NVAR          
      LWRK1  = LFREE  - KWRK1 

      IF (LWRK1 .LT. 0) CALL ERRWRK('PEGRD',-KWRK1,LWRK1)

C     1. KDENC  : Core (inactive) density matrix from fckden routine
C     2. KDENV  : Valence (active) density matrix 
C     3. KDENT  : Total density matrix (sum DC + DV)
C     4. KDENTF : Folded total density matrix 
C     ---------------------------------------------------------------
C     6. KFPE   : Polarizable Embedded (PE) Tg operator (AO basis)
C     7. KUCMO  : MO coefficients
C     8. KFPEMO : Polarizable Embedded (PE) Tg operator (MO basis)
C     9. KFPEAC : - active part
C     --------------------------------------------------------------
C    10. KGRDPE : Solvent contr. to MCSCF gradient (G) 
C                 - Output from SOLGC and SOLGO 
C    11. KDIAPE : -Output from SOLDIA (what is this??) 


      CALL DZERO(WRK(KDENC),N2BASX)
      CALL DZERO(WRK(KDENV),N2BASX)
      CALL DZERO(WRK(KDENT),N2BASX)
      CALL DZERO(WRK(KDENTF),NNBASX)
      CALL DZERO(WRK(KFPE),NNBASX)
      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KFPEMO),NNORBX)
      CALL DZERO(WRK(KFPEAC),NNASHX)
      CALL DZERO(WRK(KGRDPE),NVARH)
      CALL DZERO(WRK(KDIAPE),NVAR)

C ************* Write statements for debugging ****************
C *************************************************************

      IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) )THEN
         WRITE (LUPRI,'(/A/A,2I10)')
     *        ' >>> PEGRD - gtot (input) - non-zero elements',
     *        '     NCONF, NWOPT =',NCONF,NWOPT
         DO 40 I = 1,NCONF
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *           ' conf #',I,G(I)
 40      CONTINUE
         DO 50 I = NCONF+1,NVAR
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *           ' orb  #',I,G(I)
 50      CONTINUE
      END IF

      IF ( (IPQMMM .GE. 15 .AND. NASHT .GT. 0) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(/A)') ' >>> PEGRD - DV matrix :'
         CALL OUTPAK(DV,NASHT,1,LUPRI)
      END IF

C *************************************************************
C *************************************************************

      CALL FCKDEN((NISHT.GT.0),(NASHT.GT.0),WRK(KDENC),WRK(KDENV),
     &            CMO,DV,WRK(KWRK1),LWRK1)

      CALL DCOPY(N2BASX,WRK(KDENC),1,WRK(KDENT),1)        ! Construct DC dens. matetrix (KDENC)
      CALL DAXPY(N2BASX,1.0D0,WRK(KDENV),1,WRK(KDENT),1)  ! Add valence density matrix DV (DC + DV)

      CALL DGEFSP(NBAST,WRK(KDENT),WRK(KDENTF))           ! Fold total dens. matrix 

      call pe_master('fock', denmats=wrk(kdentf:kdentf+nnbasx-1),
     &               fckmats=wrk(kfpe:kfpe+nnbasx-1), nmats=1,
     &               Epe=wrk(kwrk1:kwrk1),
     &               dalwrk=wrk(kwrk1+1:kwrk1+1+lwrk1))
      eqmmm = wrk(kwrk1)

      IF (LOCDEB) THEN
        WRITE(LUPRI,*) 'KDENTF IN PEGRD BEFORE QMMM_FCK_AO'
        CALL OUTPAK(WRK(KDENTF),NBAST,1,LUPRI)
        CALL DCOPY(NNBASX,WRK(KDENTF),1,WRK(KDENC),1)      
      ENDIF



!      CALL QMMM_FCK_AO(WRK(KFPE),WRK(KDENTF),EQMMM,WRK(KWRK1),LWRK1, ! Gradient routine needs EQMMM
!     &                 IPQMMM)                                       ! PEFCMO should be changed  
                                                                     ! to deliver EQMMM as well
      CALL UPKCMO(CMO,WRK(KUCMO))                                    ! requires call change for other places where 
      CALL UTHU(WRK(KFPE),WRK(KFPEMO),WRK(KUCMO),WRK(KWRK1),         ! PEFCMO
     &              NBAST,NORBT)                     

!      CALL PEFCMO(WRK(KUCMO),WRK(KFPEM),DV,WRK(KWRK1),LWRK1,IPQMMM) ! edh: KFPEM is a temp. variable used to debug 
                                                                    ! and prepare this module to magnus' PE module
      WRITE (LUPRI,'(/A)') ' PE_mo matrix from PEFCMO:'             ! problem is that PEFCMO doesn't calc. EQMMM
      CALL OUTPAK(WRK(KFPEMO),NORBT,1,LUPRI)                         ! and now we get it from QMMM_FCK_AO

      IF (NASHT .GT. 0) THEN
         CALL GETAC2(WRK(KFPEMO),WRK(KFPEAC))             
      END IF 

C
C     Expextation value of FPE

      TFPEMO     = SOLELM(DV,WRK(KFPEAC),WRK(KFPEMO),TFPEAC) 

C
C ************* Write statements for debugging ****************
C *************************************************************

      IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(A,F17.8)')
     *   ' --- FPE expectation value MO :',TFPEMO
         WRITE (LUPRI,'(A,F17.8)')
     *   ' --- active part of FPE    :',TFPEAC
      ENDIF

      IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(/A)') ' PE_ao matrix in PEGRD:'
         CALL OUTPAK(WRK(KFPE),NBAST,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' PE_mo matrix in PEGRD:'
         CALL OUTPAK(WRK(KFPEMO),NORBT,1,LUPRI)
         IF (NASHT .GT. 0) THEN
         WRITE (LUPRI,'(/A)') ' PE_ac matrix in PEGRD:'
         CALL OUTPAK(WRK(KFPEAC),NASHT,1,LUPRI)
         ENDIF
      ENDIF

C *************************************************************
C *************************************************************
C
C ******* edh: SOLGC computes the solvent CI integrals *******
C ******* input: CREF(NCONF)    = CI reference state   *******
C *******        KFPEAC(NNASHX) = Solvent integrals    *******
C *******        TFPEAC         = CREF exp. value      *******
C *******        INDXCI(*)      = CI index             *******
C ******* output: GLMCI(NCONF)  = CI solv. gradient    *******       

      IF (NCONF .GT. 1) THEN
         CALL SOLGC(CREF,WRK(KFPEAC),TFPEAC,WRK(KGRDPE),INDXCI, ! NOTE: Output here is GRDPE (solv. CI PE contribution)
     &              WRK(KWRK1),LWRK1)                           ! edh: SOLGC calc. < u | Fg | 0 > + < 0 | Fg | 0 > c_u
      END IF

      IF (NWOPT .GT. 0) THEN
         CALL SOLGO(DCVAL,DV,WRK(KFPEMO),WRK(KGRDPE+NCONF))     ! edh: SOLGO calc. 2 < 0 | [Ers, Fg] | 0 > 
      END IF

      CALL SOLDIA(TFPEAC,WRK(KFPEAC),INDXCI,                    
     *            WRK(KFPEMO),DV,WRK(KDIAPE),WRK(KWRK1),LWRK1)

      DO 420 I = 0,(NVAR-1)                                    
         WRK(KDIAPE+I) = - WRK(KDIAPE+I)
  420 CONTINUE

C     
C ******************* Orthogonality test **********************
C *************************************************************
C     
      IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(/A)')' >>> PEGRD - grdj1, grdj2, diape, '//
     &                      'diape, cref'
         DO 430 I = 1,NCONF
            WRITE (LUPRI,'(A,I10,3F10.6)') ' conf #',I,
     *            WRK(KDIAPE-1+I),
     *            WRK(KDIAPE-1+I),CREF(I)
  430    CONTINUE
      END IF
C     
       TEST = DDOT(NCONF,CREF,1,WRK(KGRDPE),1)
       IF (ABS(TEST) .GT. 1.D-8) THEN
          NWARN = NWARN + 1
          WRITE (LUPRI,*) ' >>> PEGRD WARNING <<< '
          WRITE (LUPRI,*) ' <CREF | GRAD > =',TEST
       END IF

C *************************************************************
C *************************************************************
     
C     Add PE gradient contribution to MCSCF gradient
C     
      CALL DAXPY(NVARH,D1,WRK(KGRDPE),1,G,1)

      IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(/A/A,2I10)')
     *      ' >>> PEGRD - grdB, gtot (accum) - non-zero grdpe',
     *      '     NCONF, NWOPT =',NCONF,NWOPT
         DO 440 I = 1,NCONF
            IF (WRK(KGRDPE-1+I) .NE. D0)
     *         WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' conf #',I,WRK(KGRDPE-1+I),G(I)
  440    CONTINUE
         DO 450 I = NCONF+1,NVAR
            IF (WRK(KGRDPE-1+I) .NE. D0)
     *         WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' orb  #',I,WRK(KGRDPE-1+I),G(I)
  450    CONTINUE
      END IF
C     
      IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(/A/A,2I10)')
     *      ' >>> PEGRD - gtot (output) - non-zero elements',
     *      '     NCONF, NWOPT =',NCONF,NWOPT
         DO 840 I = 1,NCONF
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *      ' conf #',I,G(I)
  840    CONTINUE
         DO 850 I = NCONF+1,NVAR
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *      ' orb  #',I,G(I)
  850    CONTINUE
      END IF

      IF (LUIT2 .GT. 0) THEN
         NC4 = MAX(NCONF,4)
         NW4 = MAX(NWOPT,4)
         REWIND LUIT2
         IF (FNDLAB(EODATA,LUIT2)) BACKSPACE LUIT2
         WRITE (LUIT2) STAR8,STAR8,STAR8,SOLVDI
         IF (NCONF .GT. 1) CALL WRITT(LUIT2,NC4,WRK(KDIAPE))
         WRITE (LUIT2) STAR8,STAR8,STAR8,EODATA
      END IF

      CALL QEXIT('PEGRD')
C     end of pegrd.
      END


