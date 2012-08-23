      SUBROUTINE pe_gradient(CREF,CMO,INDXCI,DV,G,EQMMM,WRK,LFREE)

      use polarizable_embedding, only: pe_master
C
C
C     Written by Erik Donovan Hedegård (edh) based on PCMGRAD 
C                with contribution from Magnus Jógvan H.  Olsen

C     Purpose:  calculate MCSCF energy and gradient contribution
C               from a PE medium with PE module
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
      CALL QENTER('pe_gradient') 
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
      KFPEM  = KFPEMO + NNORBX ! = KFPEMO Only for test purpose!
      KFPEAC = KFPEM  + NNORBX
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
      CALL DZERO(WRK(KFPEM),NNORBX) ! Only for test!
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

!      CALL QMMM_FCK_AO(WRK(KFPE),WRK(KDENTF),EQMMM,WRK(KWRK1),LWRK1, ! Gradient routine needs EQMMM
!     &                 IPQMMM)                                       ! PEFCMO should be changed  
                                                                     ! to deliver EQMMM as well
      CALL UPKCMO(CMO,WRK(KUCMO))                                    ! requires call change for other places where 
      CALL UTHU(WRK(KFPE),WRK(KFPEMO),WRK(KUCMO),WRK(KWRK1),         ! PEFCMO
     &              NBAST,NORBT)                     

      CALL pe_fcmo(WRK(KUCMO),WRK(KFPEM),DV,WRK(KWRK1),LWRK1,IPQMMM) 

C edh: KFPEM is a temp. variable used to debug  and prepare this module to magnus' PE module
C problem is that PEFCMO doesn't calc. EQMMM  and now we get it from QMMM_FCK_AO
                                                                    
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
         WRITE (LUPRI,'(/A)') ' PE_ao matrix in KFPE:'
         CALL OUTPAK(WRK(KFPE),NBAST,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' PE_mo matrix in KFPEMO:'
         CALL OUTPAK(WRK(KFPEMO),NORBT,1,LUPRI)
         IF (NASHT .GT. 0) THEN
         WRITE (LUPRI,'(/A)') ' PE_ac matrix in PEGRD:'
         CALL OUTPAK(WRK(KFPEAC),NASHT,1,LUPRI)
         ENDIF
      ENDIF

      IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(/A)') ' PE_ao matrix from pefcmo call in PEGRD:'
         CALL OUTPAK(WRK(KFPE),NBAST,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' PE_mo matrix in KFPEM:'
         CALL OUTPAK(WRK(KFPEM),NORBT,1,LUPRI)
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

      CALL QEXIT('pe_gradient')
C     end of pegrd.
      END
#include "single.h"
C
C edh: single.h is technical. For equalence of
C compiling bewteen CRAY and other machines 
C
C  /* Deck pcmgrd */
      SUBROUTINE pe_fcmo(CMO,FSOL,DV,WRK,LFREE,IPRINT)
C
      use polarizable_embedding, only: pe_master
C
C     Written by Erik Donovan Hedegård (edh)
C
C     Purpose:  Transform (MCSCF) Fg PE operator to MO basis
C               using PE module
C     Output:
C     FSOL      Tg PE operator in MO basis
C
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
C
C
C Used from common blocks:
C   INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
C   INFORB: NNASHX, NNBASX, NNORBX, etc.
C   INFIND: IROW(*)
C   INFTAP: LUSOL,  LUIT2
C   INFPRI: IPRSOL
C

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
#include "gnrinf.h"



      DIMENSION CMO(*), FSOL(*)
      DIMENSION DV(*), WRK(*)
      PARAMETER ( D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0,
     &           DCVAL = D2, FPI = 4.0D0 * PI )

      CALL QENTER('pe_fcmo')

C     Core allocation
C
      KDENC  = 1                      
      KDENV  = KDENC  + N2BASX        
      KDENT  = KDENV  + N2BASX
      KDENTF = KDENT  + N2BASX
C     -------------------------------
      KFPE   = KDENTF + NNBASX
      KUCMO  = KFPE   + NNBASX        
      KFPEMO = KUCMO  + NORBT*NBAST 
C     ------------------------------  
      KWRK1  = KFPEMO + NNORBX
      LWRK1  = LFREE  - KWRK1

C     1. KDENC  : Core (inactive) density matrix. CALL from fckden subroutine
C     2. KDENV  : Valence (active) density matrix 
C     3. KDENT  : Total density matrix (sum DC + DA)
C     4. KDENTF : Folded total density matrix 
C     ------------------------------
C     6. KFPE   : Polarizable Embedded (PE) Tg operator (AO basis)
C     7. KUCMO  : MO coefficients
C     8. KFPEMO : Polarizable Embedded (PE) Tg operator (MO basis)


      CALL DZERO(WRK(KDENC),N2BASX)
      CALL DZERO(WRK(KDENV),N2BASX)
      CALL DZERO(WRK(KDENT),N2BASX)
      CALL DZERO(WRK(KDENTF),NNBASX)
      CALL DZERO(WRK(KFPE),NNBASX) 
      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KFPEMO),NNORBX)

      IF (LWRK1 .LT. 0) CALL ERRWRK('PEFCMO',-KWRK1,LWRK1)

      IF (IPQMMM .GE. 15 .AND. NASHT .GT. 0) THEN
         WRITE (LUPRI,'(/A)') ' >>> PEFCMO - DV matrix :'
         CALL OUTPAK(DV,NASHT,1,LUPRI)
      END IF

      CALL FCKDEN((NISHT.GT.0),(NASHT.GT.0),WRK(KDENC),WRK(KDENV),
     &            CMO,DV,WRK(KWRK1),LWRK1)

      CALL DCOPY(N2BASX,WRK(KDENC),1,WRK(KDENT),1)
      CALL DAXPY(N2BASX,1.0D0,WRK(KDENV),1,WRK(KDENT),1)

      CALL DGEFSP(NBAST,WRK(KDENT),WRK(KDENTF))

      call pe_master('fock', denmats=wrk(kdentf:kdentf+nnbasx-1),
     &               fckmats=wrk(kfpe:kfpe+nnbasx-1), nmats=1,
     &               Epe=wrk(kwrk1:kwrk1),
     &               dalwrk=wrk(kwrk1+1:kwrk1+1+lwrk1))
      eqmmm = wrk(kwrk1)


!      CALL QMMM_FCK_AO(WRK(KFPE),WRK(KDENTF),EQMMM,WRK(KWRK1),LWRK1,
!     &                 IPQMMM)

      CALL UPKCMO(CMO,WRK(KUCMO))
      CALL UTHU(WRK(KFPE),FSOL,WRK(KUCMO),WRK(KWRK1),
     &              NBAST,NORBT)

      IF (IPQMMM .GE. 15) THEN
         WRITE (LUPRI,'(/A)') ' PE_ao matrix in PEFCMO:'
         CALL OUTPAK(WRK(KFPE),NBAST,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' PE_mo matrix in PEFCMO:'
         CALL OUTPAK(FSOL,NORBT,1,LUPRI)
      END IF

      CALL QEXIT('pe_fcmo')
C     end of pefcmo.
      END

C  /* Deck pelin */
      SUBROUTINE pe_lin(NCSIM,NOSIM,BCVECS,BOVECS,CREF,CMO,INDXCI,
     &                  DV,DTV,SCVECS,SOVECS,ORBLIN,WRK,LWRK)          
C
C Written by Erik Donovan Hedegård december 2011 
C after original code by  Hans Joergen Aa. Jensen
C Common driver for pe_lnc and pe_lno
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "inflin.h"
#include "infvar.h"
C edh 13/12 2011
#include "qmmm.h"
#include "gnrinf.h"

C   Used from common blocks:
C   INFLIN : NWOPPT,NVARPT

     
      DIMENSION BCVECS(*),BOVECS(*),CREF(*),CMO(*),INDXCI(*)
      DIMENSION DV(*),DTV(*),SCVECS(*),SOVECS(*),WRK(LWRK)
      LOGICAL   ORBLIN, LOCDEB

      PRINT*, 'IN pe_lin'

      LOCDEB = .FALSE.

      CALL QENTER('pe_lin')

      IF (NCSIM .GT. 0) THEN
         IF ( (LOCDEB) .OR. (IPQMMM.GT.15) ) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED CONFIGURATION VECTOR'
            WRITE(LUPRI,*)' BEFORE pe_lnc CALL, ITERATION # '  
            CALL OUTPUT(SCVECS,1,NCONF,1,NCSIM,NCONF,NCSIM,1,LUPRI)
         END IF
         
         CALL pe_lnc(NCSIM,BCVECS,CREF,CMO,INDXCI,
     &               DV,DTV,SCVECS,WRK,LWRK)

         IF ( (LOCDEB) .OR. (IPQMMM .GT. 15) ) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED CONFIGURATION VECTOR'
            WRITE(LUPRI,*)' AFTER pe_lnc CALL, ITERATION # '
            CALL OUTPUT(SCVECS,1,NCONF,1,NCSIM,NCONF,NCSIM,1,LUPRI)
         END IF
      END IF

      IF ( NOSIM .GT.0 ) THEN
         IF ( .NOT.ORBLIN ) THEN
            NSVAR  = NVARPT
         ELSE
            NSVAR  = NWOPPT
         END IF
         IF ( (LOCDEB) .OR. (IPQMMM .GT. 15) ) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
            WRITE(LUPRI,*)' BEFORE pe_lno CALL, ITERATION # '
            CALL OUTPUT(SOVECS,1,NWOPPT,1,NOSIM,NWOPPT,NOSIM,1,LUPRI)
         END IF

         CALL pe_lno(NOSIM,BOVECS,CREF,CMO,INDXCI,
     &               DV, SOVECS,NSVAR,WRK,LWRK)

         IF ( (LOCDEB) .OR. (IPQMMM .GT. 15) ) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
            WRITE(LUPRI,*)' AFTER pe_lno, ITERATION # '
            CALL OUTPUT(SOVECS,1,NWOPPT,1,NOSIM,NWOPPT,NOSIM,1,LUPRI)
         END IF
      END IF

      CALL QEXIT('pe_lin')
      RETURN
      END

C  /* Deck pelnc */
      SUBROUTINE pe_lnc(NCSIM,BCVEC,CREF,CMO,INDXCI,    
     *                  DV,DTV,SVEC,WRK,LFREE)
C
C  Written by Erik Donovan Hedegaard Jan-03 2012
C  after original routine by Hans Jørgen Aa. Jensen
C
C  Purpose:  Calculate MCSCF Hessian contribution from a
C            surrounding PE medium to a csf trial vector.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "thrzer.h"
#include "maxash.h"
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "inflin.h"
#include "inftap.h"
#include "infpri.h"
#include "qmmm.h"
#include "qm3.h"
#include "gnrinf.h"
#include "orgcom.h"

      DIMENSION BCVEC(*),  CREF(*), CMO(*)
      DIMENSION INDXCI(*), DV(*),   DTV(NNASHX,*)
      DIMENSION SVEC(NVARPT,*),     WRK(*)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LOCDEB, FNDLAB, LPOL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , D2 = 2.0D0 )

      LPOL = .FALSE. 

      CALL QENTER('pe_lnc')
C
C  Used from common blocks:
C    INFORB : NNASHX, NNORBX, NNBASX, etc.
C    INFVAR : NWOPH
C    INFLIN : NCONST, NVARPT, NWOPPT
C

      IF (IPOLTP .GT. 0) LPOL = .TRUE.

      WRITE(LUPRI,*) 'LPOL =' ,LPOL

      XSAVE = DIPORG(1) 
      YSAVE = DIPORG(2) 
      ZSAVE = DIPORG(3) 
C
C     Core allocation
C
      KUCMO    = 1
      KFPEMO   = KUCMO    + NORBT*NBAST
      KFPEAC   = KFPEMO   + NNORBT
C     -------------------------------------------
      KINVMAT  = KFPEAC   + NNASHX
      KEFIELD  = KINVMAT  + 3*NNZAL*(3*NNZAL+1)/2
      KINDMOM  = KEFIELD  + 3*NNZAL*NCSIM
C     -------------------------------------------
      KFXC     = KINDMOM  + 3*NNZAL*NCSIM  
      KFXCAC   = KFXC     + NCSIM*NNORBT
      KTFXCAC  = KFXCAC   + NCSIM*NNASHX
C     -------------------------------------------
      KWRK1    = KTFXCAC  + NCSIM
      LWRK1    = LFREE    - KWRK1
C
C     1. KUCMO   : MO coefficients 
C     2. KFPEMO  : Fg(PE) operator mo basis
C     3. KFPEAC  : active part of Fg(PE)
C      ------------------------------------------
C     4. KINVMAT : [alpha]^(-1)
C     5. KEFIELD : Electric field on MM (polarizable) 
C        sites due to F^(1) field 
C        F(tilde) = < 0 | Fel(1) | B > (for B each state)
C     6. KINDMOM : induced moments (from NNZAL
C        polarizable sites)
C     ------------------------------------------ 
C     7. KFXC    : Fxc(PE) operator
C     8. KFXCAC  : active part of Fxc(PE) operator 
C     9. KTFXCAC : Vector of expectation values 
C        < 0 | Fxc | B >  (for B each state)

      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KFPEMO),NNORBT)
      CALL DZERO(WRK(KFPEAC),NNASHX)
      CALL DZERO(WRK(KINVMAT),3*NNZAL*(3*NNZAL+1)/2)
      CALL DZERO(WRK(KEFIELD),3*NNZAL*NCSIM)
      CALL DZERO(WRK(KINDMOM),3*NNZAL*NCSIM)
      CALL DZERO(WRK(KFXC),NCSIM*NNORBT)
      CALL DZERO(WRK(KFXCAC),NCSIM*NNASHX)
      CALL DZERO(WRK(KTFXCAC),NCSIM)     
C
      IF (LWRK1 .LT. 0) CALL ERRWRK('pe_lnc',-KWRK1,LWRK1)
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
C     ... Define F(PE) xc operator using function SOLELM (DTV is  < 0 | Fxc | B > type matrix )
C
C     1) Construct B(r) response (relay) matrix from file or solve iteratively

      N = 3*NNZAL

      IF (LPOL .AND. MMMAT) THEN
      LUQMMM = -1
        IF (LUQMMM .LT. 0) THEN
        CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
        ENDIF
      REWIND(LUQMMM)

        IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
          CALL READT(LUQMMM,N*(N+1)/2,WRK(KINVMAT))
        ELSE
          CALL QUIT('Problem reading the matrix from the QMMMIM file.')
        ENDIF

        CALL GPCLOSE(LUQMMM,'KEEP')

      ENDIF

C     2.a) ...Loop over MM sites

      IF (.NOT. LPOL) GOTO 755

      KMATAO   = KWRK1
      KMATMO   = KMATAO + 3*NNBASX
      KMATAC   = KMATMO + 3*NNORBT
      KWRK2    = KMATAC + 3*NNASHX
      LWRK2    = LFREE  - KWRK2

C     1. KMATAO  : QM dipole one-elctron integrals (Fel^(1) in ao basis)
C     2. KMATMO  : QM dipole one-elctron integrals (mo basis)
C     3. KMATAC  : Active part of F^(1)

      IF (LWRK2 .LT. 0) CALL ERRWRK('pe_lnc',-KWRK2,LWRK2)

      LRI = 0 ! counter for index in transformed electric field vector

      DO I = 1,MMCENT

        KPATOM = 0
        NCOM  = 3       ! called NOSIM occationally but denoted NCOM here 
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.

        CALL DZERO(WRK(KMATAO),3*NNBASX)
        CALL DZERO(WRK(KMATMO),3*NNORBT)
        CALL DZERO(WRK(KMATAC),3*NNASHX)

        DIPORG(1) = MMCORD(1,I)
        DIPORG(2) = MMCORD(2,I)
        DIPORG(3) = MMCORD(3,I)

C       2.b) Fel^(1) operator in AO basis

C       ...Get Fel^(1) integral: 1) x-coord. 2) y-coord. 3) z-coord.  

         RUNQM3 = .TRUE.

         CALL GET1IN(WRK(KMATAO),'NEFIELD',NCOM,WRK(KWRK2),   
     &               LWRK2,LABINT,INTREP,INTADR,I,TOFILE,        
     &               KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPQMMM)

         RUNQM3 = .FALSE.

C        2.c) Fel^(1) operator in MO basis

         CALL UTHU(WRK(KMATAO),WRK(KMATMO),WRK(KUCMO),
     &             WRK(KWRK2),NBAST,NORBT)

         CALL UTHU(WRK(KMATAO + 1*NNBASX),WRK(KMATMO + 1*NNORBT),
     &             WRK(KUCMO),WRK(KWRK2),NBAST,NORBT)

         CALL UTHU(WRK(KMATAO + 2*NNBASX),WRK(KMATMO + 2*NNORBT),
     &             WRK(KUCMO),WRK(KWRK2),NBAST,NORBT)

C        2.d) Active part of Fel(1) operator in MO basis

         IF (NASHT .GT. 0) THEN
                  
             CALL GETAC2(WRK(KMATMO),
     &                   WRK(KMATAC))

             CALL GETAC2(WRK(KMATMO + 1*NNORBT),
     &                   WRK(KMATAC + 1*NNASHX))

             CALL GETAC2(WRK(KMATMO + 2*NNORBT),
     &                   WRK(KMATAC + 2*NNASHX))
         ENDIF

C        2.e ) Make F(tilde) =  2 < 0 | F(el) | B > is the 2 included in SOLELM ???
C        for one |u > in | B > = sum(u)| u >b(u)        

         LCI = 0

         DO ICSIM = 1,NCSIM

            TXPE1 = SOLELM(DTV(1,ICSIM),WRK(KMATAC),
     &                     WRK(KMATMO),TXPEAC1)

            TXPE2 = SOLELM(DTV(1,ICSIM),WRK(KMATAC + 1*NNASHX),
     &                     WRK(KMATMO + 1*NNORBT),TXPEAC2)

            TXPE3 = SOLELM(DTV(1,ICSIM),WRK(KMATAC + 2*NNASHX),
     &                     WRK(KMATMO + 2*NNORBT),TXPEAC3)

C         ...To store the F(tilde) in dynamical memory, get
C         KEFIELD x, y, z first time loop is run. Next time  
C         Set LRI + 3 to get next MM center.         

C         x-value                                   ! This is for storage of the vector
          WRK(KEFIELD + LRI + 0 + LCI) = TXPEAC1    ! containing the expectation value of F(tilde)
C         y-value                                   ! LCI is a counter for each state in 
          WRK(KEFIELD + LRI + 1 + LCI) = TXPEAC2    ! < 0 | F(el) | B > = sum(u) < 0 | F(el) | u >
C         z-value
          WRK(KEFIELD + LRI + 2 + LCI) = TXPEAC3
C         start from x MM center of next root
          LCI = LCI + 3*NNZAL


        END DO ! NCSIM

        LRI = LRI + 3

      END DO ! MMCENT

C     ... Dot the KEFIELD matrix with B(r) matrix to get mu(tilde) for each | B > vector

      DO ICSIM = 1, NCSIM

      NDIM=3*NNZAL
      IF (MMMAT) THEN
         CALL DSPMV('L',NDIM,D1,WRK(KINVMAT),
     &              WRK(KEFIELD + (ICSIM-1)*3*NNZAL),1,D0,
     &              WRK(KINDMOM + (ICSIM-1)*3*NNZAL),1)
      ELSE IF (MMITER) THEN
         IOPT = 2 ! Do not read from file any previuos induced moments
         CALL F2QMMM(WRK(KEFIELD + 3*(ICSIM-1)*NNZAL),NDIM,
     &               WRK(KINDMOM + 3*(ICSIM-1)*NNZAL),
     &               WRK(KWRK2),LWRK2,IOPT,IPQMMM)
      ENDIF

      END DO ! ICSIM

C     3) Make F(el) and daxpy first x, then y and then z
C        for each CI B vector

       LRI = 0
         
        DO I = 1, MMCENT

        DIPORG(1) = MMCORD(1,I)
        DIPORG(2) = MMCORD(2,I)
        DIPORG(3) = MMCORD(3,I)

        CALL DZERO(WRK(KMATAO),3*NNBASX)
        CALL DZERO(WRK(KMATMO),3*NNORBT)
        CALL DZERO(WRK(KMATAC),3*NNASHX)   

C       5.a) Fel(1) operator in AO basis

        RUNQM3 = .TRUE.

        CALL GET1IN(WRK(KMATAO),'NEFIELD',NCOM,WRK(KWRK2),
     &               LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &               KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPQMMM)

        RUNQM3 = .FALSE.

C       5.b) Dipole one-electron integrals (Fel(1) operator in MO basis)

        CALL UTHU(WRK(KMATAO),WRK(KMATMO),WRK(KUCMO),
     &             WRK(KWRK2),NBAST,NORBT)

        CALL UTHU(WRK(KMATAO + 1*NNBASX),WRK(KMATMO + 1*NNORBT),
     &            WRK(KUCMO),WRK(KWRK2),NBAST,NORBT)

        CALL UTHU(WRK(KMATAO + 2*NNBASX),WRK(KMATMO + 2*NNORBT),
     &            WRK(KUCMO),WRK(KWRK2),NBAST,NORBT)

        
          LCI = 0 ! alternative set LCI = 3*NNZAL*(ICSIM-1)

          DO ICSIM = 1,NCSIM

            FACx = -WRK(KINDMOM + LRI + 0 + LCI)        

            CALL DAXPY(NNORBT,FACx,WRK(KMATMO),1,            
     &             WRK(KFXC + (ICSIM-1)*NNORBT),1)


            FACy = -WRK(KINDMOM + LRI + 1 + LCI)


            CALL DAXPY(NNORBT,FACy,WRK(KMATMO + 1*NNORBT),1, 
     &                 WRK(KFXC + (ICSIM-1)*NNORBT),1)

            FACz = -WRK(KINDMOM + LRI + 2 + LCI)       

            CALL DAXPY(NNORBT,FACz,WRK(KMATMO + 2*NNORBT),1, 
     &                 WRK(KFXC + (ICSIM-1)*NNORBT),1)

            IF (NASHT .GT. 0) THEN
              CALL GETAC2(WRK(KFXC + (ICSIM-1)*NNORBT),
     &                    WRK(KFXCAC + (ICSIM-1)*NNASHX))
            END IF

            LCI = LCI + 3*NNZAL

        END DO

      LRI = LRI + 3

      END DO

      DO ICSIM = 1,NCSIM

      TFXC = SOLELM(DV,WRK(KFXCAC + (ICSIM-1)*NNASHX),
     &              WRK(KFXC + (ICSIM-1)*NNORBT),TFXCAC)

      WRK(KTFXCAC-1+ICSIM) = TFXCAC
      
      END DO

  755 CONTINUE ! IF LPOL

      CALL PEFCMO(WRK(KUCMO),WRK(KFPEMO),DV,WRK(KWRK1),LWRK1,IPQMMM)

      IF (NASHT .GT. 0) THEN                             
          CALL GETAC2(WRK(KFPEMO),WRK(KFPEAC))
      END IF

      TFPEMO = SOLELM(DV,WRK(KFPEAC),WRK(KFPEMO),TFPEAC) 
C
C ****** edh: Write statements for debugging ******************
C *************************************************************

      IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(A,F17.8)')
     *   ' --- FPE expectation value MO :',TFPEMO
         WRITE (LUPRI,'(A,F17.8)')
     *   ' --- active part of FPE    :',TFPEAC
      
         WRITE (LUPRI,'(/A)') ' F(PE)_mo matrix in pe_lnc:'
         CALL OUTPAK(WRK(KFPEMO),  NORBT,1,LUPRI)
         IF (NASHT .GT. 0) THEN
         WRITE (LUPRI,'(/A)') ' F(PE)_ac matrix in pe_lnc:'
         CALL OUTPAK(WRK(KFPEAC),NASHT,1,LUPRI)
         END IF
      END IF

C *************************************************************
C *************************************************************

C    ...CSF part of sigma vectors    
 
      CALL SOLSC(NCSIM,0,BCVEC,CREF,SVEC,WRK(KFXCAC),WRK(KFPEAC), ! KRYCAC = KFPEAC (i.e. KRYC = KFPEMO)
     *           WRK(KTFXCAC),TFPEAC,INDXCI,WRK(KWRK1),LWRK1)      ! KRXCAC = KFXCAC (i.e. KRXC = KFXC  )  

      IF (NWOPPT .GT. 0) THEN
         MWOPH  = NWOPH
         NWOPH  = NWOPPT
C        ... tell SOLGO only to use the NWOPPT first JWOP entries
         JSVECO = 1 + NCONST
         JFXC   = KFXC
         DO ICSIM = 1,NCSIM
            IF (LPOL) CALL SOLGO(D2,DV,WRK(JFXC),SVEC(JSVECO,ICSIM))
            IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
               WRITE(LUPRI,*)
     *         ' orbital part of LINEAR TRANSFORMED CONF VEC No',ICSIM
               WRITE(LUPRI,*)' Fxc(PE) CONTRIBUTION'
               CALL OUTPUT(SVEC(JSVECO,ICSIM),1,NWOPPT,1,1,
     *                                        NWOPPT,1,1,LUPRI)
            END IF
            CALL SOLGO(D0,DTV(1,ICSIM),WRK(KFPEMO),SVEC(JSVECO,ICSIM))
            IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN 
               WRITE(LUPRI,*)
     *         ' orbital part of LINEAR TRANSFORMED CONF VEC No',ICSIM
               WRITE(LUPRI,*)' + Fg(PE) CONTRIBUTION'
               CALL OUTPUT(SVEC(JSVECO,ICSIM),1,NWOPPT,1,1,
     *                                        NWOPPT,1,1,LUPRI)
            END IF
            JFXC   = JFXC   + NNORBT 
         END DO
         NWOPH  = MWOPH
      END IF


C     ...Restore the dipole origin.

      DIPORG(1) = XSAVE
      DIPORG(1) = YSAVE
      DIPORG(1) = ZSAVE

      CALL QEXIT('pe_lnc')
      RETURN
      END
C     end of pelnc.

      SUBROUTINE pe_lno(NOSIM,BOVECS,CREF,CMO,INDXCI,    
     *                 DV,SVEC,NSVEC,WRK,LFREE)

      use polarizable_embedding, only: pe_master
C
C  Erik Donovan Hedegaard jan. 2012 
C  after original code by Hans Jorgen Aa. Jensen
C
C  Purpose:  Calculate MCSCF Hessian contribution from a
C            surrounding PE medium to an orbital trial vector.
C
C  NSVEC     may be NVAR or NWOPT, dependent on LINTRN
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "infinp.h"
#include "orgcom.h"
#include "inforb.h"
#include "infvar.h"
#include "inflin.h"
#include "inftap.h"
#include "qmmm.h"
#include "qm3.h"
#include "gnrinf.h"
C
C  Used from common blocks:
C    INFORB : NNASHX, NNORBX, NNBASX, etc.
C    INFVAR : JWOP
C    INFLIN : NWOPPT, NVARPT, NCONST, NCONRF
C
      DIMENSION BOVECS(NWOPPT,*), CREF(*), CMO(*)
      DIMENSION INDXCI(*),        DV(*)
      DIMENSION SVEC(NSVEC,*),    WRK(*)

C      DIMENSION DXCAO(N2BASX,*),DXVAO(N2BASX,*), UBO(N2ORBX,*) ! check which one - this or above !


      LOGICAL FULHES, TOFILE, TRIMAT, EXP1VL, LOCDEB, FNDLAB, LPOL

      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
C
      DOUBLE PRECISION D0, D2, D1, DP5
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0, D2 = 2.0D0, DP5 = 0.5D0 )

      PRINT*, 'IN pe_lno'

      LPOL = .FALSE. ! edh : Something needs to be done about LPOL !!

      CALL QENTER('pe_lno')

C     Determine if full Hessian or only orbital Hessian
C
      FULHES = (NSVEC .EQ. NVARPT)

      IF (IPOLTP .GT. 0) LPOL = .TRUE.

      IF (FULHES) THEN
         JSOVEC = 1 + NCONST
      ELSE
         JSOVEC = 1
      END IF
C
C *************************************************************
C *************************************************************

      IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(//A)') ' >>> TEST OUTPUT FROM pe_lno >>>'
      END IF
      IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
         IF (FULHES) THEN
            WRITE (LUPRI,'(/A)') ' >>> pe_lno - svec(ci,1) on entry'
            DO 30 I = 1,NCONST
               IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *              ' conf #',I,SVEC(I,1)
 30         CONTINUE
         END IF
         WRITE (LUPRI,'(/A)') ' >>> pe_lno - svec(orb) on entry'
         CALL OUTPUT(SVEC(JSOVEC,1),1,NWOPPT,1,NOSIM,
     *        NSVEC,NOSIM,1,LUPRI)
      END IF

C *************************************************************
C *************************************************************

C    ...Save the dipole origin

C      XSAVE = DIPORG(1)
C      YSAVE = DIPORG(2)
C      ZSAVE = DIPORG(3)
C
C     Core allocation
C
      KUCMO   = 1
      KUBO    = KUCMO   + NORBT*NBAST
C     ------------------------------------------
C      KINVMAT = KUBO    + NOSIM*N2ORBX
C      KINDMOM = KINVMAT + 3*NNZAL*(3*NNZAL+1)/2
C      KEFIEX =  KINDMOM + 3*NOSIM*NNZAL
C     ------------------------------------------
C      KFXO    = KEFIEX  + 3*NOSIM*NNZAL
C     ------------------------------------------
      KFPEMO  = KFXO    + NNORBT*NOSIM  
      KFPESQ  = KFPEMO  + NNORBX        ! still needed
      KFPXSQ  = KFPESQ  + N2ORBX
      KFPX    = KFPXSQ  + N2ORBX
      KFPXAC  = KFPX    + NOSIM*NNORBX
C     -----------------------------------------
      KFXYOA  = KFPXAC  + NOSIM*NNASHX
C     ----------------------------------------- 1 to pe_lno new module code!!
      KDXCAO  = KFXYOA  + NOSIM
      KDXVAO  = KDXCAO  + N2BASX
      KDXAO   = KDXVAO  + N2BASX
      KDFAO   = KDXAO   + NNBASX*NOSIM
      KDFMO   = KDFAO   + NNBASX*NOSIM
      KWRK1   = KDFMO   + NNORBX*NOSIM
      LWRK1   = LFREE   - KWRK1  


      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KUBO),NOSIM*N2ORBX)
C     ----------------------------------------- 1 to pe_lno new module code!!
C      CALL DZERO(WRK(KINVMAT), 3*NNZAL*(3*NNZAL+1)/2)
C      CALL DZERO(WRK(KINDMOM), 3*NOSIM*NNZAL)
C      CALL DZERO(WRK(KEFIEX), 3*NOSIM*NNZAL)
C      CALL DZERO(WRK(KFXO), NNORBT*NOSIM)
      CALL DZERO(WRK(KFPEMO), NNORBX)
      CALL DZERO(WRK(KFPESQ), N2ORBX)
      CALL DZERO(WRK(KFPXSQ), N2ORBX)
      CALL DZERO(WRK(KFPX),NOSIM*NNORBX)
      CALL DZERO(WRK(KFPXAC),NOSIM*NNASHX)
      CALL DZERO(WRK(KFXYOA),NOSIM)
C     ----------------------------------------- 1 to pe_lno new module code!!
      CALL DZERO(WRK(KDXCAO),N2BASX)
      CALL DZERO(WRK(KDXVAO),N2BASX)
      CALL DZERO(WRK(KDXAO),NNBASX*NOSIM)
      CALL DZERO(WRK(KDFAO),NNBASX*NOSIM)
      CALL DZERO(WRK(KDFMO),NNORBX)
      CALL DZERO(WRK(KUDFF),N2ORBX*NOSIM) 

      IF (LWRK1 .LT. 0) CALL ERRWRK('pe_lno',-KWRK1,LWRK1)
C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
C     Calculate unpacked orbital trial vectors in UBO 
C
      IF (NOSIM.GT.0) THEN
         DO IOSIM = 1,NOSIM
            JUBO = KUBO + (IOSIM-1)*N2ORBX
            CALL UPKWOP(NWOPPT,JWOP,BOVECS(1,IOSIM),WRK(JUBO))
            IF ( (IPQMMM .GE. 15) .OR. (LOCDEB) ) THEN
               WRITE (LUPRI,*) IOSIM,NOSIM
               CALL OUTPUT(WRK(JUBO),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                     LUPRI)
            END IF
         END DO
      END IF

!      IF (.NOT. LPOL) GOTO 755

C     1) Read B(r) response (Relay) matrix from file.

C      IF ( (LPOL) .AND. (MMMAT) ) THEN
C        N = 3*NNZAL
C        LUQMMM = -1
C        CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
C     &               'UNFORMATTED',IDUMMY,.FALSE.)
C        REWIND(LUQMMM)
C
C        IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
C          CALL READT(LUQMMM,N*(N+1)/2,WRK(KINVMAT))
C        ELSE
C          CALL QUIT('Problem reading the matrix from the QMMMIM file.')
C        ENDIF

C        CALL GPCLOSE(LUQMMM,'KEEP')

C      ENDIF

C      KPATOM = 0
C      NCOM  = 3       ! edh: sometimes called NOSIM but denoted NCOM here 
C      TOFILE = .FALSE.       
C      TRIMAT = .TRUE. 
C      EXP1VL = .FALSE.      

C     2) Construct Fxo(PE) = B(r) * < 0 | f^(1)el | 0 >  ; f(1)el is one-index transformed F(1)el

C     .. memory allocation for field matrix and one-electron transform.      

C      KMTAO   = KWRK1
C      KMTMO   = KMTAO  + 3*NNBASX
C      KMTSQ   = KMTMO  + 3*NNORBT  
C      KMTXSQ  = KMTSQ  + 3*N2ORBX
C      KMTX    = KMTXSQ + 3*N2ORBX 
C      KMTXAC  = KMTX   + 3*NOSIM*NNORBX 
C      KWRK2   = KMTXAC + 3*NOSIM*NNASHX
C      LWRK2   = LFREE   - KWRK2

C     1. KMTAO ("KMAT") : QM dipole one-elctron integrals (F^(1) in ao basis)
C     2. KMTMO    -     : QM dipole one-elctron integrals (mo basis)
C     3. KMTSQ    -     : Unpacked F^(1) (needed for one-index transform)
C     4. KMTXSQ   -     : One-index transformed F^(1)el ( =  f^(1)el )
C     5. KMTX     -     : f^(1)el triangular packed
C     6. KMTXAC   -     : Active part of f^(1)el

C      CALL DZERO(WRK(KMTAO),3*NNBASX)
C      CALL DZERO(WRK(KMTMO),3*NNORBT)
C      CALL DZERO(WRK(KMTSQ),3*N2ORBX)
C      CALL DZERO(WRK(KMTXSQ),3*N2ORBX)
C      CALL DZERO(WRK(KMTX),3*NOSIM*NNORBX)
C      CALL DZERO(WRK(KMTXAC),3*NOSIM*NNASHX)

C      IF (LWRK2 .LT. 0) CALL ERRWRK('pe_lno',-KWRK2,LWRK2)

C ********************* UNDER CONSTRUCTION ****************************

      ! workground for pe_lno !

      DO IOSIM    = 1,NOSIM
          JUBO    = KUBO     + (IOSIM - 1)*N2ORBX
          JDXAO   = KDXAO    + (IOSIM - 1)*NNBASX               ! NNBASX OR N2BASX
          CALL TR1DEN(CMO,WRK(JUBO),DV,WRK(KDXCAO),WRK(KDXVAO), ! KDXCAO : Trial vectors a.o. (inactive matrix)
     &                WRK(KWRK1),LWRK1)                         ! KDXVAO : Trial vectors a.o. (active matrix) 

          IF (NASHT .GT. 0) THEN
             CALL DAXPY(N2BASX,1.0D0,WRK(KDXVAO), 1, WRK(KDXCAO), 1)
          END IF

          CALL DGEFSP(NBAST,WRK(KDXCAO),WRK(JDXAO))             ! Symmetry pack (fold) trial vectors in a.o.

C      WRITE(LUPRI,*) 'KDXAO input to pe_master - IOSIM:', IOSIM
C      CALL OUTPAK(WRK(JDXAO),NBAST,1,LUPRI)

      END DO

      call pe_master(runtype='response',
     &              denmats=wrk(KDXAO:KDXAO+nosim*nnbasx-1), ! input for pe module "response"
     &              fckmats=wrk(KDFAO:KDFAO+nosim*nnbasx-1), ! output for pe module "response 
     &              nmats=nosim, dalwrk=wrk(KWRK1:KWRK1+LWRK1-1))
            
C    Takes fock matrices ()     
C    Returns which in MO basis was B < 0 | f^(1) | 0 > F^(1) in KDFAO 

C      CALL DGEFSP(NBAST,WRK(KDENT),WRK(KDENTF))
C
C      call pe_master('fock', denmats=wrk(kdentf:kdentf+nnbasx-1),
C     &               fckmats=wrk(kfpe:kfpe+nnbasx-1), nmats=1,
C     &               Epe=wrk(kwrk1:kwrk1),
C     &               dalwrk=wrk(kwrk1+1:kwrk1+1+lwrk1))
C      eqmmm = wrk(kwrk1)

C     Returns Fg type of contribution... (note: Fy type is obtained by PEFCMO!)

      DO IOSIM  = 1,NOSIM
         JDFMO  = KDFMO  + (IOSIM-1)*N2BASX                          ! trial vectors (final transformed to MO and folded)
         JDFAO  = KDFAO  + (IOSIM-1)*NNBASX                          ! contributions to trial vectors in a.o.
         CALL UTHU(WRK(JDFAO),WRK(KDFMO),CMO,WRK(KWRK1),NBAST,NORBT)
      END DO      

C**********************************************************************

C      LRI = 0 ! counter for index in one-index transformed electric field vector

C      DO I = 1,MMCENT

C        DIPORG(1) = MMCORD(1,I)
C        DIPORG(2) = MMCORD(2,I)
C        DIPORG(3) = MMCORD(3,I)

C       2.a Dipole one-electron integrals (Fel(1) operator in AO basis) 
C       ...Get F^(1)el integral: 1) x-coord. 2) y-coord. 3) z-coord.  

C       RUNQM3 = .TRUE.

C        CALL GET1IN(WRK(KMTAO),'NEFIELD',NCOM,WRK(KWRK2),
C     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
C     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPQMMM)

C         RUNQM3 = .FALSE.

C        2.b Dipole one-electron integrals (F^(1)el operator in MO basis)

C         CALL UTHU(WRK(KMTAO),WRK(KMTMO),WRK(KUCMO),
C     &             WRK(KWRK2),NBAST,NORBT)

C         CALL UTHU(WRK(KMTAO + 1*NNBASX),WRK(KMTMO + 1*NNORBT),
C     &             WRK(KUCMO),WRK(KWRK2),NBAST,NORBT)

C         CALL UTHU(WRK(KMTAO + 2*NNBASX),WRK(KMTMO + 2*NNORBT),
C     &             WRK(KUCMO),WRK(KWRK2),NBAST,NORBT)

C        2.c F^(1)el from packed (triangular) to unpacked (square) 

C         CALL DSPTSI(NORBT,WRK(KMTMO),WRK(KMTSQ))

C         CALL DSPTSI(NORBT,WRK(KMTMO + 1*NNORBT),
C     &               WRK(KMTSQ + 1*N2ORBX))

C         CALL DSPTSI(NORBT,WRK(KMTMO + 2*NNORBT),
C     &               WRK(KMTSQ + 2*N2ORBX))

 
C        DO IOSIM = 1, NOSIM
           
C            JUBO   = KUBO   + (IOSIM - 1) * N2ORBX         ! Unpacked orbital trial vectors
C            JMTX   = KMTX   + 3 * (IOSIM - 1) * NNORBX     ! F^(1) for each orb. trial vector          
C            JMTXAC = KMTXAC + 3 * (IOSIM - 1) * NNASHX     ! - active part         

C            CALL DZERO(WRK(KMTXSQ),3*N2ORBX) 

C            CALL TR1UH1(WRK(JUBO),WRK(KMTSQ),              !     **** x component **** 
C     &                  WRK(KMTXSQ),1)                     ! one index transform F^(1)el to f^(1)el

C            CALL DGETSP(NORBT,WRK(KMTXSQ),                 !     pack (triangular) f^(1)el
C     &                  WRK(JMTX))                         

C            CALL TR1UH1(WRK(JUBO),WRK(KMTSQ + 1*N2ORBX),   !     **** y component ****
C     &                  WRK(KMTXSQ + 1*N2ORBX),1)          ! one index transform F^(1) to f^(1) 

C            CALL DGETSP(NORBT,WRK(KMTXSQ + 1*N2ORBX),      !     pack (triangular) f^(1) 
C     &                  WRK(JMTX + 1*NNORBX))              


C            CALL TR1UH1(WRK(JUBO),WRK(KMTSQ + 2*N2ORBX),   !     **** z component **** 
C     &                  WRK(KMTXSQ + 2*N2ORBX),1)          ! one index transform F^(1) to f^(1) z

C            CALL DGETSP(NORBT,WRK(KMTXSQ + 2*N2ORBX),      !     pack (triangular) f^(1) 
C     &                  WRK(JMTX + 2*NNORBX))              


C           IF (NASHT .GT. 0) THEN
C             CALL GETAC2(WRK(JMTX),WRK(JMTXAC))
C             CALL GETAC2(WRK(JMTX + 1*NNORBX),WRK(JMTXAC + 1*NNASHX))
C             CALL GETAC2(WRK(JMTX + 2*NNORBX),WRK(JMTXAC + 2*NNASHX))
C           END IF   

C          ... Calculate < 0 | f^(1) | 0 >

C           TFX1 = SOLELM(DV,WRK(JMTXAC),                      
C     &                           WRK(JMTX),TFXAC1)
C           TFX2 = SOLELM(DV,WRK(JMTXAC + 1*NNASHX),
C     &                           WRK(JMTX + 1*NNORBX),TFXAC2)
C           TFX3 = SOLELM(DV,WRK(JMTXAC + 2*NNASHX),
C     &                           WRK(JMTX + 2*NNORBX),TFXAC3)


C        **** x-component ****                          
C         WRK(KEFIEX + 3*NNZAL*(IOSIM - 1) + LRI + 0) = TFX1 
C        **** y-component ****                         
C         WRK(KEFIEX + 3*NNZAL*(IOSIM - 1) + LRI + 1) = TFX2 
C        **** z-component ****
C         WRK(KEFIEX + 3*NNZAL*(IOSIM - 1) + LRI + 2) = TFX3
C        ... start from x of the next MM center

C           END DO ! NOSIM

C         LRI = LRI + 3

C      END DO ! MMCENT

C     ... and calculate the one-index transformed 
C     induced moment:  u = B * < 0 | f^(1) | 0 > 


C      DO IOSIM = 1, NOSIM

C        IF (IPOLTP .GT. 0) THEN
C          IF (MMMAT) THEN
C           CALL DSPMV('L',3*NNZAL,D1,WRK(KINVMAT),             ! edh: note KINVMAT is a lower triangular matrix
C     &                WRK(KEFIEX + 3*(IOSIM - 1)*NNZAL),1,D0,  
C     &                WRK(KINDMOM + 3*(IOSIM-1)*NNZAL),1)      
C          ELSE IF (MMITER) THEN
C           IOPT = 2 ! Do not read from file any previuos induced moments
C           CALL F2QMMM(WRK(KEFIEX  + 3*(IOSIM - 1)*NNZAL),NNZAL,
C     &                 WRK(KINDMOM + 3*(IOSIM-1)*NNZAL),  
C     &                 WRK(KWRK2),LWRK2,IOPT,IPQMMM)
C          ENDIF
C        END IF
C      END DO

C     3) Make F^(1)el and daxpy to get one-index transformed u; first x, then y and then z
C        
C      LRI = 0

C        DO I = 1, MMCENT

C        DIPORG(1) = MMCORD(1,I)
C        DIPORG(2) = MMCORD(2,I)
C        DIPORG(3) = MMCORD(3,I)

C        CALL DZERO(WRK(KMTAO),3*NNBASX)
C        CALL DZERO(WRK(KMTMO),3*NNORBT)

C       3.a) F^(1)el operator in AO basis

C        RUNQM3 = .TRUE.

C        CALL GET1IN(WRK(KMTAO),'NEFIELD',NCOM,WRK(KWRK2),
C     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
C     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPQMMM)

C        RUNQM3 = .FALSE.

C       3.b) Dipole one-electron integrals (F^(1)el operator in MO basis)

C        CALL UTHU(WRK(KMTAO),WRK(KMTMO),WRK(KUCMO),
C     &            WRK(KWRK2),NBAST,NORBT)

C        CALL UTHU(WRK(KMTAO + 1*NNBASX),WRK(KMTMO + 1*NNORBT),
C     &            WRK(KUCMO),WRK(KWRK2),NBAST,NORBT)

C        CALL UTHU(WRK(KMTAO + 2*NNBASX),WRK(KMTMO + 2*NNORBT),
C     &            WRK(KUCMO),WRK(KWRK2),NBAST,NORBT)

C        DO IOSIM = 1, NOSIM

C          FACx = -WRK(KINDMOM + LRI + 0 + 3*NNZAL*(IOSIM-1))        

C          CALL DAXPY(NNORBT,FACx,WRK(KMTMO),1,            
C     &               WRK(KFXO+(IOSIM-1)*NNORBT),1)

C          FACy = -WRK(KINDMOM + LRI + 1 + 3*NNZAL*(IOSIM-1))

C          CALL DAXPY(NNORBT,FACy,WRK(KMTMO + 1*NNORBT),1, 
C     &               WRK(KFXO+(IOSIM-1)*NNORBT),1)

C          FACz = -WRK(KINDMOM + LRI + 2 + 3*NNZAL*(IOSIM-1))

C          CALL DAXPY(NNORBT,FACz,WRK(KMTMO+2*NNORBT),1, 
C     &               WRK(KFXO+(IOSIM-1)*NNORBT),1)

C        END DO

C        LRI = LRI + 3

C      END DO

C     Construct Fyo(PE) (corresponds to Fg(PE) one-index transformed)

C  755 CONTINUE ! .NOT. LPOL

       CALL pe_fcmo(WRK(KUCMO),WRK(KFPEMO),DV,WRK(KWRK1),LWRK1,IPQMMM)

C      CALL PEFCMO(WRK(KUCMO),WRK(KFPEMO),DV,WRK(KWRK1),LWRK1,IPQMMM)

      CALL DSPTSI(NORBT,WRK(KFPEMO),WRK(KFPESQ))

      DO IOSIM = 1, NOSIM

         JUBO   = KUBO   + (IOSIM - 1) * N2ORBX         ! Unpacked orbital trial vectors
         JFPX   = KFPX   + (IOSIM - 1) * NNORBX         ! KFPX = Fyo(PE) 
         JFPXAC = KFPXAC + (IOSIM - 1) * NNASHX         ! - active part 

C         JTEST = KFXO    + (IOSIM - 1) * NNORBX            

C         JUDFF = KUDFF   + (IOSIM - 1) * NNORBX
          JDFMO = KDFMO   + (IOSIM - 1) * NNORBX

         WRITE(LUPRI,*) 'KDFMO output from pe_master - IOSIM:', IOSIM
         CALL OUTPAK(WRK(KDFMO),NBAST,1,LUPRI)

         CALL DZERO(WRK(KFPXSQ),N2ORBX)
         CALL DZERO(WRK(JFPX),NNORBX)
         CALL DZERO(WRK(JFPXAC),NNASHX)

         CALL TR1UH1(WRK(JUBO),WRK(KFPESQ),WRK(KFPXSQ),1)

         CALL DGETSP(NORBT,WRK(KFPXSQ),WRK(JFPX))  ! edh: JFPX = Fyo oeprator

C         IF (LPOL) CALL DAXPY(NNORBX,D1,WRK(JTEST),1,WRK(JFPX),1) ! Adds Fxo to Fyo when there are polarization contr.  

C         CALL DAXPY(NNORBX,D1,WRK(JTEST),1,WRK(JFPX),1)             ! Adds Fxo to Fyo always - if no pol then Fxo = 0. 

        CALL DAXPY(NNORBX,D1,WRK(KDFMO),1,WRK(JFPX),1)             ! Adds Fxo to Fyo always - if no pol then Fxo = 0

           IF (NASHT .GT. 0) THEN                                 ! active part of f(PE)g operator (equivalent to Tg = Ryo in RFSCF)
             CALL GETAC2(WRK(JFPX),WRK(JFPXAC))
           END IF

        FXYO = SOLELM(DV,WRK(JFPXAC),WRK(JFPX),FXYOA)

        WRK(KFXYOA + (IOSIM-1)) = FXYOA

      END DO

C     ... CSF part of sigma vectors

      IF (LSYMRF .EQ. LSYMST) THEN
         NCOLIM = 1
      ELSE
         NCOLIM = 0
      END IF
      IF (FULHES .AND. NCONST .GT. NCOLIM) THEN

        CALL SOLSC(0,NOSIM,DUMMY,CREF,SVEC,WRK(KFPXAC),DUMMY,
     &             100.D0*(WRK(KFXYOA)),DUMMY,INDXCI,WRK(KWRK1),LWRK1)
      END IF

C     ... orbital part of sigma vectors

      MWOPH  = NWOPH
      NWOPH  = NWOPPT
C    ... tell SOLGO only to use the NWOPPT first JWOP entries
      DO IOSIM = 1,NOSIM
        JFPX   =  KFPX  + (IOSIM-1)*NNORBX
        CALL SOLGO(D2,DV,WRK(JFPX),SVEC(JSOVEC,IOSIM))
      END DO
      NWOPH  = MWOPH

C     ...Restore the dipole origin.

      DIPORG(1) = XSAVE
      DIPORG(1) = YSAVE
      DIPORG(1) = ZSAVE

      CALL QEXIT('pe_lno')
      RETURN
C     ... end of pelno.
      END
