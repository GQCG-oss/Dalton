

C*****************************************************************************
      SUBROUTINE MP2SRDFTDEN(NFROZ,DONE,CMO,ORBEN,WRK,LFRSAV)
C*****************************************************************************
C
C written by Emmanuel Fromager - June 2007
C
C PURPOSE:
C calculate the second-order MP-SRDFT correction D^{(2),SRDFT} to the density matrix
C from the "long-range" second order contribution D^{(2),LR}, which is
C the standard MP2 correction to the density matrix where HF-SRDFT
C orbitals and epsilons are used and where 1/r12 is replaced by the
C long-range interaction. 
C D^{(2),SRDFT} = sum^{+infty}_{n=0} F^{n}(D^{(2),LR}).
C
C The F operator is defined as follows :
C
C F : DMO (density matrix contribution) ---> F(DMO) (new density matrix contribution) 
C
C i,j -> hole orbitals
C r,s -> virtual orbitals
C p,q -> hole or virtual orbitals
C F(DMO)_{ri}= 2/(epsilon_i-epsilon_r) * 
C int dx dx' K_{Hxc}^{sr}(x,x')*Omega^{MO}_{ri}(x)*sum_{pq}(DMO)_{pq} * Omega^{MO}_{pq}(x'). 
C F(DMO)_{ir}=F(DMO)_{ri}
C F(DMO)_{ij}=F(DMO)_{rs}=0
C K_{Hxc}^{sr}(x,x') is the second functional derivative of the SR Hxc functional 
C calculated at the HF-SRDFT density.
C 
C input  : DONE
C output : DONE = sum^{+infty}_{n=0} F^{n}(DONE)
C
C*****************************************************************************
#include "implicit.h"
C
      DIMENSION DONE(NORBT,NORBT),NFROZ(8),ORBEN(*)
      DIMENSION CMO(*),WRK(*)
      double precision getnormf 
C
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 )
C
C     DOUBLE PRECISION D2LRNORMF, FOPNORMF0, FOPNORMF1, FOPNORMF  
      DOUBLE PRECISION FACTOR, FOPNORMF   
C
#include "inforb.h"
#include "dummy.h"
#include "infpri.h"
#include "priunit.h"
#include "maxorb.h"
#include "infinp.h"
#include "scbrhf.h"
#include "inftap.h"
C
C
C     FACTOR = 0.0 ! equivalent to FNEUMANN 
C
C     FOPNORMF = 1.0 
C     FOPNORMF = 1.26 
C
      write(lupri,'(//A//A/A/A//A/)')
     & '-----------------------------------------',
     & 'Calculation of the Frobenius norm ',
     & '      of the F operator           ',
     & '(used in the auxiliary Neumann serie) ...',
     & '-----------------------------------------'
C
      FOPNORMF = GETNORMF(NFROZ,DONE,CMO,ORBEN,WRK,LFRSAV)
C     FOPNORMF = 1.0 
C     FOPNORMF = 1.26 
C
      FACTOR = D1/(D1+FOPNORMF)
C
C     FACTOR = 0.0 ! equivalent to FNEUMANN 
C     FACTOR = 0.46 
C     FACTOR = 0.44 
C     FACTOR = 0.84 
C     FACTOR = D1
C
      write(lupri,'(/A//A/A/A//A,F20.10//A/)')
     & '-----------------------------------------',
     & '   Calculation of the Neumann serie      ',
     & '   based on the auxiliary G operator     ',
     & 'and the transformed lr-MP2 density matrix',
     & 'info: factor used = ', FACTOR,
     & '-----------------------------------------'
C
      CALL GNEUMANN(FACTOR,NFROZ,DONE,CMO,ORBEN,WRK,LFRSAV)
C
      RETURN
      END


C*****************************************************************************
      FUNCTION GETNORMF(NFROZ,DONE,CMO,ORBEN,WRK,LFRSAV)
C*****************************************************************************
C
C written by Emmanuel Fromager - August 2007
C
C This function calculates the Frobenius norm of the F operator defined
C as lim_{i -> infty} ||F^{i+1}(D^{(2)lr})||_F/||F^{i}(D^{(2)lr})||_F
C*****************************************************************************
#include "implicit.h"
C
      DIMENSION DONE(NORBT,NORBT),NFROZ(8),ORBEN(*)
      DIMENSION CMO(*),WRK(*)
C
      double precision matnormf
C
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 )
C
      DOUBLE PRECISION D2LRNORMF, FOPNORMF0, FOPNORMF1, FOPNORMF  
      DOUBLE PRECISION FOPNORMF2  
C
#include "inforb.h"
#include "dummy.h"
#include "infpri.h"
#include "priunit.h"
#include "maxorb.h"
#include "infinp.h"
#include "scbrhf.h"
#include "inftap.h"
C
      KFRSAV = 1
      KFREE = KFRSAV
      LFREE = LFRSAV
C
C     Calculate the length NRFMO of the array RFMO that contains the lower
C     (virtuals-holes) rectangles of the density matrix (one for each
C     symmetry isym)
C
      NRFMO = 0
      DO ISYM = 1,NSYM
         NORBI  = NORB(ISYM)
         NRHFI  = NRHF(ISYM)
         NFROZI = NFROZ(ISYM)
         NRFMO = NRFMO + (NORBI-NRHFI)*(NRHFI-NFROZI)  
      END DO 
C
C     Max. number of iterations  
C
      MAXITER = 15 
C     MAXITER = 10 
C     MAXITER = 3 
C
C     Initialize the number of iterations
C
      ITERNORM = 0
C
C     Allocation of memory   
C
      CALL MEMGET('REAL',KDHFAO,2*N2BASX,WRK,KFREE,LFREE)
      KDAO = KDHFAO + N2BASX
      CALL MEMGET('REAL',KFAO,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRFMO,NRFMO,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KEXCMAT,N2BASX,WRK,KFREE,LFREE)
C          
C
C 
C     write (lupri,*) 'foperator in getfacopt',
C    & done(1,1),kdao,kfao,kdhfao,kexcmat,krfmo,kfree,lfree,nrfmo
C 
C     calculate the Frobenius norm of D^{(2),LR}
C
      D2LRNORMF = MATNORMF(DONE,NORBT,WRK(KFREE))
      FOPNORMF0 = D2LRNORMF 
C Manu debug b         
C     write(lupri,*) 'before applying F: D2LRNORMF =', D2LRNORMF
C     write(lupri,*) 'before applying F: D2LRNORMF mul=', 
C    &                D2LRNORMF * 1D10
C Manu debug e         
C
      CALL FOPERATOR(D1,DONE,.FALSE.,WRK(KDAO),WRK(KFAO),
     &               WRK(KDHFAO),WRK(KEXCMAT),NFROZ,WRK(KRFMO),CMO,
     &               ORBEN,WRK(KFREE),LFREE,NRFMO)
C 
C     
C     Calculate FOPNORMF2 = || F(D^{(2),LR}) ||_F / || D^{(2),LR} ||_F     
C     
      FOPNORMF1 = D2 * DDOT(NRFMO,WRK(KRFMO),1,WRK(KRFMO),1)    
      FOPNORMF1 = sqrt(FOPNORMF1)
      FOPNORMF2  = FOPNORMF1 / FOPNORMF0 
      FOPNORMF0 = FOPNORMF1
C     
C     Threshold for the convergence of the norm 
C     
      THRNORMF = 1.D-3 
C     
 3000 CALL FOPERATOR(D1,WRK(KRFMO),.TRUE.,WRK(KDAO),WRK(KFAO),
     &                   WRK(KDHFAO),WRK(KEXCMAT),NFROZ,
     &                   WRK(KRFMO),CMO,ORBEN,WRK(KFREE),LFREE,NRFMO)
C     
C     Calculate FOPNORMF = || F^{i+1}(D^{(2),LR}) ||_F / || F^{i}(D^{(2),LR}) ||_F, i>= 1     
C     
      FOPNORMF1 = D2 * DDOT(NRFMO,WRK(KRFMO),1,WRK(KRFMO),1)    
      FOPNORMF1 = sqrt(FOPNORMF1)
      FOPNORMF  = FOPNORMF1 / FOPNORMF0 
      FOPNORMF0 = FOPNORMF1
C Manu debug b      
C     write(lupri,*) 'print FOPNORMF before dividing by D2LRNORMF,
C    &                ', FOPNORMF
C     write(lupri,*) 'print D2LRNORMF before dividing by D2LRNORMF,
C    &                ', D2LRNORMF
C Manu debug e      
C     
C     If not converged ...    
C     
      IF (ABS(FOPNORMF-FOPNORMF2) .GE. THRNORMF) THEN  
         WRITE(LUPRI,'(/A)') '------------------------------------' 
         WRITE(LUPRI,*) 'Iteration ',ITERNORM 
         WRITE(LUPRI,'(A12,F17.12)') ' FOPNORMF = ', FOPNORMF   
         WRITE(LUPRI,*) 'The norm of the F operator' 
         WRITE(LUPRI,*) 'is not converged yet !' 
         WRITE(LUPRI,'(A/)') '------------------------------------' 
C
C
         ITERNORM = ITERNORM + 1 
C
C        If the max. number of iterations is reached ...
C
         IF (ITERNORM .EQ. MAXITER) THEN
            WRITE(LUPRI,*) '**************    info    **************' 
            WRITE(LUPRI,*) 'The norm of F is not converged !!!'
            WRITE(LUPRI,*) 'Increase MAXITER (currently', MAXITER, ')'
            WRITE(LUPRI,*) 'Anyway, let us use that one ...'
            WRITE(LUPRI,*) '****************************************' 
C
            GO TO 4000 
C
         ELSE    
C
             FOPNORMF2=FOPNORMF ! save the norm to check convergence at the next iteration
C            
C            Ready to apply the F operator once again ! 
C 
C
             GO TO 3000          
C
         END IF
C
C     else, the norm is converged ...
C
      ELSE
         WRITE(LUPRI,*) '***********************************' 
         WRITE(LUPRI,*) '***********************************' 
         WRITE(LUPRI,*) 'Norm of the F operator converged'
         WRITE(LUPRI,*) 'after', ITERNORM, ' iterations !' 
         WRITE(LUPRI,*) 'final norm: ', FOPNORMF
         WRITE(LUPRI,*) 'Information about the convergence:'
         WRITE(LUPRI,'(A18,F14.12)') ' diff. FOPNORMF = ',
     &                                FOPNORMF-FOPNORMF2   
         WRITE(LUPRI,*) '***********************************' 
         WRITE(LUPRI,*) '***********************************' 
C
C
         GO TO 4000
C
      END IF 
C
C     return the result and free the memory 
C
 4000 GETNORMF = FOPNORMF  
      CALL MEMREL('Releasing mem. used for DHFAO,DAO,FAO,RFMO,EXCMAT',
     &             WRK,KFRSAV,KDHFAO,KFREE,LFREE)  
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE GNEUMANN(FACTOR,NFROZ,DONE,CMO,ORBEN,WRK,LFRSAV)
C*****************************************************************************
C
C written by Emmanuel Fromager - August 2007
C
C This subroutine calculates the Neumann serie  
C DONE = sum^{+infty}_{n=0} G^{n}(DONE1) where G = a * F^2 +(1-a) *F
C and DONE1 = a * F(DONE) + DONE
C The F operator is defined in the FOPERATOR subroutine
C The a coeff (FACTOR) is given in input
C*****************************************************************************
#include "implicit.h"
C
      DIMENSION DONE(N2ORBX),NFROZ(8),ORBEN(*)
      DIMENSION CMO(*),WRK(*)
C
      double precision matnormf
      DOUBLE PRECISION FACTOR
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 , D2 = 2.0D0 )
C
      DOUBLE PRECISION TD2LRNORMF, GOPNORMF0, GOPNORMF1, GOPNORMF  
      DOUBLE PRECISION AF2DLRNORMF, D2LRNORMF1, DONENORMF, DONENORMF1 
      DOUBLE PRECISION RFMOTOTNORMF, SUPRFMOTOT, SCOVERLR
C
#include "inforb.h"
#include "dummy.h"
#include "infpri.h"
#include "priunit.h"
#include "maxorb.h"
#include "infinp.h"
#include "scbrhf.h"
#include "inftap.h"
C
      KFRSAV = 1
      KFREE = KFRSAV
      LFREE = LFRSAV
C
C     Calculate the length NRFMO of the array RFMO that contains the lower
C     (virtuals-holes) rectangles of the density matrix (one for each
C     symmetry isym)
C
      NRFMO = 0
      DO ISYM = 1,NSYM
         NORBI  = NORB(ISYM)
         NRHFI  = NRHF(ISYM)
         NFROZI = NFROZ(ISYM)
         NRFMO = NRFMO + (NORBI-NRHFI)*(NRHFI-NFROZI)  
      END DO 
C
C     Max. number of iterations  
C
      MAXITER = 15 
C     MAXITER = 10 
C     MAXITER = 3 
C
C     Initialize ... 
C
      ITERDEN = 0
      SUPRFMO = D0
      SUPRFMOTOT = D0
      AF2DLRNORMF = D0
      D2LRNORMF1 = D0
      DONENORMF = D0
      DONENORMF1 = D0
      RFMOTOTNORMF = D0
      SCOVERLR = D0 
      SCOVERLRINFTY = D0 
C
C     Allocation of memory   
C
      CALL MEMGET('REAL',KDHFAO,2*N2BASX,WRK,KFREE,LFREE)
      KDAO = KDHFAO + N2BASX
      CALL MEMGET('REAL',KFAO,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRFMO,NRFMO,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KEXCMAT,N2BASX,WRK,KFREE,LFREE)
C
C     self-consistency contributions will be separated from D2lr for
C     analysis purposes ...
      CALL MEMGET('REAL',KRFMOTOT,NRFMO,WRK,KFREE,LFREE)
      CALL DZERO(WRK(KRFMOTOT),NRFMO) ! initialization to zero 
C
C          
C     For analysis purposes: get the Frobenius norm of D^{(2),lr} ...
C
      D2LRNORMF1 = MATNORMF(DONE,NORBT,WRK(KFREE))
C      
C     ... and largest matrix element of D^{(2),lr} 
C
      IMAXDONE = IDAMAX(NORBT*NORBT,DONE,1)
      SUPDONE  = ABS(DONE(IMAXDONE)) 
C 
C     write (lupri,*) 'foperator in gneumann',
C    & done(1,1),kdao,kfao,kdhfao,kexcmat,krfmo,kfree,lfree,nrfmo
C 
C     If FACTOR is not zero calculate FACTOR * F(DONE) and ...
C 
C 
      IF ( FACTOR.NE.D0 ) THEN
C 
C     print *, 'we are in GNEUMANN and FACTOR differs from 0'

         CALL FOPERATOR(FACTOR,DONE,.FALSE.,WRK(KDAO),WRK(KFAO),
     &                  WRK(KDHFAO),WRK(KEXCMAT),NFROZ,WRK(KRFMO),CMO,
     &                  ORBEN,WRK(KFREE),LFREE,NRFMO)
C 
C     For analysis purposes: get max ( FACTOR * F(D^{(2),lr})) in
C     absolute value and the Frobenius norm of FACTOR * F(D^{(2),lr})
C 
         IMAXRFMO = IDAMAX(NRFMO,WRK(KRFMO),1)
         SUPRFMO  = WRK(KRFMO-1+IMAXRFMO) 
         IMAXAFD2LR = IMAXRFMO
         AF2DLRNORMF = D2 * DDOT(NRFMO,WRK(KRFMO),1,WRK(KRFMO),1)    
         AF2DLRNORMF = SQRT(AF2DLRNORMF)
C 
C     ... add to DONE 
C 
         CALL ADDVH(DONE,WRK(KRFMO),NFROZ) 
C             ADDVH(DONE,RFMO,NFROZ) 
C     
C        For analysis purposes: sc contribution a * F (D2lr) added to
C        WRK(KRFMOTOT)    
C
         CALL DAXPY(NRFMO,D1,WRK(KRFMO),1,WRK(KRFMOTOT),1)
C        write(lupri,*) 'print WRK(KRFMOTOT) after 1st addvh'
C        CALL OUTPUT(WRK(KRFMOTOT),1,NRFMO,1,1,NRFMO,1,1,LUPRI)
C
C
      ELSE 
C     
C        print *, 'we are in GNEUMANN and FACTOR equals 0'
C     
      END IF
C
C
C     For analysis purposes : calculate the Frobenius norm of 
C     TD^{(2),LR} = FACTOR * F(D^{(2),LR}) + D^{(2),LR} 
C
      TD2LRNORMF = MATNORMF(DONE,NORBT,WRK(KFREE))
      GOPNORMF0 =  TD2LRNORMF 
C     D2LRNORMF = D2LRNORMF * 1.D8 ! (to avoid precision problems) 
C Manu debug b         
C     write(lupri,*) 'before applying F: D2LRNORMF =', D2LRNORMF
C     write(lupri,*) 'before applying F: D2LRNORMF mul=', 
C    &                D2LRNORMF * 1.D10
C Manu debug e         
C
      CALL GOPERATOR(FACTOR,DONE,.FALSE.,WRK(KDAO),WRK(KFAO),
     &               WRK(KDHFAO),WRK(KEXCMAT),NFROZ,WRK(KRFMO),CMO,
     &               ORBEN,WRK(KFREE),LFREE,NRFMO)
C 
C     Check convergence 
C     
      THRDEN = 1.D-8 ! convergence threshold for the density matrix 
C     THRDEN = 1.D-15 ! just for analysis ... 
C     THRDEN = 1.D-20 ! just for analysis ... 
C     THRDEN = 1.D-6 ! just for debugging ... 
C     
 3000 IMAXRFMO = IDAMAX(NRFMO,WRK(KRFMO),1)
      SUPRGMO  = WRK(KRFMO-1+IMAXRFMO) 
C     
C     Calculate GOPNORMF = || G^{i+1}(TD^{(2),LR}) ||_F / || G^{i}(TD^{(2),LR}) ||_F, i>= 0     
C     
      GOPNORMF1 = D2 * DDOT(NRFMO,WRK(KRFMO),1,WRK(KRFMO),1)    
      GOPNORMF1 = sqrt(GOPNORMF1)
C Manu debug b      
C     write(lupri,*) 'GOPNORMF0 =', GOPNORMF0
C     write(lupri,*) 'GOPNORMF1 =', GOPNORMF1
C Manu debug e      
      GOPNORMF  = GOPNORMF1 / GOPNORMF0 
      GOPNORMF0 = GOPNORMF1
C     
C     If not converged ...    
C     
      IF (ABS(SUPRGMO) .GE. THRDEN) THEN  
         WRITE(LUPRI,'(/A)')
     &    '---------------------------------------------' 
         WRITE(LUPRI,*) 'MP2-SRDFT density matrix: Iteration ',ITERDEN 
         WRITE(LUPRI,'(A11,F14.12)') ' SUPRGMO = ', SUPRGMO   
         WRITE(LUPRI,'(A12,F17.12)') ' GOPNORMF = ', GOPNORMF   
C        WRITE(LUPRI,'(A13,F14.12)') ' SUPaFD2LR = ', SUPRFMO   
C        WRITE(LUPRI,'(A14,F14.12)') ' TD2LRNORMF = ', TD2LRNORMF   
         WRITE(LUPRI,*) 'The density matrix is not converged yet !' 
         WRITE(LUPRI,'(A/)')
     &    '---------------------------------------------' 
C
C        "Add" RGMO to DONE      
C
C Manu debug b
C          write(lupri,*) 'print NORB'
C          write(lupri,*) 'NSYM = ', NSYM
C          DO ISYM = 1, NSYM
C          write(lupri,*) 'NORB(', ISYM, ') =', NORB(ISYM)
C          END DO
C          write(lupri,*) 'print NRHF'
C          write(lupri,*) 'NSYM = ', NSYM
C          DO ISYM = 1, NSYM
C          write(lupri,*) 'NRHF(', ISYM, ') =', NRHF(ISYM)
C          END DO
C          write(lupri,*) 'print NFROZ'
C          write(lupri,*) 'NSYM = ', NSYM
C          DO ISYM = 1, NSYM
C          write(lupri,*) 'NFROZ(', ISYM, ') =', NFROZ(ISYM)
C          END DO
C          write(lupri,*) 'print NBAS'
C          write(lupri,*) 'NSYM = ', NSYM
C          DO ISYM = 1, NSYM
C          write(lupri,*) 'NBAS(', ISYM, ') =', NBAS(ISYM)
C          END DO
C        write(lupri,*) 'print WRK(KRFMO) before addvh'
C        CALL OUTPUT(WRK(KRFMO),1,NRFMO,1,1,NRFMO,1,1,LUPRI)
C        WRITE(LUPRI,*) 'print DONE before addvh' 
C        CALL OUTPUT(DONE,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      call flshfo(lupri)
C Manu debug e
         CALL ADDVH(DONE,WRK(KRFMO),NFROZ) 
C             ADDVH(DONE,RFMO,NFROZ) 
C
C     We accumulate self-consistency contributions in WRK(KRFMOTOT) for
C     analysis purposes
C
         CALL DAXPY(NRFMO,D1,WRK(KRFMO),1,WRK(KRFMOTOT),1)
C        write(lupri,*) 'print WRK(KRFMOTOT) after 2nd addvh'
C        CALL OUTPUT(WRK(KRFMOTOT),1,NRFMO,1,1,NRFMO,1,1,LUPRI)
C     
C
C Manu debug b
C        WRITE(LUPRI,*) 'print DONE after addvh' 
C        CALL OUTPUT(DONE,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
C     call flshfo(lupri)
C Manu debug e
C
         ITERDEN = ITERDEN + 1 
C
C        If the max. number of iterations is reached ...
C
         IF (ITERDEN .EQ. MAXITER) THEN
            WRITE(LUPRI,*) '>>>>>>>>>>    MP2-SRDFT info   <<<<<<<<<<' 
            WRITE(LUPRI,*) 'The density matrix is not converged !!!'
            WRITE(LUPRI,*) 'Increase MAXITER (currently', MAXITER, ')'
            WRITE(LUPRI,*) 'Anyway, get the natural orbitals ...'
            WRITE(LUPRI,*) '>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<' 
C
            GO TO 4000 
C
         ELSE    
C
C          Ready to apply the G operator once again ! 
C 
           CALL GOPERATOR(FACTOR,WRK(KRFMO),.TRUE.,WRK(KDAO),WRK(KFAO),
     &                   WRK(KDHFAO),WRK(KEXCMAT),NFROZ,
     &                   WRK(KRFMO),CMO,ORBEN,WRK(KFREE),LFREE,NRFMO)
C 
C
           GO TO 3000          
C
         END IF
C
C     else, the density matrix is converged ...
C
      ELSE
C 
C        For analysis purposes: get the Frobenius norm of DONE ...
C   
         DONENORMF = MATNORMF(DONE,NORBT,WRK(KFREE))
C
C        ... and its ratio and difference with the norm of D2LR
C
         DONENORMF1 = ABS(DONENORMF-D2LRNORMF1)
         IF (D2LRNORMF1.GT.(1.D-7)) THEN
            DONENORMF = (DONENORMF/D2LRNORMF1) -1
            DONENORMF = ABS(DONENORMF)
         END IF
C
         WRITE(LUPRI,'(/A)') '***********************************' 
         WRITE(LUPRI,*) '***********************************' 
         WRITE(LUPRI,*) 'MP2-SRDFT density matrix converged'
         WRITE(LUPRI,*) 'after', ITERDEN, ' iterations !' 
         WRITE(LUPRI,*) 'Information about the convergence:'
         WRITE(LUPRI,'(A11,F14.12)') ' SUPRGMO = ', SUPRGMO   
         WRITE(LUPRI,'(A13,F14.12)') ' SUPaFD2LR = ', SUPRFMO   
C         
         IF (.NOT.(FACTOR.EQ.0)) THEN
            WRITE(LUPRI,*) 'ISUPaFD2LR = ', IMAXAFD2LR 
         END IF   
C         
         WRITE(LUPRI,'(A15,F14.12)') ' aF2DLRNORMF = ', AF2DLRNORMF 
         WRITE(LUPRI,'(A14,F14.12)') ' TD2LRNORMF = ', TD2LRNORMF   
         IF (D2LRNORMF1.GT.(1.D-7)) THEN
            WRITE(LUPRI,'(A28,F14.12)') ' | D2NORMF/D2LRNORMF -1 | = ', 
     &                                      DONENORMF  
         ENDIF
         WRITE(LUPRI,'(A25,F14.12)') ' | D2NORMF-D2LRNORMF | = ', 
     &                                   DONENORMF1  
         WRITE(LUPRI,'(A13,F14.12)') ' D2LRNORMF = ', D2LRNORMF1   
         WRITE(LUPRI,'(A11,F15.12)') ' SUPD2LR = ', SUPDONE   
         WRITE(LUPRI,*) '***********************************' 
         WRITE(LUPRI,*) '***********************************' 
C
C        Compute the Frobenius norm of the self-consistency contribution
C
         RFMOTOTNORMF = D2 * 
     &                  DDOT(NRFMO,WRK(KRFMOTOT),1,WRK(KRFMOTOT),1) 
         RFMOTOTNORMF = sqrt(RFMOTOTNORMF)
C
C        ... and its ratio with the Frobenius norm of the purely
C        long-range contribution  
C
         IF (D2LRNORMF1.GT.(1.D-7)) THEN
            SCOVERLR = (RFMOTOTNORMF/D2LRNORMF1) 
         END IF
C
C        Compute the maximum norm of the self-consistency contribution
C
         IMAXRFMOTOT=IDAMAX(NRFMO,WRK(KRFMOTOT),1)
         SUPRFMOTOT=ABS(WRK(KRFMOTOT-1+IMAXRFMOTOT))

C        ... and its ratio with the maximum norm of the purely
C        long-range contribution  
C
         IF (SUPDONE.GT.(1.D-7)) THEN
            SCOVERLRINFTY = (SUPRFMOTOT/SUPDONE) 
         END IF
C
         WRITE(LUPRI,'(/A)') '***********************************' 
         WRITE(LUPRI,*) '***********************************' 
         WRITE(LUPRI,*) 'Frobenius/infinite norm of the self-consistent'
         WRITE(LUPRI,*) 'contribution to the MP2-srDFT density matrix'
         WRITE(LUPRI,'(A26,F15.12)') ' ||D2-D2lr||_F= ', RFMOTOTNORMF 
         WRITE(LUPRI,'(A26,F15.12)') ' ||D2lr||_F= ', D2LRNORMF1 
         WRITE(LUPRI,'(A26,F15.12)') ' ||D2-D2lr||_infty= ', SUPRFMOTOT 
         WRITE(LUPRI,'(A26,F15.12)') ' ||D2lr||_infty= ', SUPDONE 
C
         IF (D2LRNORMF1.GT.(1.D-7)) THEN
            WRITE(LUPRI,'(A26,F7.4)') ' ||D2-D2lr||_F/||D2lr||_F='
     &                                  , SCOVERLR 
         END IF
         IF (SUPDONE.GT.(1.D-7)) THEN
            WRITE(LUPRI,'(A34,F7.4)') 
     &      '||D2-D2lr||_infty/||D2lr||_infty='
     &                                  , SCOVERLRINFTY 
         END IF
         WRITE(LUPRI,'(/A)') '***********************************' 
         WRITE(LUPRI,*) '***********************************' 
C
C
C        write(LUPRI,*) 'print the converged DONE '
C        call OUTPUT(DONE,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
C
         GO TO 4000
C
      END IF 
C
C     Free the memory and return
C
 4000 CALL MEMREL('Releasing mem. used for DHFAO,DAO,FAO,RFMO,EXCMAT
     &             and RFMOTOT',
     &             WRK,KFRSAV,KDHFAO,KFREE,LFREE)  
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE GOPERATOR(FACTOR,DMO,RECT,DAO,FAO,DHFAO,EXCMAT,NFROZ,
     &                     RGMO,CMO,ORBEN,WRK,LFRSAV,NRFMO)
C*****************************************************************************
C
C     written by Emmanuel Fromager - August 2007 (based on FOPERATOR)
C
C     This subroutine calculates the matrix RGMO = G(DMO)      
C     where G = FACTOR * F^2 + (1-FACTOR) * F 
C     The F operator is defined in the subroutine FOPERATOR
C
C*****************************************************************************
#include "implicit.h"
C 
      DIMENSION DMO(*),NFROZ(8),ORBEN(*)
      DIMENSION CMO(*),WRK(*),FAO(*),DAO(*),DHFAO(*),EXCMAT(*),RGMO(*)
      LOGICAL   RECT
      DOUBLE PRECISION FACTOR
C
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 )
C
#include "inforb.h"
#include "dummy.h"
#include "infpri.h"
#include "priunit.h"
#include "maxorb.h"
#include "infinp.h"
#include "scbrhf.h"
#include "inftap.h"
C
C 
C     If FACTOR is not zero then ...
C 
      IF (.NOT.(FACTOR.EQ.0)) THEN
C 
C     print *, 'we are in GOPERATOR and FACTOR differs from 0'
          KFRSAV = 1
          KFREE = KFRSAV
          LFREE = LFRSAV
C
          CALL MEMGET('REAL',KRFMO,NRFMO,WRK,KFREE,LFREE)
C
C         Calculate F(DMO)  
C
          CALL FOPERATOR(D1,DMO,RECT,DAO,FAO,DHFAO,EXCMAT,NFROZ,
     &                   WRK(KRFMO),CMO,ORBEN,WRK(KFREE),LFREE,NRFMO)
C
C         Calculate FACTOR * F^2(DMO) ...  
C
          CALL FOPERATOR(FACTOR,WRK(KRFMO),.TRUE.,DAO,FAO,DHFAO,
     &                   EXCMAT,NFROZ,RGMO,CMO,ORBEN,WRK(KFREE),
     &                   LFREE,NRFMO)
C
C         ... and add (1-FACTOR) * F(DMO)    
C     
          CALL DAXPY(NRFMO,D1-FACTOR,WRK(KRFMO),1,RGMO,1)
C
          CALL MEMREL('Releasing mem. used for RFMO in GOPERATOR',
     &               WRK,KFRSAV,KRFMO,KFREE,LFREE)  
C
      ELSE
C
C     print *, 'we are in GOPERATOR and FACTOR equals 0'
C
          CALL FOPERATOR(D1,DMO,RECT,DAO,FAO,DHFAO,EXCMAT,NFROZ,
     &                  RGMO,CMO,ORBEN,WRK,LFRSAV,NRFMO)

      END IF
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE FOPERATOR(FACTOR,DMO,RECT,DAO,FAO,DHFAO,EXCMAT,NFROZ,
     &                     RFMO,CMO,ORBEN,WRK,LFRSAV,NRFMO)
C*****************************************************************************
C
C written by Emmanuel Fromager - August 2007
C
C This subroutine calculates the matrix RFMO = FACTOR * F(DMO) IN THE MO BASIS. 
C If DMO is a packed rectangular matrix (we only keep the lower
C hole-particle block of the density matrix contribution DMO in that case) then RECT should be set to .TRUE.
C If DMO is an unpacked symmetric matrix RECT should be set to .FALSE.
C 
C The F operator is defined as follows :
C
C F : DMO (density matrix contribution) ---> F(DMO) (new density matrix contribution) 
C
C i,j -> hole orbitals
C r,s -> virtual orbitals
C p,q -> hole or virtual orbitals
C F(DMO)_{ri}= 2/(epsilon_i-epsilon_r) * 
C int dx dx' K_{Hxc}^{sr}(x,x')*Omega^{MO}_{ri}(x)*sum_{pq}(DMO)_{pq} * Omega^{MO}_{pq}(x'). 
C F(DMO)_{ir}=F(DMO)_{ri}
C F(DMO)_{ij}=F(DMO)_{rs}=0
C K_{Hxc}^{sr}(x,x') is the second functional derivative of the SR Hxc functional 
C calculated at the HF-SRDFT density.
C 
C input :  DMO
C output : RFMO = FACTOR * F(DMO)
C
C*****************************************************************************
#include "implicit.h"
C
      DIMENSION DMO(*),NFROZ(8),ORBEN(*)
      DIMENSION CMO(*),WRK(*),FAO(*),DAO(*),DHFAO(*),EXCMAT(*),RFMO(*)
      LOGICAL   RECT
      DOUBLE PRECISION FACTOR
C
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 )
C
C
#include "inforb.h"
#include "dummy.h"
#include "infpri.h"
#include "priunit.h"
#include "maxorb.h"
#include "infinp.h"
#include "scbrhf.h"
#include "inftap.h"
C
      KFRSAV = 1
      KFREE = KFRSAV
      LFREE = LFRSAV
C
C     Write the HF-SRDFT density matrix in the AO basis (contravariant transformation) 
C     We need it to calculate the second order functional derivative of
C     the SR xc functional at the HF-SRDFT density
C
      CALL MEMGET('REAL',KDHFMO,N2ORBX,WRK,KFREE,LFREE)
      CALL GETDHFAO(DHFAO,WRK(KDHFMO),CMO,WRK(KFREE)) 
C          GETDHFAO(DHFAO,DHFMO,CMO,WRK) 
C Manu debug b         
C     WRITE(LUPRI,*) 'print DHFMO'
C     CALL OUTPUT(WRK(KDHFMO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
C Manu debug e         
C          
      CALL MEMREL('Releasing mem. used for DHFMO',
     &             WRK,KFRSAV,KDHFMO,KFREE,LFREE)  

C     WRITE(LUPRI,*) 'print DHFAO'
C     CALL OUTPUT(DHFAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
C Manu debug b         
C                              T
C Let us calculate 2*CMO(1)*CMO(1) to check uhutb on H2 ...     
C
C     CALL MEMGET('REAL',KCCT,N2BASX,WRK,KFREE,LFREE)
C     CALL DGEMM('N','T',NBAST,NBAST,1,1.D0,
C    &                 CMO(1),NBAST,
C    &                 CMO(1),NBAST,0.D0,
C    &                 WRK(KCCT),NBAST)
C     CALL DSCAL(N2BASX,D2,WRK(KCCT),1)
C     WRITE(LUPRI,*) 'print WRK(KCCT)'
C     CALL OUTPUT(WRK(KCCT),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
C      
C Manu debug e         
C
C     For analysis purposes : calculate the Frobenius norm of D^{(2),LR} 
C
C     D2LRNORMF = MATNORMF(DONE,NORBT,WRK(KFREE))
C     FOPNORMF0 = D2LRNORMF 
C     D2LRNORMF = D2LRNORMF * 1D8 ! (to avoid precision problems) 
C Manu debug b         
C     write(lupri,*) 'before applying F: D2LRNORMF =', D2LRNORMF
C     write(lupri,*) 'before applying F: D2LRNORMF mul=', 
C    &                D2LRNORMF * 1D10
C Manu debug e         
C
C
C     Transform DMO from MO to AO basis (contravariant transformation) 
C
      IF (.NOT.RECT) THEN
C
C        WRITE(LUPRI,*) 'print DMO before applying F'
C        CALL OUTPUT(DMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
C
c        WRITE(LUPRI,*) 'print DONE (sym. 1) before uhutb'
c        CALL OUTPUT(DONE,1,NORB(1),1,NORB(1),NORBT,NORBT,1,LUPRI)
C
C
         CALL UHUTB(DAO,NBAST,DMO,NORBT,CMO,NSYM,NBAS,NORB,
     &              WRK(KFREE))
C
C
C          UHUTB(HAO,NAOT,HMO,NMOT,U,NSYM,NAO,NMO,WRK)
c     JCMOi = 1
c     do isym = 1,nsym  
c        nbasi = nbas(isym)
c        norbi = norb(isym)
c        write(lupri,*) 'print CMO after uhutb, sym =', isym
c        CALL OUTPUT(CMO(jcmoi),1,nbasi,1,norbi,nbasi,norbi,-1,LUPRI)
c        jcmoi = jcmoi + norbi*nbasi
c     end do 
      ELSE
C
c          write(lupri,*) 'print NORB before urutb'
c          write(lupri,*) 'NSYM = ', NSYM
c          DO ISYM = 1, NSYM
c          write(lupri,*) 'NORB(', ISYM, ') =', NORB(ISYM)
c          END DO
c          write(lupri,*) 'print NRHF before urutb'
c          write(lupri,*) 'NSYM = ', NSYM
c          DO ISYM = 1, NSYM
c          write(lupri,*) 'NRHF(', ISYM, ') =', NRHF(ISYM)
c          END DO
c          write(lupri,*) 'print NFROZ before urutb'
c          write(lupri,*) 'NSYM = ', NSYM
c          DO ISYM = 1, NSYM
c          write(lupri,*) 'NFROZ(', ISYM, ') =', NFROZ(ISYM)
c          END DO
c          write(lupri,*) 'print NBAS before urutb'
c          write(lupri,*) 'NSYM = ', NSYM
c          DO ISYM = 1, NSYM
c          write(lupri,*) 'NBAS(', ISYM, ') =', NBAS(ISYM)
c          END DO
C
C
C          WRITE(LUPRI,*) 'entering urutb'
C
           CALL URUTB(DAO,NBAST,DMO,CMO,
     &                NSYM,NBAS,NORB,NRHF,NFROZ,WRK(KFREE))
C             URUTB(RAO,NAOT,RMO,U,WRK,NSYM,NAO,NMO,NMOCC,NFROZ,WRK)
C
C          WRITE(LUPRI,*) 'print RFMO in the AO basis (DAO)' 
C          CALL OUTPUT(DAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF

C
C     Calculate the SR-Hartree contribution in the AO basis FAO_{pq}=
C     int dx dx' * w_{ee}^{sr}(x,x')*Omega^{AO}_{pq}(x)*sum_{mn}(DAO)_{mn} * Omega^{AO}_{mn}(x') 
C     = sum_{mn} (pq|mn)^{sr} * (DAO)_{mn}   
C
C
      CALL DZERO(FAO,N2BASX)
C
C     WRITE(LUPRI,*) 'print DAO in foperator'
C     CALL OUTPUT(DAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
C
      ISYMDM = 1  ! DAO of symmetry 1 
      IFCTYP = 11 ! Calculate Hartree term only 
      NDAO  = 1   ! only one density matrix
C
C     write(LUPRI,*) ' i am entering sirfck2'
C     ADDSRI = .FALSE. 
C     CALL SIRFCK(WRK(KFAO),WRK(DAO),NDAO,ISYMDM,IFCTYP,DIRECT,
C    &                  WRK(KFREE),LFREE)
      CALL SIRFCK2(LUSRINT,'AOSR2INT',FAO,DAO,
     &             NDAO,ISYMDM,IFCTYP,WRK(KFREE),LFREE)
C Manu debug b
C     WRITE(LUPRI,*) 'print FAO before calling srdft'
C     CALL OUTPUT(FAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      call flshfo(lupri)
C Manu debug e
C     write(LUPRI,*) ' i am leaving sirfck2'
C     ADDSRI = .TRUE. 
C          SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,DIRECT,
C    &                  WRK,LWORK)
C
C     Calculate the SR-xc contribution in the AO basis EXCMAT{pq}=
C     int dx dx' * K_{xc}^{sr}(x,x')*Omega^{AO}_{pq}(x)*sum_{mn}(DAO)_{mn} * Omega^{AO}_{mn}(x'). 
C     K_{xc}^{sr}(x,x') is the second functional derivative of the SR xc functional
C     calculated at the HF-SRDFT density.
C
      CALL DZERO(EXCMAT,N2BASX)
C     WRITE(LUPRI,*) 'print EXCMAT before call srdft'
C     CALL OUTPUT(WRK(KEXCMAT),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
C
!!!!!! are the AOs in the "right order" (by symmetry) ?  
C     write(LUPRI,*) ' i am entering srdft'
C
      CALL SRDFT(1,EXCMAT,DHFAO,VDUMMY,
     &           .FALSE.,.FALSE.,.TRUE.,.FALSE.,
     &           WRK(KFREE),LFREE,IPRFCK)
C          SRDFT(ND_SIM,EXCMAT,DMAT,EDFT(1:3),DOERG,DOGRD,DOATR,WORK,LWORK,IPRINT)
C
C     WRITE(LUPRI,*) 'print EXCMAT after calling srdft'
C     CALL OUTPUT(WRK(KEXCMAT),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
C     write(LUPRI,*) ' i am leaving srdft'
C
C     Add SR xc contributions to FAO  
C
      CALL DAXPY(N2BASX,D1,EXCMAT,1,FAO,1)
C     WRITE(LUPRI,*) 'print FAO+EXCAO'
C     CALL OUTPUT(FAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      call flshfo(lupri)
C
C     Transformation from AO to MO basis (FAO --> RFMO) :
C     we only keep the lower virtuals-holes rectangle RFMO(isym) 
C     for each symmetry isym 
C    
C
C Manu debug b
C     print *, 'ENTERING UTRUB'
C     print *, 'nbast,nsym,kfree',nbast,nsym,kfree
C Manu debug e
      CALL UTRUB(FAO,NBAST,RFMO,CMO,NSYM,
     &           NBAS,NORB,NRHF,NFROZ,WRK(KFREE))
C          UTRUB(RAO,NAOT,RMO,U,NSYM,NAO,NMO,NMOCC,NFROZ,WRK)
C
C     
C Manu debug b
C     print *, 'LEAVING UTRUB'
C Manu debug e
C
C     Scale RFMO with factor 2 * FACTOR 
C
      CALL DSCAL(NRFMO,D2*FACTOR,RFMO,1)
C
C     Multiply RFMO_{ri} by (epsilon_i-epsilon_r)
C
C     write(lupri,*) 'calling resolvh factor,nrfmo',factor,nrfmo
C     write(lupri,*) 'RFMO matrix before ORBEN in FOPERATOR'
C     call output(rfmo,1,1,1,nrfmo,1,nrfmo,1,lupri)
      CALL RESOLVH(RFMO,ORBEN,NFROZ) 
C          RESOLVH(RFMO,ORBEN,NFROZ) 
C     
      RETURN
      END


C*****************************************************************************
      SUBROUTINE UHUTB(HAO,NAOT,HMO,NMOT,U,NSYM,NAO,NMO,WRK) 
C*****************************************************************************
C
C     SUBROUTINE UHUTB(HAO,NAOT,HMO,NMOT,U,NAO,NMO,WRK)
C
C May-2007 Emmanuel Fromager (based on UTHUB)
C
C
C This subroutine transforms the NSYM symmetric, UNPACKED matrix HMO
C of dimension NMO(isym) to the symmetric, UNPACKED HAO of dimension
C NAO(isym).
C U(NAO(isym),NMO(isym)) is the transformation matrix.
C
C                                                       T
C  HAO( ij ) =  SUM (k,l = 1,...,NMO)  U (i,k) HMO(kl) U(l,j)
C
C This routine is optimized with respect to vector operations.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C#include "inforb.h"
      DIMENSION HAO(NAOT,NAOT),HMO(NMOT,NMOT),U(*),WRK(*)
      DIMENSION NAO(8),NMO(8)
C
      CALL DZERO(HAO,NAOT*NAOT) 
C
      IUOFF  = 1
      IHAOFF = 1
      IHMOFF = 1 ! commented if iorb is used
      DO 600 ISYM = 1,NSYM
         NMOI = NMO(ISYM)
      IF (NMOI .eq. 0) go to 600
         NAOI = NAO(ISYM)
C        IORBI = IORB(ISYM) ! just to try ...
C
            CALL DGEMM('N','T',NMOI,NAOI,NMOI,1.D0,
C    &                 HMO(IORBI+1,IORBI+1),NMOT,
     &                 HMO(IHMOFF,IHMOFF),NMOT, ! commented if iorb is used
     &                 U(IUOFF),NAOI,0.D0,
     &                 WRK(1),NMOI)
C Manu debug b
c           write(lupri,*) 'print U(iuoff) in uhutb, sym = ', isym
c           CALL OUTPUT(U(IUOFF),1,NAOI,1,nmoi,NAOI,nmoi,-1,LUPRI)
c           write(lupri,*) 'print WRK(1),isym in uhutb', WRK(1), isym
c           CALL OUTPUT(wrk,1,NmOI,1,naoi,NmOI,naoi,-1,LUPRI)
C           CALL OUTPUT(WRK,1,1,1,NMOI*NAOI,1,NMOI*NAOI,1,LUPRI)
C Manu debug e
            CALL DGEMM('N','N',NAOI,NAOI,NMOI,1.D0,
     &                 U(IUOFF),NAOI,
     &                 WRK(1),NMOI,0.D0,
     &                 HAO(IHAOFF,IHAOFF),NAOT) ! initial code
C    &                 HAO(IHAOFF,IHAOFF),NAOI) ! just for debugging
C
C Manu debug b
c           hone = 0
c           do i = 1,nmoi
c             do j =1,nmoi
c                hone = hone + U(IUOFF + (i-1)*naoi) * 
c    &                  HMO(IHMOFF -1 +i,IHMOFF -1+j) *
c    &                  U(IUOFF + (j-1)*naoi) 
c             end do
c           end do
c           write(lupri,*) 'HAO(isym =',isym,')(1,1) in uhutb =', hone
C
c           write(lupri,*) 'HAO(isym =',isym,') in uhutb'
c           CALL OUTPUT(HAO(IHAOFF,IHAOFF),1,NAOI,
c    &                  1,NAOI,NAOT,NAOT,1,LUPRI)
C Manu debug e
         IUOFF  = IUOFF + NAOI*NMOI
         IHAOFF = IHAOFF + NAOI
         IHMOFF = IHMOFF + NMOI ! commented if iorb is used
  600 CONTINUE
C Manu debug b
C     WRITE(LUPRI,*) 'print '
C     CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
C Manu debug e
C
      RETURN
      END
C  /* Deck urutb */
      SUBROUTINE URUTB(RAO,NAOT,RMO,U,NSYM,NAO,NMO,NMOCC,NFROZ,WRK)
C
C May-2007 Emmanuel Fromager (based on UHUTB)
C
C
C This subroutine transforms the RECTANGULAR, PACKED matrix RMO
C of dimension (NMO(isym) - NMOCC(isym)) * (NMOCC(isym)-NFROZ(isym))
C to the SYMMETRIC, UNPACKED RAO of dimension NAO(isym).
C U(NAO(isym),NMO(isym)) is the transformation matrix written as
C 
C             NFROZ(isym)   NMOCC(isym)-NFROZ(isym)      NMO(isym) - NMOCC(isym)    
C            <---------->   <---------------------->   <---------------------->       
C U(isym) = [ ---------- |       UH(isym)            |          UV(isym)       ]    
C
C
C                         T          T    T               
C  and RAO = UV * RMO * UH + UH * RMO * UV  
C
C This routine is optimized with respect to vector operations.
C
#include "implicit.h"
      DIMENSION RAO(NAOT,NAOT),RMO(*),U(*),WRK(*)
C     DIMENSION NAO(NSYM),NMO(NSYM)
C     DIMENSION NMOCC(NSYM),NFROZ(NSYM)
      DIMENSION NAO(8),NMO(8)
      DIMENSION NMOCC(8),NFROZ(8)
C
      CALL DZERO(RAO,NAOT*NAOT) 
C
      IUHOFF = 1 
      IUVOFF = 1
      IRAOFF = 1 
      IRMOFF = 1  
C      
C          print *, 'print NAO in urutb'
C          print *, 'NSYM = ', NSYM
C          DO ISYM = 1, NSYM
C          print *, 'NAO(', ISYM, ') =', NAO(ISYM)
C          END DO
C          write(lupri,*) 'print NAO in urutb'
C          write(lupri,*) 'NSYM = ', NSYM
C          DO ISYM = 1, NSYM
C          write(lupri,*) 'NAO(', ISYM, ') =', NAO(ISYM)
C          END DO
C          write(lupri,*) 'print NMO in urutb'
C          write(lupri,*) 'NSYM = ', NSYM
C          DO ISYM = 1, NSYM
C          write(lupri,*) 'NMO(', ISYM, ') =', NMO(ISYM)
C          END DO
C          write(lupri,*) 'print NMOCC in urutb'
C          write(lupri,*) 'NSYM = ', NSYM
C          DO ISYM = 1, NSYM
C          write(lupri,*) 'NMOCC(', ISYM, ') =', NMOCC(ISYM)
C          END DO
C          write(lupri,*) 'print NFROZ in urutb'
C          write(lupri,*) 'NSYM = ', NSYM
C          DO ISYM = 1, NSYM
C          write(lupri,*) 'NFROZ(', ISYM, ') =', NFROZ(ISYM)
C          END DO
C      
C      
C     NMOCCI = 0
C      
C      
      DO 600 ISYM = 1,NSYM
         NMOI   = NMO(ISYM)
         NAOI   = NAO(ISYM)
         NMOCCI = NMOCC(ISYM)
         NFROZI = NFROZ(ISYM)
C
         IF (NMOI .eq. 0) go to 600
         IF (NMOI.eq.NMOCCI .OR. (NMOCCI-NFROZI).eq.0) THEN ! no virtual or hole orbitals
            IUHOFF = IUHOFF + NAOI*NMOI 
            IUVOFF = IUVOFF + NAOI*NMOI 
            IRAOFF = IRAOFF + NAOI
            go to 600
         END IF   
C
         IUHOFF = IUHOFF + NAOI*NFROZI
         IUVOFF = IUVOFF + NAOI*NMOCCI
C
C                               T
C           WRK = RMO(isym) * UH(isym)   
C                          
C           print *, 'call the 1st dgemm' 
C           print *, 'NMO(', ISYM, ') = ', NMOI
C           print *, 'NMOCCI = ', NMOCCI
C           print *, 'NMOCC(', ISYM, ') = ', NMOCC(ISYM)
C           print *, 'NORB(', ISYM, ') = ', NORB(ISYM)
C           print *, 'NAO(', ISYM, ') = ', NAOI
C           print *, 'NFROZ(', ISYM, ') = ', NFROZI
C
            CALL DGEMM('N','T',NMOI-NMOCCI,NAOI,NMOCCI-NFROZI,1.D0,
     &                 RMO(IRMOFF),NMOI-NMOCCI,
     &                 U(IUHOFF),NAOI,0.D0,
     &                 WRK(1),NMOI-NMOCCI)
C     
C           RAO(isym) = UV(isym) * WRK           
C                          
C           print *, 'call the 2nd dgemm' 
C
            CALL DGEMM('N','N',NAOI,NAOI,NMOI-NMOCCI,1.D0,
     &                 U(IUVOFF),NAOI,
     &                 WRK(1),NMOI-NMOCCI,0.D0,
C    &                 RAO(IRAOFF,IRAOFF),NAOI) ! just for debugging
     &                 RAO(IRAOFF,IRAOFF),NAOT) ! initial code
C
C                                      T   
C           RAO(isym) = RAO(isym) + RAO(isym)
C
            DO I = IRAOFF,(IRAOFF+NAOI-1) 
               DO J = IRAOFF,I
                  RAO(I,J) = RAO(I,J) + RAO(J,I) 
                  RAO(J,I) = RAO(I,J)
               END DO
            END DO
C   
         IUHOFF = IUHOFF + NAOI*(NMOI-NFROZI) 
         IUVOFF = IUVOFF + NAOI*(NMOI-NMOCCI) 
         IRAOFF = IRAOFF + NAOI
         IRMOFF = IRMOFF + (NMOI-NMOCCI)*(NMOCCI-NFROZI)
  600 CONTINUE
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE UTRUB(RAO,NAOT,RMO,U,NSYM,NAO,NMO,NMOCC,NFROZ,WRK)
C*****************************************************************************
C
C May-2007 Emmanuel Fromager (based on URUTB)
C
C
C This subroutine transforms the SYMMETRIC, UNPACKED matrix RAO 
C of dimension NAO(isym) to the RECTANGULAR, PACKED matrix RMO
C of dimension (NMO(isym) - NMOCC(isym)) * (NMOCC(isym)-NFROZ(isym)).
C U(NAO(isym),NMO(isym)) is the transformation matrix written as
C 
C             NFROZ(isym)   NMOCC(isym)-NFROZ(isym)      NMO(isym) - NMOCC(isym)    
C            <---------->   <---------------------->   <---------------------->       
C U(isym) = [ ---------- |       UH(isym)            |          UV(isym)       ]    
C
C
C              T
C  and RMO = UV * RAO * UH 
C
C This routine is optimized with respect to vector operations.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION RAO(NAOT,NAOT),RMO(*),U(*),WRK(*)
C     DIMENSION NAO(NSYM),NMO(NSYM)
C     DIMENSION NMOCC(NSYM),NFROZ(NSYM)
      DIMENSION NAO(8),NMO(8)
      DIMENSION NMOCC(8),NFROZ(8)
C
C     Calculate the length NRMO of RMO  
C
      NRMO = 0
      DO ISYM = 1,NSYM
         NMOI   = NMO(ISYM)
         NMOCCI = NMOCC(ISYM)
         NFROZI = NFROZ(ISYM)
         NRMO   = NRMO + (NMOI-NMOCCI)*(NMOCCI-NFROZI)
      END DO 
C
      CALL DZERO(RMO,NRMO)
C
      IUHOFF = 1 
      IUVOFF = 1
      IRAOFF = 1 
      IRMOFF = 1  
C      
      DO 600 ISYM = 1,NSYM
         NMOI   = NMO(ISYM)
         NAOI   = NAO(ISYM)
         NMOCCI = NMOCC(ISYM)
         NFROZI = NFROZ(ISYM)
C        print *,'isym...',isym,nmoi,naoi,nmocci,nfrozi
C
         IF (NMOI .eq. 0) go to 600
         IF (NMOI.eq.NMOCCI .OR. (NMOCCI-NFROZI).eq.0) THEN ! no virtual or hole orbitals
            IUHOFF = IUHOFF + NAOI*NMOI 
            IUVOFF = IUVOFF + NAOI*NMOI 
            IRAOFF = IRAOFF + NAOI
            go to 600
         END IF   
C
         IUHOFF = IUHOFF + NAOI*NFROZI
         IUVOFF = IUVOFF + NAOI*NMOCCI
C
C                          
C           WRK = RAO(isym) * UH(isym)   
C
C Manu debug b
C           write(lupri,*) 'symmetry nr.', isym
C           write(lupri,*) 'print UH(',isym,') in utrub'
C           CALL OUTPUT(U(IUHOFF),1,NAOI,1,NMOCCI-NFROZI,
C    &                  NAOI,NMOCCI-NFROZI,-1,LUPRI)
C           write(lupri,*) 'print RAO(',isym,') in utrub'
C           CALL OUTPUT(RAO(IRAOFF,IRAOFF),1,NAOI,1,NAOI,
C    &                  NAOT,NAOT,1,LUPRI)
C Manu debug e
C
            CALL DGEMM('N','N',NAOI,NMOCCI-NFROZI,NAOI,1.D0,
     &                 RAO(IRAOFF,IRAOFF),NAOT,
     &                 U(IUHOFF),NAOI,0.D0,
     &                 WRK(1),NAOI)
C     
C           
C Manu debug b
C           write(lupri,*) 'print WRK(', isym,'),isym in utrub'
C           CALL OUTPUT(wrk,1,NAOI,1,NMOCCI-NFROZI,NAOI,
C    &                  NMOCCI-NFROZI,-1,LUPRI)
C Manu debug e
C                               T     
C           RMO(isym) = UV(isym) * WRK           
C
C Manu debug b
C           write(lupri,*) 'print UV(',isym,') in utrub'
C           CALL OUTPUT(U(IUVOFF),1,NAOI,1,NMOI-NMOCCI,
C    &                  NAOI,NMOI-NMOCCI,-1,LUPRI)
C Manu debug e
C
            CALL DGEMM('T','N',NMOI-NMOCCI,NMOCCI-NFROZI,NAOI,1.D0,
     &                 U(IUVOFF),NAOI,
     &                 WRK(1),NAOI,0.D0,
     &                 RMO(IRMOFF),NMOI-NMOCCI)
C
C Manu debug b
C           write(lupri,*) 'print RMO(',isym,') in utrub'
C           CALL OUTPUT(RMO(IRMOFF),1,NMOI-NMOCCI,1,NMOCCI-NFROZI,
C    &                  NMOI-NMOCCI,NMOCCI-NFROZI,-1,LUPRI)
C Manu debug e
C
         IUHOFF = IUHOFF + NAOI*(NMOI-NFROZI) 
         IUVOFF = IUVOFF + NAOI*(NMOI-NMOCCI) 
         IRAOFF = IRAOFF + NAOI
         IRMOFF = IRMOFF + (NMOI-NMOCCI)*(NMOCCI-NFROZI)
  600 CONTINUE
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE RESOLVH(RFMO,ORBEN,NFROZ) 
C*****************************************************************************
C
C June-2007 Emmanuel Fromager 
C
C
C This subroutine multiplies the elements RFMO_{ri} 
C of the RECTANGULAR PACKED matrix RFMO by 1/(epsilon_i - epsilon_r)
C r --> virtuals, i --> holes
C 
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION RFMO(*),ORBEN(*),NFROZ(8)   
C 
#include "inforb.h"
C
      IRFMO = 0 
C
      DO 1500 ISYM = 1,NSYM
         IORBI  = IORB(ISYM)
         NRHFI  = NRHF(ISYM)
         NORBI  = NORB(ISYM)
         NFROZI = NFROZ(ISYM)
         NVIRI  = NORBI - NRHFI 
         NHOLI  = NRHFI - NFROZI 
         if (nviri .eq. 0 .or. nholi.eq.0) goto 1500
C        write(lupri,*) 'RESOLV isym',isym
C        call output(RFMO(IRFMO+1),1,nviri,1,nholi,nviri,nholi,1,lupri)
         DO 1600 I = 1,NVIRI 
            DO 1700 J = 1,NHOLI
               DENOM = ORBEN(J+IORBI+NFROZI) - ORBEN(I+IORBI+NRHFI)
C              write(lupri,*) 'i,j,denom',i,j,denom
               RFMO(IRFMO+I+(J-1)*NVIRI) = 
     &         RFMO(IRFMO+I+(J-1)*NVIRI)/DENOM 
 1700       CONTINUE
 1600    CONTINUE
C        call output(RFMO(IRFMO+1),1,nviri,1,nholi,nviri,nholi,1,lupri)
         IRFMO = IRFMO + NVIRI*NHOLI
 1500 CONTINUE
C     
      RETURN
      END
C
C  /* Deck getdhfao */
      SUBROUTINE GETDHFAO(DHFAO,DHFMO,CMO,WRK) 
C
C June-2007 Emmanuel Fromager 
C
C This subroutine transforms the HF density matrix in the AO basis 
C
#include "implicit.h"
C
      DIMENSION DHFAO(NBAST,NBAST),DHFMO(NORBT,NORBT)  
      DIMENSION CMO(*),WRK(*)
      PARAMETER (D2 = 2.0D0)
C
#include "inforb.h"
#include "priunit.h"
C
      KFRSAV = 1
      KFREE = KFRSAV
C     LFREE = LFRSAV
C
      CALL DZERO(DHFAO,N2BASX)
      CALL DZERO(DHFMO,N2ORBX)
C
      DO 1500 ISYM = 1,NSYM
         IORBI = IORB(ISYM)
         NRHFI = NRHF(ISYM)
         NORBI = NORB(ISYM)
C
         DO 1600 I = (IORBI+1),(IORBI+NRHFI)
C
            DHFMO(I,I) = D2
C
 1600    CONTINUE
 1500 CONTINUE
C
C     write(lupri,*) 'print DHFMO in getdhfao', DHFMO
      CALL UHUTB(DHFAO,NBAST,DHFMO,NORBT,CMO,NSYM,NBAS,NORB,
     &           WRK(KFREE))
C          UHUTB(HAO,NAOT,HMO,NMOT,U,NSYM,NAO,NMO,WRK)
C     write(lupri,*) 'print DHFAO in getdhfao', DHFAO
C
      RETURN
      END
C
C  /* Deck addzeroth */
      SUBROUTINE ADDZEROTH(DONE,DOADD) 
C
C June-2007 Emmanuel Fromager 
C
C This subroutine adds (if doadd = true) or substracts (if doadd = false) 
C the zeroth order contribution to the density matrix 
C
#include "implicit.h"
C
      DIMENSION DONE(NORBT,NORBT)
      LOGICAL   DOADD
      PARAMETER (D2 = 2.0D0)
C
#include "inforb.h"
C
      DO 1500 ISYM = 1,NSYM
         IORBI = IORB(ISYM)
         NRHFI = NRHF(ISYM)
         NORBI = NORB(ISYM)
C
         DO 1600 I = (IORBI+1),(IORBI+NRHFI)
C
         IF (DOADD) THEN 
            DONE(I,I) = D2 + DONE(I,I)
         ELSE   
            DONE(I,I) = DONE(I,I) - D2
         END IF    
C
 1600    CONTINUE
 1500 CONTINUE
C
      RETURN
      END
C
C  /* Deck addvh */
      SUBROUTINE ADDVH(DONE,RFMO,NFROZ) 
C
C June-2007 Emmanuel Fromager 
C     
C This subroutine adds the RECTANGULAR PACKED matrix RFMO to the
C virtuals-holes and holes-virtuals rectangles of the
C SYMMETRIC UNPACKED matrix DONE 
C
#include "implicit.h"
      DIMENSION DONE(NORBT,NORBT),RFMO(*),NFROZ(8)   
C 
#include "inforb.h"
C 
         IRFMO = 0 
         DO 2500 ISYM = 1,NSYM
            IORBI  = IORB(ISYM)
            NRHFI  = NRHF(ISYM)
            NORBI  = NORB(ISYM)
            NFROZI = NFROZ(ISYM)
            NVIRI  = NORBI - NRHFI 
            NHOLI  = NRHFI - NFROZI 
            DO 2600 I = 1,NVIRI 
               DO 2700 J = 1,NHOLI
                  DONE(I+IORBI+NRHFI,J+IORBI+NFROZI) = 
     &            DONE(I+IORBI+NRHFI,J+IORBI+NFROZI) +
     &            RFMO(IRFMO+I+(J-1)*NVIRI) 
C
                  DONE(J+IORBI+NFROZI,I+IORBI+NRHFI) = 
     &            DONE(I+IORBI+NRHFI,J+IORBI+NFROZI) 
 2700          CONTINUE
 2600       CONTINUE
            IRFMO = IRFMO + NVIRI*NHOLI
 2500    CONTINUE
C
      RETURN
      END
C
C  /* Function matnormf */
      DOUBLE PRECISION FUNCTION MATNORMF(DONE,N,WRK)  
C
C August-2007 Emmanuel Fromager 
C
C This function calculates the Frobenius norm of the UNPACKED DONE
C matrix:
C                                          T    
C        ||DONE||_F = sqrt [ tr (DONE * DONE) ]
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION DONE(N,N), WRK(*) 
C
      XXX = 0.0D0
C      
C Manu debug b      
C     WRITE(LUPRI,*) 'in MATNORMF: print DONE before calculating DONE2 '
C     CALL OUTPUT(DONE,1,N,1,N,N,N,1,LUPRI)
C Manu debug e     
      CALL DGEMM('N','T',N,N,N,1.D0,DONE,N,DONE,N,0.D0,WRK,N)
C Manu debug b      
C     WRITE(LUPRI,*) 'in MATNORMF: print DONE2'
C     CALL OUTPUT(WRK,1,N,1,N,N,N,1,LUPRI)
C Manu debug e     
C
      XXX = TRACEMAT(WRK,N) 
C      
C Manu debug b      
C     write(lupri,*) 'in MATNORMF : trace of DONE2 =', XXX
C Manu debug e      
C      
      XXX = sqrt(XXX) 
C Manu debug b      
C     write(lupri,*) 'in MATNORMF : sqrt of trace of DONE2 =', XXX
C     write(lupri,*) 'in MATNORMF : sqrt of trace of DONE2 mul =', 
C    &                XXX * 1D10
C Manu debug e      
      MATNORMF = XXX
      RETURN
      END
C
C  /* Function tracemat */
      DOUBLE PRECISION FUNCTION TRACEMAT(DONE,N)  
C
C August-2007 Emmanuel Fromager 
C
C This function calculates the trace of the UNPACKED DONE
C matrix:
C
#include "implicit.h"
      DIMENSION DONE(N,N) 
C
      XXX = 0.0D0
      DO I = 1,N
         XXX= XXX + DONE(I,I)
      END DO
      TRACEMAT = XXX
      RETURN
      END
C
C --- end of sirius/sirmp2.F ---
