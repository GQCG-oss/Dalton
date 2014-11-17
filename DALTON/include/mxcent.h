!
!     File: mxcent.h
!
!     MXCENT_QM = max number of QM nuclei + point charges + ghost orbital centers
!     MXCENT = MXCENT_QM + max number of MM nuclei in QM3 model
!              (number of nuclei in QMMM model is separately allocated in qmmm.h)
!
!     IF you change MXCENT you need to rebuild with "make".
!
!     In case of QM3 MXCENT_QM is used to allocate memory in herrdn.F.
!     To run a QM3 calculation in most cases MXCENT will have to be
!     around 2000 - 3000!!! Remember to set MXQM3 = MXCENT_QM in qm3.h!!!
!
      INTEGER MXCENT_QM, MXCENT, MXCOOR
      PARAMETER (MXCENT_QM = 500, MXCENT = 500, MXCOOR = 3*MXCENT)
! -- end of mxcent.h --
