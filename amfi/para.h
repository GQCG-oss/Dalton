      parameter (Lmax=4)   ! max. angular momentum of basis functions  DO NOT INCREASE
cbs                          TO MORE THAN FOUR  !!!!!!!!!!!!!!!!!!!
cbs                          if you do, you will have to edit the ixyzpow array by hand...............
cbs                          (in datapow.h)
cbs                          and the overlap coefficients for HERMIT in 
cbs                          contandmult.f
      parameter (MxprimL=30) ! max. of primitives per angular momentum
      parameter (MxcontL=14) ! max. of contracted functions per angular momentum
      parameter (Mxcart=180) ! max. of contracted cartesian functions        
      parameter (ndfmx=4*Lmax+4) ! dimension of precomputed double factorials     
      parameter (ncffline=9) ! max number of contraction coeffs in input line
      parameter (icont4= 5000000) ! size of array with contracted integrals
      parameter (mxangint=900000) ! size of array with atomic spin orbit integrals
      parameter (mxicart=400000) ! size of array with cartesian integrals
      parameter (mxfrac= 4000000) ! maximum size of precomputed fractions for cfunct 
      parameter (inpfile= 38) ! the input-file MNF.INP  
