cbs   irreducible representation of the cartesian functions in D2H
cbs   order taken from Tables for Group theory:
cbs   Atkins, Child and Phillips   Oxford University Press 1970
cbs   1. AG: only even powers (0,0,0)
cbs   2. B1G: (1,1,0)    L_z
cbs   3. B2G: (1,0,1)    L_y
cbs   4. B3G: (0,1,1)    L_x
cbs   5. AU:  (1,1,1)
cbs   6. B1U: (0,0,1)
cbs   7. B2U: (0,1,0)
cbs   8. B3U: (1,0,0)
      integer ipow2ired,iredorder,incrLM,shiftIRED,shiftIRIR
      common /ireduceD2H/ numbofsym, ! number of symmetries
     *ipow2ired(0:1,0:1,0:1),! gives IR by checking powers
     *iredorder(8),  ! maybe reordering  of IRs is necessary
     *iredorderinv(8),
     *nfunctions(8,0:Lmax), !number of functions per L and IR
     *nfunctperIRED(8),   ! number of functions per IR
     *incrLM(-Lmax:Lmax,0:Lmax), ! shift of orbitalnumber in IR for L,M
     *shiftIRED(8), !shift to get to absolute number from relative number in IR
     *iredLM(-Lmax:Lmax,0:Lmax), ! IR for L and M
     *shiftIRIR(36),  ! shift for (IR1,IR2)-block (IR1<=IR2)
     *Loffunction(Mxcart), !gives L value of cartesian function
     *Moffunction(Mxcart), !gives M value of cartesian function
     *Iredoffunction(Mxcart),! give IRED of cartesian function
     *Iredoffunctnew(Mxcart),! give IRED of cartesian function incl. add. functions
     *IRwithLX(8), ! gives IR interacting by L_X
     *IRwithLY(8), ! gives IR interacting by L_y
     *IRwithLZ(8), ! gives IR interacting by L_z
     *mult(8,8),
     *itotalperIR(8), ! total number of functions per IR
     *nmbMperIRL(8,0:Lmax), !number of M-values in an IR for an L-value
     *numbofcart  ! number of cartesian functions
cbs
