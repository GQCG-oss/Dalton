      parameter (Lpowmax=4) 
      dimension ixyzpow(3*(Lpowmax+1)*(Lpowmax+1)) !  
      data ixyzpow /
cbs   the ones and zeros stand four odd and even powers of x,y,z 
cbs   if you want to go higher than l=4, you have to look up 
cbs   the powers yourself, and add them to the table 
     *0,0,0,                                 ! s-functions
     *0,1,0, 0,0,1, 1,0,0,                   ! p-functions
     *1,1,0, 0,1,1, 0,0,0,  1,0,1, 0,0,0,    ! d-functions
     *0,1,0, 1,1,1, 0,1,0,  0,0,1, 1,0,0,    ! f-functions
     *0,0,1, 1,0,0,                          ! f-functions
     *1,1,0, 0,1,1, 1,1,0,  0,1,1, 0,0,0,    ! g-functions
     *1,0,1, 0,0,0, 1,0,1,  0,0,0            ! g-functions
     */

