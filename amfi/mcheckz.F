      logical function mcheckz(m1,m2,m3,m4)
cbs   makes a check, if there is an interaction inbetween cartesian functions 
cbs   with m-values m1-m4 
      integer*4 m1,m2,m3,m4,int12a,int12b,
     *int34a,int34b
      mcheckz=.true.   
      int12a=m1+m2
      int12b=-m1+m2
      int34a=m3+m4
      int34b=-m3+m4
cbs   lots of checks 
      if (iabs(int12a+int34a).eq.0) return
      if (iabs(int12a-int34a).eq.0) return
      if (iabs(int12b+int34b).eq.0) return
      if (iabs(int12b-int34b).eq.0) return
      if (iabs(int12a+int34b).eq.0) return
      if (iabs(int12a-int34b).eq.0) return
      if (iabs(int12b+int34a).eq.0) return
      if (iabs(int12b-int34a).eq.0) return
      mcheckz=.false.
      return 
      end 
