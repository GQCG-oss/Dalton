*
*     Parameter definitions
*     Note, these parameters define the size of the problems which
*     can be treated by the RASSCF program. Changing any of these 
*     requires, taht the program is totaly recompiled and linked.
*
      Parameter ( mxSym  = 8       )
*...  max number of irreps
      Parameter ( mxBas  = 255     )
*...  max number of basis func. per irrep      
      Parameter ( mxOrb  = 2*mxBas )
*...  max total number of basis func.          
      Parameter ( mxAct  = 50      )
*...  max number of active orbitals      
      Parameter ( mxRoot = 25      )
*...  max number of roots                      
      Parameter ( mxRef  = 5       )
*...  max number of reference configurations in root selectioning
      Parameter ( mxIter = 50      )
*...  max number of macro iterations           
      Parameter ( mxCiIt = 202     )
*...  max number of micro iterations for the CI section
      Parameter ( mxSxIt = 50      )
*...  max number of micro iterations for the SX section
      Parameter ( mxTit  = 10      )
*...  max number of title lines
*
