#define MAX_NO_OF_PRIMITIVES_PER_BASIS_FUNC 44

struct DistributionSpecStruct_{
  real coeff;           /* Coefficient A */
  real exponent;        /* exponent alfa */
  real extent;
  real centerCoords[3]; /* x0, y0, z0    */
  integer monomialInts[3];  /* nx, ny, nz    */
};
typedef struct DistributionSpecStruct_ DistributionSpecStruct;

#define MAX_NO_OF_CONTR_GAUSSIANS 20

struct ShellSpecStruct_{
  integer noOfContr;
  real coeffList[MAX_NO_OF_CONTR_GAUSSIANS];
  real exponentList[MAX_NO_OF_CONTR_GAUSSIANS];
  real extent;
  real centerCoords[3]; /* x0, y0, z0 */
  integer shellType; /* 1 <-> 's', 2 <-> 'p', 3 <-> 'd' etc */
  integer noOfBasisFuncs;
  integer startIndexInMatrix; /* start index in density matrix  */
};
typedef struct ShellSpecStruct_ ShellSpecStruct; 

struct BasisFuncStruct_{
  integer noOfContr;
  real coeffList[MAX_NO_OF_CONTR_GAUSSIANS];
  real exponentList[MAX_NO_OF_CONTR_GAUSSIANS];
  real extent;
  real centerCoords[3]; /* x0, y0, z0 */
  integer shellType; /* 0 <-> 's', 1 <-> 'p', 2 <-> 'd' etc */
  integer functionNumber; /* -1,0,1 for 'p', -2,-1,0,1,2 for 'd', etc */
  integer noOfSimplePrimitives;
  integer simplePrimitiveIndex;
};
typedef struct BasisFuncStruct_ BasisFuncStruct; 

struct BasisInfoStruct_{
  integer noOfShells;
  ShellSpecStruct* shellList;
  integer noOfBasisFuncs;
  BasisFuncStruct* basisFuncList;
  integer noOfSimplePrimitives;
  DistributionSpecStruct* simplePrimitiveList;
};
typedef struct BasisInfoStruct_ BasisInfoStruct;

integer get_shells(BasisInfoStruct* basisInfo);
integer get_basis_funcs(BasisInfoStruct* basisInfo);
integer get_simple_primitives_all(BasisInfoStruct* basisInfo);
integer get_product_simple_primitives(BasisInfoStruct* basisInfoA, integer iA,
				  BasisInfoStruct* basisInfoB, integer iB,
				  DistributionSpecStruct resultList[],
				  integer maxCount);

