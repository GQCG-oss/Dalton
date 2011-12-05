#define MAX_NO_OF_PRIMITIVES_PER_BASIS_FUNC 44

struct DistributionSpecStruct_{
  real coeff;           /* Coefficient A */
  real exponent;        /* exponent alfa */
  real extent;
  real centerCoords[3]; /* x0, y0, z0    */
  int monomialInts[3];  /* nx, ny, nz    */
};
typedef struct DistributionSpecStruct_ DistributionSpecStruct;

#define MAX_NO_OF_CONTR_GAUSSIANS 20

struct ShellSpecStruct_{
  int noOfContr;
  real coeffList[MAX_NO_OF_CONTR_GAUSSIANS];
  real exponentList[MAX_NO_OF_CONTR_GAUSSIANS];
  real extent;
  real centerCoords[3]; /* x0, y0, z0 */
  int shellType; /* 1 <-> 's', 2 <-> 'p', 3 <-> 'd' etc */
  int noOfBasisFuncs;
  int startIndexInMatrix; /* start index in density matrix  */
};
typedef struct ShellSpecStruct_ ShellSpecStruct; 

struct BasisFuncStruct_{
  int noOfContr;
  real coeffList[MAX_NO_OF_CONTR_GAUSSIANS];
  real exponentList[MAX_NO_OF_CONTR_GAUSSIANS];
  real extent;
  real centerCoords[3]; /* x0, y0, z0 */
  int shellType; /* 0 <-> 's', 1 <-> 'p', 2 <-> 'd' etc */
  int functionNumber; /* -1,0,1 for 'p', -2,-1,0,1,2 for 'd', etc */
  int noOfSimplePrimitives;
  int simplePrimitiveIndex;
};
typedef struct BasisFuncStruct_ BasisFuncStruct; 

struct BasisInfoStruct_{
  int noOfShells;
  ShellSpecStruct* shellList;
  int noOfBasisFuncs;
  BasisFuncStruct* basisFuncList;
  int noOfSimplePrimitives;
  DistributionSpecStruct* simplePrimitiveList;
};
typedef struct BasisInfoStruct_ BasisInfoStruct;

int get_shells(BasisInfoStruct* basisInfo);
int get_basis_funcs(BasisInfoStruct* basisInfo);
int get_simple_primitives_all(BasisInfoStruct* basisInfo);
int get_product_simple_primitives(BasisInfoStruct* basisInfoA, int iA,
				  BasisInfoStruct* basisInfoB, int iB,
				  DistributionSpecStruct resultList[],
				  int maxCount);

