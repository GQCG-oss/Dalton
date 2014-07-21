
#ifndef __DALTON_CONSTS__
#define __DALTON_CONSTS__

//DALTON PARAMETERS
//static const int MXSHEL = 750;
static const int MXSHEL = 1000;
static const int MXPRIM = 8000;
static const int MXCONT = 35;
//static const int MXPOINTS = 65536;

static const int maxJ = 36; // according to DALTON


// two extra restrictions in this code which are likely to hold true in any realistic calculation
static const int LMAX = 5;
static const int maxK = 32;


//OTHER CONSTANTS
static const int nmS[] = {1,3,5,7,9,11,   4};
static const int nmC[] = {1,3,6,10,15,21, 4};

static const double PI = 3.1415926535897932384626433832;
static const int CACHE_LINE_SIZE = 64;

static const int MMAX = (LMAX+1)*(LMAX+2)/2; //(h+1)*;
static const int LSP  = LMAX + 1; //code for SP GTOs

static const int DOUBLES_PER_CACHE_LINE = 8;

typedef long int lint;

#endif
