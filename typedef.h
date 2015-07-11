#define DIM 3

// maximal length of string
#define MAX_STR_LEN 256

// vector
typedef struct { 
	double x, y, z; 
} VecR;

// tensor
typedef struct { 
	double xx,xy,xz, yx,yy,yz, zx,zy,zz; 
} TensR;

// atomic information
typedef struct { 
	VecR r, v, f, a; 
	double m, q, sig, eps;
	long int resID, atomID, typeID;
	char resName[5], atomName[5], typeName[5]; 
} Atom;

