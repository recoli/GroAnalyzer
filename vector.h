double inn_prod(int dimension, double* a, double* b);
void   out_prod(int dimension, double* a, double* b, double* c);
double dist2(int dimension, double* a);
void   norm(int dimension, double* a);
double calcDih(int dimension, double* a, double* b, double* c, double* d);
void   rotate_axis(int dimension, double* r, double* u,
				   double theta, double* newxyz, double** R);
void   rotate_mat(int dimension, double* r, double** R, double* newxyz);
