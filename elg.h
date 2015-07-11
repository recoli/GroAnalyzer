typedef struct {
   double liquid, gas, r0, ksi;
} Vec_elg;

double elg (double x);
Vec_elg elg_fit (double delr, int maxbin, double r[], double pn[],
                 double b_0, double b_1, double b_2, double b_3);
double integr_r3 ( double dR, int maxBin, double* r, double* pn );
