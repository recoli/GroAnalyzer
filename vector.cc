// vector operations
#include <cmath>

double inn_prod(int DIM, double* a, double* b)
{
   double c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
   return c;
}

void out_prod(int DIM, double* a, double* b, double* c)
{
   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];
}

double dist2(int DIM, double* a)
{
   double r2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
   return r2;
}

void norm(int DIM, double* a)
{
   double r = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	for(int i=0; i < DIM; i++)
      a[i] /= r;
}

double calcDih(int DIM, double* a, double* b, double* c, double* d)
{
	double v1[DIM], v2[DIM], v3[DIM];
	for(int i=0; i < DIM; i++)
   {
      v1[i] = b[i] - a[i];
      v2[i] = c[i] - b[i];
      v3[i] = d[i] - c[i];
	}

	double v1x2[DIM], v2x3[DIM];
	out_prod(DIM, v1, v2, v1x2);
	out_prod(DIM, v2, v3, v2x3);

	double ytan = sqrt(dist2(DIM, v2)) * inn_prod(DIM, v1, v2x3);
	double xtan = inn_prod(DIM, v1x2, v2x3);
	double dih = atan2(ytan, xtan) * 180.0 / M_PI;

	return dih;
}

void rotate_axis(int DIM, double* r, double* u,
				 double theta, double* newxyz, double** R)
{
	double rx = r[0];
	double ry = r[1];
	double rz = r[2];

	double ux = u[0];
	double uy = u[1];
	double uz = u[2];

	double cc = cos(theta);
	double ss = sin(theta);

	R[0][0] = cc + ux*ux * (1.0 - cc);
	R[0][1] = ux * uy * (1.0 - cc) - uz * ss;
	R[0][2] = ux * uz * (1.0 - cc) + uy * ss;

	R[1][0] = uy * ux * (1.0 - cc) + uz * ss;
	R[1][1] = cc + uy*uy * (1.0 - cc);
	R[1][2] = uy * uz * (1.0 - cc) - ux * ss;

	R[2][0] = uz * ux * (1.0 - cc) - uy * ss;
	R[2][1] = uz * uy * (1.0 - cc) + ux * ss;
	R[2][2] = cc + uz*uz * (1.0 - cc);

	newxyz[0] = R[0][0] * rx + R[0][1] * ry + R[0][2] * rz;
	newxyz[1] = R[1][0] * rx + R[1][1] * ry + R[1][2] * rz;
	newxyz[2] = R[2][0] * rx + R[2][1] * ry + R[2][2] * rz;
}

void rotate_mat(int DIM, double* r, double** R, double* newxyz)
{
	double rx = r[0];
	double ry = r[1];
	double rz = r[2];

	newxyz[0] = R[0][0] * rx + R[0][1] * ry + R[0][2] * rz;
	newxyz[1] = R[1][0] * rx + R[1][1] * ry + R[1][2] * rz;
	newxyz[2] = R[2][0] * rx + R[2][1] * ry + R[2][2] * rz;
}

