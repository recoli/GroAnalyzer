#include <cstdio>  
#include <cstdlib>  

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>

#include "typedef.h"  
#include "vector.h"  
#include "file.h"  
#include "elg.h"  

#define MAXBIN 500

using namespace std;
  
int main()  
{  
	double* avR = new double [MAXBIN];
	double* avDens = new double [MAXBIN];
	double* inv_dV = new double [MAXBIN];

	// count frames
	double  count = 0;

	double* box = new double [DIM];

	double dR = 0.002; // nm
	for(int bin=0; bin < MAXBIN; bin++)
	{
		avR[bin] = ((double)bin + 0.5) * dR;
		avDens[bin] = 0.0;
		inv_dV[bin] = 1.0 / ( M_PI * ( pow((double)(bin+1) * dR, 2.0) -
										pow((double)bin * dR, 2.0) ) );
	}

	FILE *lenFile;
	lenFile = fopen("length.txt", "w") ;
	if (NULL == lenFile)
	{
		printf("Error: cannot write to file length.txt!\n");
		exit(1);
	}

	FILE *xyzFile;
	xyzFile = fopen("new-traj.xyz", "w") ;
	if (NULL == xyzFile)
	{
		printf("Error: cannot write to file new-traj.xyz!\n");
		exit(1);
	}


	// get number of atoms and frames
	long int nAtoms, nFrames;
	readNAtomsFrames("traj.gro", &nAtoms, &nFrames);


	// read atomic mass
	double *mass = new double [nAtoms];
	FILE *massFile;
	massFile = fopen("masses.dat", "r") ;
	if (NULL == massFile)
	{
		printf("Error: cannot open file masses.dat!\n");
		exit(1);
	}

	int imass = 0;
	while(fscanf(massFile, "%lf", &mass[imass])!=EOF)
	{
		imass++;
	}
	fclose(massFile);


	// read atomic information from gro file
	Atom *atom = new Atom [nAtoms];

	FILE *groFile;
	groFile = fopen("traj.gro", "r") ;
	if (NULL == groFile)
	{
		printf("Error: cannot open file traj.gro!\n");
		exit(1);
	}

	for (long int frame=0; frame < nFrames; frame ++)
	{
		readGRO(groFile, nAtoms, atom, box);

		// get number of molecules
		// and number of atoms per molecule
		long int nMols, apm;
		for (long int i = 0; i < nAtoms - 1; ++ i)
		{
			if (atom[i].resID != atom[i+1].resID)
			{
				apm = i + 1;
				break;
			}
		}

		if(nAtoms % apm != 0)
		{
			printf("Error: nAtoms(%ld) %% apm(%ld) is not zero!\n",
					nAtoms, apm);
			exit(1);
		}
		nMols = nAtoms / apm;


		// apply PBC
		for (long int i=1; i < nAtoms; i++)
		{
			if      (atom[i].r.x - atom[0].r.x >  0.5 * box[0]) { atom[i].r.x -= box[0]; }
			else if (atom[i].r.x - atom[0].r.x < -0.5 * box[0]) { atom[i].r.x += box[0]; }
			if      (atom[i].r.y - atom[0].r.y >  0.5 * box[1]) { atom[i].r.y -= box[1]; }
			else if (atom[i].r.y - atom[0].r.y < -0.5 * box[1]) { atom[i].r.y += box[1]; }
			if      (atom[i].r.z - atom[0].r.z >  0.5 * box[2]) { atom[i].r.z -= box[2]; }
			else if (atom[i].r.z - atom[0].r.z < -0.5 * box[2]) { atom[i].r.z += box[2]; }

			long int mol   = i / apm;
			long int start = mol * apm;
			long int end   = (mol+1) * apm - 1;
		}
		for (long int mol = 1; mol < nMols; ++ mol)
		{
			long int start = mol * apm;
			long int end   = (mol+1) * apm - 1;
			for (long int i = start+1; i <= end; ++i)
			{
				if      (atom[i].r.x - atom[start].r.x >  0.5 * box[0]) { atom[i].r.x -= box[0]; }
				else if (atom[i].r.x - atom[start].r.x < -0.5 * box[0]) { atom[i].r.x += box[0]; }
				if      (atom[i].r.y - atom[start].r.y >  0.5 * box[1]) { atom[i].r.y -= box[1]; }
				else if (atom[i].r.y - atom[start].r.y < -0.5 * box[1]) { atom[i].r.y += box[1]; }
				if      (atom[i].r.z - atom[start].r.z >  0.5 * box[2]) { atom[i].r.z -= box[2]; }
				else if (atom[i].r.z - atom[start].r.z < -0.5 * box[2]) { atom[i].r.z += box[2]; }
			}
		}


		// get long axis of the assembly according to M.R.Wilson
		// Journal of Molecular Liquids 68 (1996) 23-31

		double *axis = new double [nMols];

		// compute center-of-mass
		double *com = new double [DIM];
		com[0] = 0.0;
		com[1] = 0.0;
		com[2] = 0.0;
		double msum = 0.0;

		for (long int i=0; i < nAtoms; i++)
		{
			double m = mass[i % apm];
			com[0] += atom[i].r.x * m;
			com[1] += atom[i].r.y * m;
			com[2] += atom[i].r.z * m;
			msum += m;
		}
		com[0] /= msum;
		com[1] /= msum;
		com[2] /= msum;

		// shift the system to COM
		for (long int i=0; i < nAtoms; i++)
		{
			atom[i].r.x -= com[0];
			atom[i].r.y -= com[1];
			atom[i].r.z -= com[2];
		}

		// compute inertia tensor
		double **I = new double* [DIM];
		for (int row=0; row < DIM; row++)
		{
			I[row] = new double [DIM];
			for (int col=0; col < DIM; col++)
			{
				I[row][col] = 0.0;
			}
		}

		for (long int i=0; i < nAtoms; i++)
		{
			// Note: coordinates have been shifted to COM
			double *r = new double [DIM];
			r[0] = atom[i].r.x;
			r[1] = atom[i].r.y;
			r[2] = atom[i].r.z;

			double r2 = dist2(DIM, r);
			double m = mass[i % apm];

			for (int row=0; row < DIM; row++)
			{
				I[row][row] += m * r2;
				for (int col=0; col < DIM; col++)
				{
					I[row][col] += -m * r[row] * r[col];
				}
			}
		}

		// diagonalize inertia tensor
		// get long axis as the eigenvector associated
		// with the smallest eigenvalue
		gsl_matrix *I_mat = gsl_matrix_alloc(DIM, DIM);
		for (int row=0; row < DIM; row++)
		{
			for (int col=0; col < DIM; col++)
			{
				gsl_matrix_set(I_mat, row, col, I[row][col]);
			}
		}

		gsl_vector *eval = gsl_vector_alloc (DIM);
		gsl_matrix *evec = gsl_matrix_alloc (DIM, DIM);

		gsl_matrix *I_mat_cp = gsl_matrix_alloc (DIM, DIM);
		gsl_eigen_symmv_workspace *wks = gsl_eigen_symmv_alloc (DIM);

		gsl_matrix_memcpy(I_mat_cp, I_mat);
		gsl_eigen_symmv (I_mat_cp, eval, evec, wks);

		gsl_eigen_symmv_free (wks);
		gsl_matrix_free (I_mat_cp);
 
		gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);

		axis[0] = gsl_matrix_get(evec, 0, 0);
		axis[1] = gsl_matrix_get(evec, 1, 0);
		axis[2] = gsl_matrix_get(evec, 2, 0);

		// rotate the long axis to x-axis
		// u_axis: unit vector to rotate around
		// theta: angle to rotate

		double *x_axis = new double [DIM];
		x_axis[0] = 1.0;
		x_axis[1] = 0.0;
		x_axis[2] = 0.0;

		double *u_axis = new double [DIM];
		out_prod(DIM, axis, x_axis, u_axis);
		norm(DIM, u_axis);

		double theta = acos(inn_prod(DIM, axis, x_axis) / 
							sqrt(dist2(DIM, axis)));

		double *new_axis = new double [DIM];

		double **R = new double* [DIM];
		for (int row=0; row < DIM; row++)
		{
			R[row] = new double [DIM];
		}

		rotate_axis(DIM, axis, u_axis, theta, new_axis, R);

		if(fabs(new_axis[1]) > 1.0e-06 ||
			fabs(new_axis[2]) > 1.0e-06)
		{
			printf("Error in rotating long axis!\n");
			exit(1);
		}

		// also rotate the assembly

		double xmin = 0.0;
		double xmax = 0.0;

		double* r = new double [DIM];
		double* new_r = new double [DIM];

		fprintf(xyzFile, "%ld\n\n", nAtoms);
		for (long int i=0; i < nAtoms; i++)
		{
			r[0] = atom[i].r.x;
			r[1] = atom[i].r.y;
			r[2] = atom[i].r.z;

			rotate_mat(DIM, r, R, new_r);
			fprintf(xyzFile, "%s %f %f %f\n", 
					atom[i].atomName, new_r[0]*10.0, new_r[1]*10.0, new_r[2]*10.0);

			if(new_r[0] < xmin) xmin = new_r[0];
			if(new_r[0] > xmax) xmax = new_r[0];
		}

		fprintf(lenFile, "%ld %f %f %f\n", frame, xmax-xmin, xmin, xmax);

		// radial density
		// only calculate from the last 40% trajectory
		if(frame >= nFrames*3/5)
		{
			for (long int i=0; i < nAtoms; i++)
			{
				r[0] = atom[i].r.x;
				r[1] = atom[i].r.y;
				r[2] = atom[i].r.z;

				rotate_mat(DIM, r, R, new_r);

				double radius = sqrt(new_r[1]*new_r[1] + new_r[2]*new_r[2]);
				int bin = (int)(radius / dR);

				double m = mass[i % apm];
				avDens[bin] += inv_dV[bin] * m / (6.022141 * 100.0) / (xmax-xmin);
			}

			++ count;
		}

		delete(r);
		delete(new_r);
	}

	fclose(groFile);
	fclose(xyzFile);
	fclose(lenFile);

	FILE *densFile;
	densFile = fopen("dens.txt", "w") ;
	if (NULL == densFile)
	{
		printf("Error: cannot write to file dens.txt!\n");
		exit(1);
	}

	for(int bin=0; bin < MAXBIN; bin++)
	{
		double r = ((double)bin + 0.5) * dR;
		avDens[bin] /= (double)count;
		fprintf(densFile, "%e %e\n", r, avDens[bin]);
	}
	fclose(densFile);

	Vec_elg  elg_param;
	elg_param = elg_fit(dR, MAXBIN-10, avR, avDens,
						2.5, 0.0, 0.4, 0.05);
	printf("Radius = %f\n", elg_param.r0);

	return 0;
}  

