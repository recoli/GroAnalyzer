#include <cmath>  
#include <cstdio>  
#include <cstdlib>  

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>

#include "typedef.h"  
#include "vector.h"  
#include "file.h"  

using namespace std;
  
int main()  
{  
	// histogram for the distribution of angles
	int *hist = new int [361];
	int maxbin = 360;
	int count  = 0;
	for (int bin = 0; bin <= maxbin; ++ bin)
	{
		hist[bin] = 0;
	}

	// pbc box
	double *box = new double [DIM];


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
		++ imass;
	}
	fclose(massFile);


	// read atomic information from gro file
	Atom *atom = new Atom [nAtoms];

	FILE *groFile;
	groFile = fopen("traj.gro", "r") ;
	//long int nFrames = 2501;
	for (long int frame = 0; frame < nFrames; ++ frame)
	{
		readGRO(groFile, nAtoms, atom, box);

		// only calculate distribution from the last 40% trajectory
		if (frame < nFrames*3/5)
		{
			continue;
		}

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

		if (frame == nFrames - 1)
		{
			printf("#nAtoms=%ld, nMols=%ld, apm=%ld\n",
					nAtoms, nMols, apm);
		}

		// compute order parameter according to M.R.Wilson
		// Journal of Molecular Liquids 68 (1996) 23-31

		// long-axis directions
		// normal directions
		// center of mass
		double **axis = new double* [nMols];
		double **norm = new double* [nMols];
		double **cent = new double* [nMols];

		// compute long axis for each molecule
		for (long int mol=0; mol < nMols; mol++)
		{
			axis[mol] = new double [DIM];
			norm[mol] = new double [DIM];
			cent[mol] = new double [DIM];

			// compute center-of-mass
			cent[mol][0] = 0.0;
			cent[mol][1] = 0.0;
			cent[mol][2] = 0.0;
			double msum = 0.0;

			for (long int at = 0; at < apm; ++ at)
			{
				double m = mass[at];
				cent[mol][0] += atom[at + mol * apm].r.x * m;
				cent[mol][1] += atom[at + mol * apm].r.y * m;
				cent[mol][2] += atom[at + mol * apm].r.z * m;
				msum += m;
			}
			cent[mol][0] /= msum;
			cent[mol][1] /= msum;
			cent[mol][2] /= msum;

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

			for (long int at=0; at < apm; at++)
			{
				double *r = new double [DIM];
				r[0] = atom[at+mol*apm].r.x - cent[mol][0];
				r[1] = atom[at+mol*apm].r.y - cent[mol][1];
				r[2] = atom[at+mol*apm].r.z - cent[mol][2];
				double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
				double m = mass[at];

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

			gsl_matrix *I_mat_cp = gsl_matrix_alloc (DIM, DIM);
			gsl_matrix_memcpy(I_mat_cp, I_mat);
			gsl_vector *eval = gsl_vector_alloc (DIM);
			gsl_matrix *evec = gsl_matrix_alloc (DIM, DIM);
 
			gsl_eigen_symmv_workspace *wks = gsl_eigen_symmv_alloc (DIM);
			gsl_eigen_symmv (I_mat_cp, eval, evec, wks);
			gsl_eigen_symmv_free (wks);
			gsl_matrix_free (I_mat_cp);
 
			gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);

			// long-axis direction: 
			// eigenvector associated with the smallest eigenvalue
			axis[mol][0] = gsl_matrix_get(evec, 0, 0);
			axis[mol][1] = gsl_matrix_get(evec, 1, 0);
			axis[mol][2] = gsl_matrix_get(evec, 2, 0);

			// normal direction: 
			// eigenvector associated with the largest eigenvalue
			norm[mol][0] = gsl_matrix_get(evec, 0, 2);
			norm[mol][1] = gsl_matrix_get(evec, 1, 2);
			norm[mol][2] = gsl_matrix_get(evec, 2, 2);
		}

		// calculate slip angles between adjacent molecules
		for (long int imol=0; imol < nMols; imol++)
		{
			for (long int kmol=imol+1; kmol < nMols; kmol++)
			{
				// vector connecting center-of-mass
				double *vec = new double [DIM];
				vec[0] = cent[imol][0] - cent[kmol][0];
				vec[1] = cent[imol][1] - cent[kmol][1];
				vec[2] = cent[imol][2] - cent[kmol][2];

				// apply PBC
				if (vec[0] >  0.5 * box[0]) { vec[0] -= box[0]; }
				if (vec[0] < -0.5 * box[0]) { vec[0] += box[0]; }
				if (vec[1] >  0.5 * box[1]) { vec[1] -= box[1]; }
				if (vec[1] < -0.5 * box[1]) { vec[1] += box[1]; }
				if (vec[2] >  0.5 * box[2]) { vec[2] -= box[2]; }
				if (vec[2] < -0.5 * box[2]) { vec[2] += box[2]; }

				// apply PBC
				double dvec = sqrt(dist2(DIM, vec));
				double d_i  = fabs(inn_prod(DIM, norm[imol], vec));
				double d_k  = fabs(inn_prod(DIM, norm[kmol], vec));

				// criteria as adjacent molecules: 
				// the vertial distance between two molecules
				// is smaller than 5 Angstrom (0.5 nm)
				if (d_i > 0.5 || d_k > 0.5) { continue; }

				double theta_i = acos(d_i / dvec) / M_PI * 180.0;
				double theta_k = acos(d_k / dvec) / M_PI * 180.0;

				int bin_i = (int)(90.0 - theta_i + 0.5);
				int bin_k = (int)(90.0 - theta_k + 0.5);
				++ hist[bin_i];
				++ hist[bin_k];
				++ count;
				++ count;
			}
		}
	}
	fclose(groFile);

	for (int bin = 0; bin <= 360; ++ bin)
	{
		printf("%d  %f\n", bin, (double)hist[bin] / (double)count);
	}

	return 0;
}  

