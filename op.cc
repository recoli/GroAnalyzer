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
	for (long int frame = 0; frame < nFrames; ++ frame)
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


		// compute order parameter according to M.R.Wilson
		// Journal of Molecular Liquids 68 (1996) 23-31

		// long-axis directions
		// normal directions
		// center of mass
		double **axis = new double* [nMols];
		double **norm = new double* [nMols];
		double **cent = new double* [nMols];

		// initialize the ordering tensor Q
		// for long-axis(a) and normal(n) directions
		double **Qa = new double* [DIM];
		double **Qn = new double* [DIM];
		for (int row=0; row < DIM; row++)
		{
			Qa[row] = new double [DIM];
			Qn[row] = new double [DIM];
			for (int col=0; col < DIM; col++)
			{
				Qa[row][col] = 0.0;
				Qn[row][col] = 0.0;
			}
		}

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

			// add to ordering tensor Q
			for (int row=0; row < DIM; row++)
			{
				Qa[row][row] += -0.5;
				Qn[row][row] += -0.5;
				for (int col=0; col < DIM; col++)
				{
					Qa[row][col] += 1.5 * axis[mol][row] * axis[mol][col];
					Qn[row][col] += 1.5 * norm[mol][row] * norm[mol][col];
				}
			}
		}

		// average Q over molecules
		for (int row=0; row < DIM; row++)
		{
			for (int col=0; col < DIM; col++)
			{
				Qa[row][col] /= (double)nMols;
				Qn[row][col] /= (double)nMols;
			}
		}

		// diagonalize ordering tensor Q
		// get order parameter as the largest eigenvalue
		gsl_matrix *Qa_mat = gsl_matrix_alloc(DIM, DIM);
		gsl_matrix *Qn_mat = gsl_matrix_alloc(DIM, DIM);
		for (int row=0; row < DIM; row++)
		{
			for (int col=0; col < DIM; col++)
			{
				gsl_matrix_set(Qa_mat, row, col, Qa[row][col]);
				gsl_matrix_set(Qn_mat, row, col, Qn[row][col]);
			}
		}


		gsl_vector *eval = gsl_vector_alloc (DIM);
		gsl_matrix *evec = gsl_matrix_alloc (DIM, DIM);
 
		gsl_matrix *Q_mat_cp = gsl_matrix_alloc (DIM, DIM);
		gsl_eigen_symmv_workspace *wks = gsl_eigen_symmv_alloc (DIM);


		gsl_matrix_memcpy (Q_mat_cp, Qa_mat);
		gsl_eigen_symmv (Q_mat_cp, eval, evec, wks);
		gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
		double order_param_a = gsl_vector_get(eval, 2);


		gsl_matrix_memcpy (Q_mat_cp, Qn_mat);
		gsl_eigen_symmv (Q_mat_cp, eval, evec, wks);
		gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
		double order_param_n = gsl_vector_get(eval, 2);


		gsl_eigen_symmv_free (wks);
		gsl_matrix_free (Q_mat_cp);


		printf("%ld  %f  %f\n", frame, order_param_a, order_param_n);

	}
	fclose(groFile);

	return 0;
}  

