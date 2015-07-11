/* 
   fit to the logistic function
   p(r) = p0 + (p1 - p0) * elg((r - r0) / ksi)

   elg(x) =  1.0 / (1.0 + exp(x))
   [elg(x)]' =  -elg(x) * (1.0 - elg(x))

   beta(1) = p1
   beta(2) = p0
   beta(3) = r0
   beta(4) = ksi
   iter:  iterate or not
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "elg.h"

double elg (double x)
{
   double elg_x;
   elg_x = 1.0 / ( 1.0 + exp(x) ) ;
   return elg_x;
}

Vec_elg elg_fit (double delr, int maxbin, double r[], double pn[],
                 double b_0, double b_1, double b_2, double b_3)
{
   Vec_elg elg_param;
   const double sqrtPI = sqrt(M_PI);

   // damping factor for iteration

   const double damp = 0.8; 

   // variables

   int bin, count;
   int iter;
   double beta[4];

   int i, j, k;
   double a[4][4], b[4], aii, aji;

   // allocate arrays

   double* y = new double [maxbin];
   double** jaco = new double* [maxbin];
   for (bin=0; bin<maxbin; bin++) {
 		jaco[bin] = new double [4];
   }

   // initial guess for beta vector  

   beta[0] = b_0;
   beta[1] = b_1;
   beta[2] = b_2;
   beta[3] = b_3;

   // iteratively solve the non-linear least square fitting  
 
   count = 0;
   iter = 1;
   while (iter==1) 
   {

      // calculate y(bin) and jacobian matrices 
      // pres(r) = pres_g + ( pres_l - pres_g ) * elg( ( r - r0 ) / ksi ) 
      
      // beta[0] ---> pres_l
      // beta[1] ---> pres_g
      // beta[2] ---> r0
      // beta[3] ---> ksi

      for ( bin=0; bin<maxbin; bin++ ) 
	  {
         y[bin]       = beta[1]
                      + ( beta[0] - beta[1] )  
                      * elg( ( r[bin]-beta[2] ) / beta[3] );
         jaco[bin][0] = elg( ( r[bin]-beta[2] ) / beta[3] );
         jaco[bin][1] = 1.0 - elg( ( r[bin]-beta[2] ) / beta[3] );
         jaco[bin][2] = ( beta[0] - beta[1] )  
                      * (-1.0) * exp( ( r[bin]-beta[2] ) / beta[3] )
                      / pow( 1.0 + exp( ( r[bin]-beta[2] ) / beta[3] ), 2.0 )
                      * ( -1.0 / beta[3] );
         jaco[bin][3] = ( beta[0] - beta[1] )  
                      * (-1.0) * exp( ( r[bin]-beta[2] ) / beta[3] )
                      / pow( 1.0 + exp( ( r[bin]-beta[2] ) / beta[3] ), 2.0 )
                      * ( r[bin]-beta[2] ) * ( -1.0/pow(beta[3],2.0) );
      }

      // determine coefficients of the linear equation system  

      for ( i=0; i<4; i++ ) 
      {
         for ( j=0; j<4; j++ ) 
         {
            a[i][j] = 0.0;
         }
         b[i] = 0.0;
      }

      // summation starts from bin=10 to avoid bad statistics near the center 

      for ( bin=10; bin<maxbin; bin++ ) 
      {
         for ( i=0; i<4; i++ ) 
         {
            for ( j=0; j<4; j++ ) 
            {
               a[i][j] += jaco[bin][i] * jaco[bin][j];
            }
            b[i] += jaco[bin][i] * ( pn[bin] - y[bin] );
         }
      }

      // use Gaussian elimination  

      for ( i=0; i<4; i++ ) 
      {
         aii = a[i][i];
         for ( j=0; j<4; j++ ) 
         {
            a[i][j] /= aii;
         }
         b[i] /= aii;
         for ( j=i+1; j<4; j++ ) 
         {
            aji = a[j][i];
            for ( k=0; k<4; k++ ) 
            {
               a[j][k] -= aji * a[i][k];
            }
            b[j] -= aji * b[i];
         }
      }
      for ( i=4-1; i>=0; i-- ) 
      {
         b[i] /= a[i][i];
         for ( j=i-1; j>=0; j-- ) 
         {
            b[j] -= a[j][i] * b[i];
         }
      }

      // generate new beta vector and check the change  
 
      for ( i=0; i<4; i++ ) 
      {
         beta[i] += b[i] * damp;
         if ( fabs( b[i]/beta[i] ) > 1.0e-4 ) 
         {
            iter = 1;
            break;
         } 
         else 
         {
            iter = 0;
            continue;
         }
      }

      count++;
      if ( count > 1000 ) 
      {
         printf("<> Error: Iteration exceeds 1000 in elg_fit !\n");
         exit(1);
      }

   // end of loop (iter)  

   }

   // free arrays

   free(y);
   free(jaco);

   elg_param.liquid = beta[0];
   elg_param.gas = beta[1];
   elg_param.r0 = beta[2];
   elg_param.ksi = beta[3];

   return elg_param;
}

double integr_r3 ( double dR, int maxBin, double* r, double* pn )
{
   double *dpn = new double [maxBin];

   double integr;
   integr = 0.0;

   int bin;

   for ( bin=0; bin<maxBin; bin++ ) 
   {

      //Numerical differentiation
      //http://en.wikipedia.org/wiki/Finite_difference_coefficients
      //                           0     1     2     3     4     5     6
      //Forward coefficients:   -49/20   6  -15/2  20/3 -15/4   6/5  -1/6
      //Backward coefficients:   49/20  -6   15/2 -20/3  15/4  -6/5   1/6
      //                          -3    -2    -1     0     1     2     3
      //Central coefficients:    -1/60  3/20 -3/4    0    3/4  -3/20  1/60

      if ( bin < 3 )
      {
         dpn[bin] = -49.0/20.0*pn[bin]
                    + 6.0*pn[bin+1] - 15.0/2.0*pn[bin+2] + 20.0/3.0*pn[bin+3]
                    - 15.0/4.0*pn[bin+4] + 6.0/5.0*pn[bin+5] - 1.0/6.0*pn[bin+6];
         dpn[bin] /= dR;
      }
      else if ( bin > maxBin-4 )
      {
         dpn[bin] = 49.0/20.0*pn[bin]
                    - 6.0*pn[bin-1] + 15.0/2.0*pn[bin-2] - 20.0/3.0*pn[bin-3]
                    + 15.0/4.0*pn[bin-4] - 6.0/5.0*pn[bin-5] + 1.0/6.0*pn[bin-6];
         dpn[bin] /= dR;
      }
      else
      {
         dpn[bin] = -1.0/60.0*pn[bin-3] + 3.0/20.0*pn[bin-2] - 3.0/4.0*pn[bin-1]
                    + 3.0/4.0*pn[bin+1] - 3.0/20.0*pn[bin+2] + 1.0/60.0*pn[bin+3];
         dpn[bin] /= dR;
      }

      integr += pow(r[bin],3.0) * dpn[bin] * dR;
   }

   return integr;
}

