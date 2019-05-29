



#include <mpi.h>
#include <omp.h>
#include <math.h>

#include "Resolution.h"

extern int X_np;

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

void Flux_X
(
// ============================================================================ //
int myid,

double Ep,

double e,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*MR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ML1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*J_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*inFx1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*EpX)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	double rho,U,V,W,VV,P,C,T,h,H;
	double u,v,w;
	double temp;

	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
	double XIx,XIy,XIz,ETx,ETy,ETz,ZTx,ZTy,ZTz;
	double _rho,_u,_v,_w,iU,_V,_W,__U,__V,__W,_VV,iP,_T,iC,_H;
	double rho_,u_,v_,w_,U_,V_,W_,U__,V__,W__,VV_,P_,T_,C_,H_;
	double _U_,_V_,_W_;
	double dU1,dU2,dU3,dU4,dU5;
	double thetap,Unp,Cnp,theta,Cnt,deltaU,deltaP,beta;
	double Fav1,Fav2,Fav3,Fav4,Fav5;

/**** MUSCL 5th-order ****/
// ============================================================================ //


//// ============================================ ////
		if (myid ==0) istart = 4;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;				  ////
//// ============================================ ////



	for (i = istart; i <= iend; i++) {
	
	#pragma omp parallel for private(k,temp)
		for ( j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {
				
				temp = 1./(2/J[i-2][j][k]-13/J[i-1][j][k]+47/J[i][j][k]+27/J[i+1][j][k]-3/J[i+2][j][k]);

				ML1[i][j-1][k-1] = temp*(2*U1_[i-2][j][k]-13*U1_[i-1][j][k]+47*U1_[i][j][k]+27*U1_[i+1][j][k]-3*U1_[i+2][j][k]);
				ML2[i][j-1][k-1] = temp*(2*U2_[i-2][j][k]-13*U2_[i-1][j][k]+47*U2_[i][j][k]+27*U2_[i+1][j][k]-3*U2_[i+2][j][k]);
				ML3[i][j-1][k-1] = temp*(2*U3_[i-2][j][k]-13*U3_[i-1][j][k]+47*U3_[i][j][k]+27*U3_[i+1][j][k]-3*U3_[i+2][j][k]);
				ML4[i][j-1][k-1] = temp*(2*U4_[i-2][j][k]-13*U4_[i-1][j][k]+47*U4_[i][j][k]+27*U4_[i+1][j][k]-3*U4_[i+2][j][k]);
				ML5[i][j-1][k-1] = temp*(2*U5_[i-2][j][k]-13*U5_[i-1][j][k]+47*U5_[i][j][k]+27*U5_[i+1][j][k]-3*U5_[i+2][j][k]);

				temp = 1./(-3/J[i-2][j][k]+27/J[i-1][j][k]+47/J[i][j][k]-13/J[i+1][j][k]+2/J[i+2][j][k]);
						
				MR1[i-1][j-1][k-1] = temp*(-3*U1_[i-2][j][k]+27*U1_[i-1][j][k]+47*U1_[i][j][k]-13*U1_[i+1][j][k]+2*U1_[i+2][j][k]);
				MR2[i-1][j-1][k-1] = temp*(-3*U2_[i-2][j][k]+27*U2_[i-1][j][k]+47*U2_[i][j][k]-13*U2_[i+1][j][k]+2*U2_[i+2][j][k]);
				MR3[i-1][j-1][k-1] = temp*(-3*U3_[i-2][j][k]+27*U3_[i-1][j][k]+47*U3_[i][j][k]-13*U3_[i+1][j][k]+2*U3_[i+2][j][k]);
				MR4[i-1][j-1][k-1] = temp*(-3*U4_[i-2][j][k]+27*U4_[i-1][j][k]+47*U4_[i][j][k]-13*U4_[i+1][j][k]+2*U4_[i+2][j][k]);
				MR5[i-1][j-1][k-1] = temp*(-3*U5_[i-2][j][k]+27*U5_[i-1][j][k]+47*U5_[i][j][k]-13*U5_[i+1][j][k]+2*U5_[i+2][j][k]);
					
				
					
				//if (k == 24 & j == 24)
					//printf("%d\t%f\t%f\t%f\t%f\n",myid,ML1[i][j-1][k-1],MR1[i-1][j-1][k-1],ML2[i][j-1][k-1],MR2[i-1][j-1][k-1]);
					
			}
		}
	}

	
	if (myid == 0) {
		//#pragma omp parallel for private(k,k_)
		for ( j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {
		
				i = 2;
				
				temp = 1./(1/J[i][j][k]+1/J[i+1][j][k]);
			
				ML1[i][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i+1][j][k]);
				ML2[i][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i+1][j][k]);
				ML3[i][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i+1][j][k]);
				ML4[i][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i+1][j][k]);
				ML5[i][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i+1][j][k]);

				
// ====================================================================================//
		
				// i = 3; 
				
				// temp = 1./(1/J[i][j][k]+1/J[i+1][j][k]);
			
				// ML1[i][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i+1][j][k]);
				// ML2[i][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i+1][j][k]);
				// ML3[i][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i+1][j][k]);
				// ML4[i][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i+1][j][k]);
				// ML5[i][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i+1][j][k]);
				
				
				// MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i+1][j][k]);
				// MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i+1][j][k]);
				// MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i+1][j][k]);
				// MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i+1][j][k]);
				// MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i+1][j][k]);
				
// ====================================================================================//
						
				i = 3; 
				
				temp = 1./(-1/J[i-1][j][k]+5/J[i][j][k]+2/J[i+1][j][k]);

				ML1[i][j-1][k-1] = temp*(-U1_[i-1][j][k]+5*U1_[i][j][k]+2*U1_[i+1][j][k]);
				ML2[i][j-1][k-1] = temp*(-U2_[i-1][j][k]+5*U2_[i][j][k]+2*U2_[i+1][j][k]);
				ML3[i][j-1][k-1] = temp*(-U3_[i-1][j][k]+5*U3_[i][j][k]+2*U3_[i+1][j][k]);
				ML4[i][j-1][k-1] = temp*(-U4_[i-1][j][k]+5*U4_[i][j][k]+2*U4_[i+1][j][k]);
				ML5[i][j-1][k-1] = temp*(-U5_[i-1][j][k]+5*U5_[i][j][k]+2*U5_[i+1][j][k]);
				
				temp = 1./(2/J[i-1][j][k]+5/J[i][j][k]-1/J[i+1][j][k]);
						
			    
				MR1[i-1][j-1][k-1] = temp*(2*U1_[i-1][j][k]+5*U1_[i][j][k]-U1_[i+1][j][k]);
				MR2[i-1][j-1][k-1] = temp*(2*U2_[i-1][j][k]+5*U2_[i][j][k]-U2_[i+1][j][k]);
				MR3[i-1][j-1][k-1] = temp*(2*U3_[i-1][j][k]+5*U3_[i][j][k]-U3_[i+1][j][k]);
				MR4[i-1][j-1][k-1] = temp*(2*U4_[i-1][j][k]+5*U4_[i][j][k]-U4_[i+1][j][k]);
				MR5[i-1][j-1][k-1] = temp*(2*U5_[i-1][j][k]+5*U5_[i][j][k]-U5_[i+1][j][k]);
				
				}
			}
		}
		
		
	
	// if (3+gstart[myid] <= X_out-nx_outlet+3 && 3+gend0[myid] >= X_out-nx_outlet+3 ) {
			
		// for (j = ny-ny_abs+1; j <= nyy; j++) {
			// for (k = 2; k < nzz; k++) {

				// ii = X_out-nx_outlet+2-gstart[myid];
			
				// temp = 1./(1/J[ii][j][k]+1/J[ii+1][j][k]);
			
				// ML1[ii][j-1][k-1] = temp*(U1_[ii][j][k]+U1_[ii+1][j][k]);
				// ML2[ii][j-1][k-1] = temp*(U2_[ii][j][k]+U2_[ii+1][j][k]);
				// ML3[ii][j-1][k-1] = temp*(U3_[ii][j][k]+U3_[ii+1][j][k]);
				// ML4[ii][j-1][k-1] = temp*(U4_[ii][j][k]+U4_[ii+1][j][k]);
				// ML5[ii][j-1][k-1] = temp*(U5_[ii][j][k]+U5_[ii+1][j][k]);


// // ====================================================================================//
		
				
				// ii = X_out-nx_outlet+3-gstart[myid];
				
				
				// temp = 1./(1/J[ii][j][k]+1/J[ii+1][j][k]);
			
				// ML1[ii][j-1][k-1] = temp*(U1_[ii][j][k]+U1_[ii+1][j][k]);
				// ML2[ii][j-1][k-1] = temp*(U2_[ii][j][k]+U2_[ii+1][j][k]);
				// ML3[ii][j-1][k-1] = temp*(U3_[ii][j][k]+U3_[ii+1][j][k]);
				// ML4[ii][j-1][k-1] = temp*(U4_[ii][j][k]+U4_[ii+1][j][k]);
				// ML5[ii][j-1][k-1] = temp*(U5_[ii][j][k]+U5_[ii+1][j][k]);
				
				// MR1[ii-1][j-1][k-1] = temp*(U1_[ii][j][k]+U1_[ii+1][j][k]);
				// MR2[ii-1][j-1][k-1] = temp*(U2_[ii][j][k]+U2_[ii+1][j][k]);
				// MR3[ii-1][j-1][k-1] = temp*(U3_[ii][j][k]+U3_[ii+1][j][k]);
				// MR4[ii-1][j-1][k-1] = temp*(U4_[ii][j][k]+U4_[ii+1][j][k]);
				// MR5[ii-1][j-1][k-1] = temp*(U5_[ii][j][k]+U5_[ii+1][j][k]);
				

// // ====================================================================================//
						
				// ii = X_out-nx_outlet+4-gstart[myid];
				
				
				// temp = 1./(-1/J[ii-1][j][k]+5/J[ii][j][k]+2/J[ii+1][j][k]);

				// ML1[ii][j-1][k-1] = temp*(-U1_[ii-1][j][k]+5*U1_[ii][j][k]+2*U1_[ii+1][j][k]);
				// ML2[ii][j-1][k-1] = temp*(-U2_[ii-1][j][k]+5*U2_[ii][j][k]+2*U2_[ii+1][j][k]);
				// ML3[ii][j-1][k-1] = temp*(-U3_[ii-1][j][k]+5*U3_[ii][j][k]+2*U3_[ii+1][j][k]);
				// ML4[ii][j-1][k-1] = temp*(-U4_[ii-1][j][k]+5*U4_[ii][j][k]+2*U4_[ii+1][j][k]);
				// ML5[ii][j-1][k-1] = temp*(-U5_[ii-1][j][k]+5*U5_[ii][j][k]+2*U5_[ii+1][j][k]);
				
				// temp = 1./(2/J[ii-1][j][k]+5/J[ii][j][k]-1/J[ii+1][j][k]);
							
				// MR1[ii-1][j-1][k-1] = temp*(2*U1_[ii-1][j][k]+5*U1_[ii][j][k]-U1_[ii+1][j][k]);
				// MR2[ii-1][j-1][k-1] = temp*(2*U2_[ii-1][j][k]+5*U2_[ii][j][k]-U2_[ii+1][j][k]);
				// MR3[ii-1][j-1][k-1] = temp*(2*U3_[ii-1][j][k]+5*U3_[ii][j][k]-U3_[ii+1][j][k]);
				// MR4[ii-1][j-1][k-1] = temp*(2*U4_[ii-1][j][k]+5*U4_[ii][j][k]-U4_[ii+1][j][k]);
				// MR5[ii-1][j-1][k-1] = temp*(2*U5_[ii-1][j][k]+5*U5_[ii][j][k]-U5_[ii+1][j][k]);
					
					
			// }
		// }
		
	// }
		
		
		
		
		
	// if (3+gstart[myid] <= nx_inlet+2 && 3+gend0[myid] >= nx_inlet+2 ) {
			
		// for (j = ny-ny_abs+1; j <= nyy; j++) {
			// for (k = 2; k < nzz; k++) {
			
			
				// ii = nx_inlet+2-gstart[myid]-1;
				
				// temp = 1./(-1/J[ii-1][j][k]+5/J[ii][j][k]+2/J[ii+1][j][k]);

				
				// ML1[ii][j-1][k-1] = temp*(-U1_[ii-1][j][k]+5*U1_[ii][j][k]+2*U1_[ii+1][j][k]);
				// ML2[ii][j-1][k-1] = temp*(-U2_[ii-1][j][k]+5*U2_[ii][j][k]+2*U2_[ii+1][j][k]);
				// ML3[ii][j-1][k-1] = temp*(-U3_[ii-1][j][k]+5*U3_[ii][j][k]+2*U3_[ii+1][j][k]);
				// ML4[ii][j-1][k-1] = temp*(-U4_[ii-1][j][k]+5*U4_[ii][j][k]+2*U4_[ii+1][j][k]);
				// ML5[ii][j-1][k-1] = temp*(-U5_[ii-1][j][k]+5*U5_[ii][j][k]+2*U5_[ii+1][j][k]);
							
				// temp = 1./(2/J[ii-1][j][k]+5/J[ii][j][k]-1/J[ii+1][j][k]);
							
				// MR1[ii-1][j-1][k-1] = temp*(2*U1_[ii-1][j][k]+5*U1_[ii][j][k]-U1_[ii+1][j][k]);
				// MR2[ii-1][j-1][k-1] = temp*(2*U2_[ii-1][j][k]+5*U2_[ii][j][k]-U2_[ii+1][j][k]);
				// MR3[ii-1][j-1][k-1] = temp*(2*U3_[ii-1][j][k]+5*U3_[ii][j][k]-U3_[ii+1][j][k]);
				// MR4[ii-1][j-1][k-1] = temp*(2*U4_[ii-1][j][k]+5*U4_[ii][j][k]-U4_[ii+1][j][k]);
				// MR5[ii-1][j-1][k-1] = temp*(2*U5_[ii-1][j][k]+5*U5_[ii][j][k]-U5_[ii+1][j][k]);

				
				// //printf("%f\t%f\t%f\n",MR3[ii][j-1][k-1],U3_[ii][j][k]*J[ii][j][k],U3_[ii+1][j][k]*J[ii+1][j][k]);


// // ====================================================================================//
								
				// ii = nx_inlet+2-gstart[myid];
				
				// temp = 1./(1/J[ii][j][k]+1/J[ii-1][j][k]);
					
				// ML1[ii][j-1][k-1] = temp*(U1_[ii][j][k]+U1_[ii-1][j][k]);
				// ML2[ii][j-1][k-1] = temp*(U2_[ii][j][k]+U2_[ii-1][j][k]);
				// ML3[ii][j-1][k-1] = temp*(U3_[ii][j][k]+U3_[ii-1][j][k]);
				// ML4[ii][j-1][k-1] = temp*(U4_[ii][j][k]+U4_[ii-1][j][k]);
				// ML5[ii][j-1][k-1] = temp*(U5_[ii][j][k]+U5_[ii-1][j][k]);		
				
				
				// temp = 1./(1/J[ii][j][k]+1/J[ii-1][j][k]);
					
				// MR1[ii-1][j-1][k-1] = temp*(U1_[ii][j][k]+U1_[ii-1][j][k]);
				// MR2[ii-1][j-1][k-1] = temp*(U2_[ii][j][k]+U2_[ii-1][j][k]);
				// MR3[ii-1][j-1][k-1] = temp*(U3_[ii][j][k]+U3_[ii-1][j][k]);
				// MR4[ii-1][j-1][k-1] = temp*(U4_[ii][j][k]+U4_[ii-1][j][k]);
				// MR5[ii-1][j-1][k-1] = temp*(U5_[ii][j][k]+U5_[ii-1][j][k]);		
						


// // ====================================================================================//
										
						
				// ii = nx_inlet+2-gstart[myid];
				
				// temp = 1./(1/J[ii][j][k]+1/J[ii+1][j][k]);
					
				// MR1[ii][j-1][k-1] = temp*(U1_[ii][j][k]+U1_[ii+1][j][k]);
				// MR2[ii][j-1][k-1] = temp*(U2_[ii][j][k]+U2_[ii+1][j][k]);
				// MR3[ii][j-1][k-1] = temp*(U3_[ii][j][k]+U3_[ii+1][j][k]);
				// MR4[ii][j-1][k-1] = temp*(U4_[ii][j][k]+U4_[ii+1][j][k]);
				// MR5[ii][j-1][k-1] = temp*(U5_[ii][j][k]+U5_[ii+1][j][k]);
				
				
				
			// }
		// }
		
	// }
	
	
		
		
		
		

	if (myid == nproc-1) {
	
		for ( j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {
			
				
				i = gend[myid];
							
				temp = 1./(-1/J[i-1][j][k]+5/J[i][j][k]+2/J[i+1][j][k]);

				
				ML1[i][j-1][k-1] = temp*(-U1_[i-1][j][k]+5*U1_[i][j][k]+2*U1_[i+1][j][k]);
				ML2[i][j-1][k-1] = temp*(-U2_[i-1][j][k]+5*U2_[i][j][k]+2*U2_[i+1][j][k]);
				ML3[i][j-1][k-1] = temp*(-U3_[i-1][j][k]+5*U3_[i][j][k]+2*U3_[i+1][j][k]);
				ML4[i][j-1][k-1] = temp*(-U4_[i-1][j][k]+5*U4_[i][j][k]+2*U4_[i+1][j][k]);
				ML5[i][j-1][k-1] = temp*(-U5_[i-1][j][k]+5*U5_[i][j][k]+2*U5_[i+1][j][k]);
							
				temp = 1./(2/J[i-1][j][k]+5/J[i][j][k]-1/J[i+1][j][k]);
							
				MR1[i-1][j-1][k-1] = temp*(2*U1_[i-1][j][k]+5*U1_[i][j][k]-U1_[i+1][j][k]);
				MR2[i-1][j-1][k-1] = temp*(2*U2_[i-1][j][k]+5*U2_[i][j][k]-U2_[i+1][j][k]);
				MR3[i-1][j-1][k-1] = temp*(2*U3_[i-1][j][k]+5*U3_[i][j][k]-U3_[i+1][j][k]);
				MR4[i-1][j-1][k-1] = temp*(2*U4_[i-1][j][k]+5*U4_[i][j][k]-U4_[i+1][j][k]);
				MR5[i-1][j-1][k-1] = temp*(2*U5_[i-1][j][k]+5*U5_[i][j][k]-U5_[i+1][j][k]);



// ====================================================================================//
				
				
				// i = gend[myid];
				
				// temp = 1./(1/J[i][j][k]+1/J[i-1][j][k]);
				
				// ML1[i][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i-1][j][k]);
				// ML2[i][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i-1][j][k]);
				// ML3[i][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i-1][j][k]);
				// ML4[i][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i-1][j][k]);
				// ML5[i][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i-1][j][k]);
				
				
				// temp = 1./(1/J[i][j][k]+1/J[i-1][j][k]);
				
				// MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i-1][j][k]);
				// MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i-1][j][k]);
				// MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i-1][j][k]);
				// MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i-1][j][k]);
				// MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i-1][j][k]);
				
				
				
				

// ====================================================================================//
				
				
				i = gend[myid];
				
				temp = 1./(1/J[i][j][k]+1/J[i+1][j][k]);
					
				MR1[i][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i+1][j][k]);
				MR2[i][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i+1][j][k]);
				MR3[i][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i+1][j][k]);
				MR4[i][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i+1][j][k]);
				MR5[i][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i+1][j][k]);
				
				}
			}
		}
// ============================================================================ //
/**** MUSCL 5th-order-end****/


		
//// ============================================ ////
		istart = 2;							      ////
//// ============================================ ////
		iend = gend[myid];			     		  ////
//// ============================================ ////

	/*---X fluxes---*/
	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,XIx,XIy,XIz,\
	_rho,_u,_v,_w,iU,_V,_W,__U,_VV,iP,_T,iC,_H,\
	rho_,u_,v_,w_,U_,V_,W_,U__,VV_,P_,T_,C_,H_,\
	rho,u,v,w,U,V,W,_U_,VV,H,C,P,T,\
	dU1,dU2,dU3,dU4,dU5,\
	beta,S,_S_,\
	temp,temp1,temp2,temp3,temp4,temp5,tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2,\
	d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35,\
	d41,d42,d43,d44,d45,d51,d52,d53,d54,d55,\
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	k\
	)

		for (j = 1; j < ny; j++) {
			for (k = 1; k < nz; k++) {

				xix=xidx_u[i][j][k];
				xiy=xidy_u[i][j][k];
				xiz=xidz_u[i][j][k];
				etx=etdx_u[i][j][k];
				ety=etdy_u[i][j][k];
				etz=etdz_u[i][j][k];          
				ztx=ztdx_u[i][j][k];
				zty=ztdy_u[i][j][k];
				ztz=ztdz_u[i][j][k];
				XIx=xix/(sqrt(xix*xix+xiy*xiy+xiz*xiz));
				XIy=xiy/(sqrt(xix*xix+xiy*xiy+xiz*xiz));
				XIz=xiz/(sqrt(xix*xix+xiy*xiy+xiz*xiz));

				/* lefr parameter */
				_rho = ML1[i][j][k];
				_u = ML2[i][j][k]/_rho;
				_v = ML3[i][j][k]/_rho;
				_w = ML4[i][j][k]/_rho;
				

				iU = xix*_u+xiy*_v+xiz*_w;
				_V = etx*_u+ety*_v+etz*_w;
				_W = ztx*_u+zty*_v+ztz*_w;

				__U = XIx*_u+XIy*_v+XIz*_w;

				_VV = _u*_u+_v*_v+_w*_w;
				iP = (ML5[i][j][k]-0.5*_rho*_VV)*(K-1);

				_T = iP/_rho;
				iC = K*iP/_rho; /**** iC = sqrt(K*iP/_rho); ****/
				_H = 0.5*_VV+iC/(K-1);

				/* right parameter */
				rho_ = MR1[i][j][k];
				u_ = MR2[i][j][k]/rho_;
				v_ = MR3[i][j][k]/rho_;
				w_ = MR4[i][j][k]/rho_;
				
				U_ = xix*u_+xiy*v_+xiz*w_;
				V_ = etx*u_+ety*v_+etz*w_;
				W_ = ztx*u_+zty*v_+ztz*w_;

				U__ = XIx*u_+XIy*v_+XIz*w_;  

				VV_ = u_*u_+v_*v_+w_*w_;
				P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);

				T_ = P_/rho_;
				C_ = K*P_/rho_; /**** C_ = sqrt(K*P_/rho_); ****/
				H_ = 0.5*VV_+C_/(K-1);
				
				/* Low Mach correction */
				beta = sqrt(max(VV_/C_,_VV/iC));
				
				_u = 0.5*(_u+u_)+0.5*beta*(_u-u_);
				_v = 0.5*(_v+v_)+0.5*beta*(_v-v_);
				_w = 0.5*(_w+w_)+0.5*beta*(_w-w_);
				iU = 0.5*(iU+U_)+0.5*beta*(iU-U_);
				__U = 0.5*(__U+U__)+0.5*beta*(__U-U__);
				
				u_ = 0.5*(u_+_u)+0.5*beta*(u_-_u);
				v_ = 0.5*(v_+_v)+0.5*beta*(v_-_v);
				w_ = 0.5*(w_+_w)+0.5*beta*(w_-_w);
				U_ = 0.5*(U_+iU)+0.5*beta*(U_-iU);
				U__ = 0.5*(U__+__U)+0.5*beta*(U__-__U);
				
				/* average flux varibale */
				rho = sqrt(_rho*rho_);
				u = (_u*sqrt(_rho)+u_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				v = (_v*sqrt(_rho)+v_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				w = (_w*sqrt(_rho)+w_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				
				U = (iU*sqrt(_rho)+U_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				_U_ = (__U*sqrt(_rho)+U__*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				
				H = (_H*sqrt(_rho)+H_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				
				VV = u*u+v*v+w*w;
				
				C = (H-0.5*VV)*(K-1);
				P = rho*C/K;
				T = P/rho;
				
				/* jump dU */
				dU1 = rho_-_rho;
				dU2 =rho_*u_-_rho*_u;
				dU3 = rho_*v_-_rho*_v;
				dU4 = rho_*w_-_rho*_w;
				dU5 = (P_/(K-1)+0.5*rho_*(u_*u_+v_*v_+w_*w_))-(iP/(K-1)+0.5*_rho*(_u*_u+_v*_v+_w*_w));
								
				thetap = min(VV/C,1.0);
				Unp = 0.5*(1.0+thetap)*_U_;
				Cnp = 0.5*sqrt(4.0*C*thetap+(1.0-thetap)*(1.0-thetap)*_U_*_U_);
				
				theta = max(VV/C,e);
				Cnt = sqrt(max(VV/C,e))*sqrt(C);
				
				deltaU = (Cnp-0.5*(1.0-thetap)*_U_*Unp/Cnt-thetap*fabs(_U_))*(P_-iP)/(rho*theta*C)+Unp*(U__-__U)/Cnt;
				deltaP = Unp*(P_-iP)/Cnt+(Cnp-fabs(_U_)+0.5*(1.0-thetap)*_U_*Unp/Cnt)*rho*(U__-__U);
				
				/* artificial viscosity */
				Fav1 = fabs(_U_)*dU1+deltaU*rho;
				Fav2 = fabs(_U_)*dU2+deltaU*rho*u+deltaP*XIx;
				Fav3 = fabs(_U_)*dU3+deltaU*rho*v+deltaP*XIy;
				Fav4 = fabs(_U_)*dU4+deltaU*rho*w+deltaP*XIz;
				Fav5 = fabs(_U_)*dU5+deltaU*rho*H;
				
				/* inviscid fluxes */
				inFx1[i][j][k] = 0.5*((_rho*iU+rho_*U_)-Ep*Fav1*(sqrt(xix*xix+xiy*xiy+xiz*xiz)))/J_u[i][j][k];
				inFx2[i][j][k] = 0.5*((_rho*_u*iU+rho_*u_*U_+xix*(iP+P_))-Ep*Fav2*(sqrt(xix*xix+xiy*xiy+xiz*xiz)))/J_u[i][j][k];
				inFx3[i][j][k] = 0.5*((_rho*_v*iU+rho_*v_*U_+xiy*(iP+P_))-Ep*Fav3*(sqrt(xix*xix+xiy*xiy+xiz*xiz)))/J_u[i][j][k];
				inFx4[i][j][k] = 0.5*((_rho*_w*iU+rho_*w_*U_+xiz*(iP+P_))-Ep*Fav4*(sqrt(xix*xix+xiy*xiy+xiz*xiz)))/J_u[i][j][k];
				inFx5[i][j][k] = 0.5*((iU*(3.5*iP+0.5*_rho*_VV)+U_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5*(sqrt(xix*xix+xiy*xiy+xiz*xiz)))/J_u[i][j][k];

				
				// if(myid == 0) {
				
					// inFx1[i-1][j][k] = 0;
					// inFx2[i-1][j][k] = 0.5*(xix*(iP+P_))/J_u[i][j][k];
					// inFx3[i-1][j][k] = 0;
					// inFx4[i-1][j][k] = 0;
					// inFx5[i-1][j][k] = 0;
				
				// }
				
				// if(myid == np-1) {
				
					// inFx1[i+1][j][k] = 0;
					// inFx2[i+1][j][k] = 0.5*(xix*(iP+P_))/J_u[i][j][k];
					// inFx3[i+1][j][k] = 0;
					// inFx4[i+1][j][k] = 0;
					// inFx5[i+1][j][k] = 0;
				
				// }

			}
		}
	}
	

}