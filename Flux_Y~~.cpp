



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "Resolution.h"

extern int X_np;

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

void Flux_Y
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
double (*J_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*inFy1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*EpY)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
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
	double thetap,Vnp,Cnp,theta,Cnt,deltaU,deltaP,beta;
	double Fav1,Fav2,Fav3,Fav4,Fav5;

/**** MUSCL 5th-order ****/
// ============================================================================ //
	
//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 3;						  ////	
//// ============================================ ////
		iend = gend[myid];						  ////
//// ============================================ ////
	
	

	for (i = istart; i <= iend; i++) {
#pragma omp parallel for private(k,temp)
		for ( j = 4; j <= ny-2; j++) {
			for (k = 2; k <= nz; k++) {
			
				temp = 1./(2/J[i][j-2][k]-13/J[i][j-1][k]+47/J[i][j][k]+27/J[i][j+1][k]-3/J[i][j+2][k]);

				ML1[i-1][j][k-1] = temp*(2*U1_[i][j-2][k]-13*U1_[i][j-1][k]+47*U1_[i][j][k]+27*U1_[i][j+1][k]-3*U1_[i][j+2][k]);
				ML2[i-1][j][k-1] = temp*(2*U2_[i][j-2][k]-13*U2_[i][j-1][k]+47*U2_[i][j][k]+27*U2_[i][j+1][k]-3*U2_[i][j+2][k]);
				ML3[i-1][j][k-1] = temp*(2*U3_[i][j-2][k]-13*U3_[i][j-1][k]+47*U3_[i][j][k]+27*U3_[i][j+1][k]-3*U3_[i][j+2][k]);
				ML4[i-1][j][k-1] = temp*(2*U4_[i][j-2][k]-13*U4_[i][j-1][k]+47*U4_[i][j][k]+27*U4_[i][j+1][k]-3*U4_[i][j+2][k]);
				ML5[i-1][j][k-1] = temp*(2*U5_[i][j-2][k]-13*U5_[i][j-1][k]+47*U5_[i][j][k]+27*U5_[i][j+1][k]-3*U5_[i][j+2][k]);

				temp = 1./(-3/J[i][j-2][k]+27/J[i][j-1][k]+47/J[i][j][k]-13/J[i][j+1][k]+2/J[i][j+2][k]);
						
				MR1[i-1][j-1][k-1] = temp*(-3*U1_[i][j-2][k]+27*U1_[i][j-1][k]+47*U1_[i][j][k]-13*U1_[i][j+1][k]+2*U1_[i][j+2][k]);
				MR2[i-1][j-1][k-1] = temp*(-3*U2_[i][j-2][k]+27*U2_[i][j-1][k]+47*U2_[i][j][k]-13*U2_[i][j+1][k]+2*U2_[i][j+2][k]);
				MR3[i-1][j-1][k-1] = temp*(-3*U3_[i][j-2][k]+27*U3_[i][j-1][k]+47*U3_[i][j][k]-13*U3_[i][j+1][k]+2*U3_[i][j+2][k]);
				MR4[i-1][j-1][k-1] = temp*(-3*U4_[i][j-2][k]+27*U4_[i][j-1][k]+47*U4_[i][j][k]-13*U4_[i][j+1][k]+2*U4_[i][j+2][k]);
				MR5[i-1][j-1][k-1] = temp*(-3*U5_[i][j-2][k]+27*U5_[i][j-1][k]+47*U5_[i][j][k]-13*U5_[i][j+1][k]+2*U5_[i][j+2][k]);
					
				
			}
		}
	}
	

	
	
	for (i = istart; i <= iend; i++) {
		for (k = 2; k <= nz; k++) {
		
			
		
			j = 1;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			

// ====================================================================================//
		
			
			j = 2;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
					
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
// ====================================================================================//		
			j = 3; 
			
			temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

			ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
						
			temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
			
// ====================================================================================//	
		
			j = ny-1;
			
			temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

			ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
						
			temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
			
// ====================================================================================//



			j = ny;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
					
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
			
			
// ====================================================================================//

			j = ny;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
				
			MR1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			MR2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			MR3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			MR4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			MR5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
				
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	for (i = istart; i <= iend; i++) {
	
		if (i+gstart[myid] > nx_inlet+2 && i+gstart[myid] < (X_out-nx_outlet)+3) {
		
			for (k = 2; k <= nz; k++) {
		
			
			j = Y_out-ny_abs+1;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			

// ====================================================================================//
		
			
			j = Y_out-ny_abs+2;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
					
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
// ====================================================================================//		
			j = Y_out-ny_abs+3;
			
			temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

			ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
						
			temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
			
// ====================================================================================//	
				
		}
		
		
		
		
		
			for (k = 2; k <= nz; k++) {
			
	// ------------------------------------------------------------------------------------ //			
	// ---------------------------------- adiabatic wall ---------------------------------- //

	
	// ====================================================================================//		
				j = Y_out-ny_abs; 
				
				temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

				ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
				ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
				ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
				ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
				ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
							
				temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
							
				MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
				MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
				MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
				MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
				MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
				
// ====================================================================================//	
							
				
				j = Y_out-ny_abs+1;
				
				temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
						
				ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
				ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
				ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
				ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
				ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
				
				
				temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
				
				MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
				MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
				MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
				MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
				MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
				
				
	// ====================================================================================//	
			
				j = Y_out-ny_abs+1;
				
				temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
				
				MR1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
				MR2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
				MR3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
				MR4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
				MR5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
				
				
			}
			
		}
	}
		
	
	
	



//// ============================================ ////
			istart = 2;							  ////	
//// ============================================ ////
			iend = gend[myid]-1;	     		  ////
//// ============================================ ////

	/*---Y fluxes---*/
	for (i = istart; i <= iend; i++) {

	#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,ETx,ETy,ETz,\
	_rho,_u,_v,_w,iU,_V,_W,__V,_VV,iP,_T,iC,_H,\
	rho_,u_,v_,w_,U_,V_,W_,V__,VV_,P_,T_,C_,H_,\
	rho,u,v,w,U,V,W,_V_,VV,H,C,P,T,\
	dU1,dU2,dU3,dU4,dU5,\
	beta,S,_S_,\
	temp,temp1,temp2,temp3,temp4,temp5,tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2,\
	d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35,\
	d41,d42,d43,d44,d45,d51,d52,d53,d54,d55,\
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	k\
	)

		for (j = 1; j < nyy; j++) {
			for (k = 1; k < nz; k++) {

				xix=xidx_v[i][j][k];
				xiy=xidy_v[i][j][k];
				xiz=xidz_v[i][j][k];
				etx=etdx_v[i][j][k];
				ety=etdy_v[i][j][k];
				etz=etdz_v[i][j][k];          
				ztx=ztdx_v[i][j][k];
				zty=ztdy_v[i][j][k];
				ztz=ztdz_v[i][j][k];
				ETx=etx/(sqrt(etx*etx+ety*ety+etz*etz));
				ETy=ety/(sqrt(etx*etx+ety*ety+etz*etz));
				ETz=etz/(sqrt(etx*etx+ety*ety+etz*etz));

				/* lefr parameter */
				_rho = ML1[i][j][k];
				_u = ML2[i][j][k]/_rho;
				_v = ML3[i][j][k]/_rho;
				_w = ML4[i][j][k]/_rho;

				iU = xix*_u+xiy*_v+xiz*_w;
				_V = etx*_u+ety*_v+etz*_w;

				_W = ztx*_u+zty*_v+ztz*_w;


				__V = ETx*_u+ETy*_v+ETz*_w;


				_VV = _u*_u+_v*_v+_w*_w;
				iP = (ML5[i][j][k]-0.5*_rho*_VV)*(K-1);
				_T = iP/_rho;
				iC = K*iP/_rho;
				_H = 0.5*_VV+iC/(K-1);

				/* right parameter */
				rho_ = MR1[i][j][k];
				u_ = MR2[i][j][k]/rho_;
				v_ = MR3[i][j][k]/rho_;
				w_ = MR4[i][j][k]/rho_;

				U_ = xix*u_+xiy*v_+xiz*w_;
				V_ = etx*u_+ety*v_+etz*w_;
				W_ = ztx*u_+zty*v_+ztz*w_;

				V__ = ETx*u_+ETy*v_+ETz*w_;

				VV_ = u_*u_+v_*v_+w_*w_;
				P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);
				T_ = P_/rho_;
				C_ =K*P_/rho_;
				H_ = 0.5*VV_+C_/(K-1);
				
				/* Low Mach correction */
				beta = sqrt(max(VV_/C_,_VV/iC));
				
				_u = 0.5*(_u+u_)+0.5*beta*(_u-u_);
				_v = 0.5*(_v+v_)+0.5*beta*(_v-v_);
				_w = 0.5*(_w+w_)+0.5*beta*(_w-w_);
				_V = 0.5*(_V+V_)+0.5*beta*(_V-V_);
				__V = 0.5*(__V+V__)+0.5*beta*(__V-V__);
				
				u_ = 0.5*(u_+_u)+0.5*beta*(u_-_u);
				v_ = 0.5*(v_+_v)+0.5*beta*(v_-_v);
				w_ = 0.5*(w_+_w)+0.5*beta*(w_-_w);
				V_ = 0.5*(V_+_V)+0.5*beta*(V_-_V);
				V__ = 0.5*(V__+__V)+0.5*beta*(V__-__V);
				
				/* average flux varibale */
				
				rho = sqrt(_rho*rho_);
				u = (_u*sqrt(_rho)+u_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				v = (_v*sqrt(_rho)+v_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				w = (_w*sqrt(_rho)+w_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				
				V = (_V*sqrt(_rho)+V_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				
				_V_ = (__V*sqrt(_rho)+V__*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				
				VV = u*u+v*v+w*w;
				H = (_H*sqrt(_rho)+H_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
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
				Vnp = 0.5*(1.0+thetap)*_V_;
				Cnp = 0.5*sqrt(4.0*C*thetap+(1.0-thetap)*(1.0-thetap)*_V_*_V_);
				
				theta = max(VV/C,e);
				Cnt = sqrt(max(VV/C,e))*sqrt(C);
				
				deltaU = (Cnp-0.5*(1.0-thetap)*_V_*Vnp/Cnt-thetap*fabs(_V_))*(P_-iP)/(rho*theta*C)+Vnp*(V__-__V)/Cnt;
				deltaP = Vnp*(P_-iP)/Cnt+(Cnp-fabs(_V_)+0.5*(1.0-thetap)*_V_*Vnp/Cnt)*rho*(V__-__V);
				
				/* artificial viscosity */
				Fav1 = fabs(_V_)*dU1+deltaU*rho;
				Fav2 = fabs(_V_)*dU2+deltaU*rho*u+deltaP*ETx;
				Fav3 = fabs(_V_)*dU3+deltaU*rho*v+deltaP*ETy;
				Fav4 = fabs(_V_)*dU4+deltaU*rho*w+deltaP*ETz;
				Fav5 = fabs(_V_)*dU5+deltaU*rho*H;
				
				/* inviscid fluxes */
				
				inFy1[i][j][k] = 0.5*((_rho*_V+rho_*V_)-Ep*Fav1*(sqrt(etx*etx+ety*ety+etz*etz)))/J_v[i][j][k];
				inFy2[i][j][k] = 0.5*((_rho*_u*_V+rho_*u_*V_+etx*(iP+P_))-Ep*Fav2*(sqrt(etx*etx+ety*ety+etz*etz)))/J_v[i][j][k];
				inFy3[i][j][k] = 0.5*((_rho*_v*_V+rho_*v_*V_+ety*(iP+P_))-Ep*Fav3*(sqrt(etx*etx+ety*ety+etz*etz)))/J_v[i][j][k];
				inFy4[i][j][k] = 0.5*((_rho*_w*_V+rho_*w_*V_+etz*(iP+P_))-Ep*Fav4*(sqrt(etx*etx+ety*ety+etz*etz)))/J_v[i][j][k];
				inFy5[i][j][k] = 0.5*((_V*(3.5*iP+0.5*_rho*_VV)+V_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5*(sqrt(etx*etx+ety*ety+etz*etz)))/J_v[i][j][k];
				

				if (j==1 | j==ny) {

					inFy1[i][j][k] = 0;
					inFy2[i][j][k] = 0;
					inFy3[i][j][k] = 0.5*(ety*(iP+P_)-Ep*Fav3)/J_v[i][j][k];
					inFy4[i][j][k] = 0;
					inFy5[i][j][k] = 0;

				}
				
				
			}
		}
	}
	

}