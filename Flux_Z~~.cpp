



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "Resolution.h"

extern int X_np;

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

void Flux_Z
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
double (*J_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*inFz1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*EpZ)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
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
	double thetap,Wnp,Cnp,theta,Cnt,deltaU,deltaP,beta;
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
		for ( j = 2; j <= ny; j++) {
			for (k = 3; k <= nz-1; k++) {
			
				temp = 1./(2/J[i][j][k-2]-13/J[i][j][k-1]+47/J[i][j][k]+27/J[i][j][k+1]-3/J[i][j][k+2]);

				ML1[i-1][j-1][k] = temp*(2*U1_[i][j][k-2]-13*U1_[i][j][k-1]+47*U1_[i][j][k]+27*U1_[i][j][k+1]-3*U1_[i][j][k+2]);
				ML2[i-1][j-1][k] = temp*(2*U2_[i][j][k-2]-13*U2_[i][j][k-1]+47*U2_[i][j][k]+27*U2_[i][j][k+1]-3*U2_[i][j][k+2]);
				ML3[i-1][j-1][k] = temp*(2*U3_[i][j][k-2]-13*U3_[i][j][k-1]+47*U3_[i][j][k]+27*U3_[i][j][k+1]-3*U3_[i][j][k+2]);
				ML4[i-1][j-1][k] = temp*(2*U4_[i][j][k-2]-13*U4_[i][j][k-1]+47*U4_[i][j][k]+27*U4_[i][j][k+1]-3*U4_[i][j][k+2]);
				ML5[i-1][j-1][k] = temp*(2*U5_[i][j][k-2]-13*U5_[i][j][k-1]+47*U5_[i][j][k]+27*U5_[i][j][k+1]-3*U5_[i][j][k+2]);

				temp = 1./(-3/J[i][j][k-2]+27/J[i][j][k-1]+47/J[i][j][k]-13/J[i][j][k+1]+2/J[i][j][k+2]);
						
				MR1[i-1][j-1][k-1] = temp*(-3*U1_[i][j][k-2]+27*U1_[i][j][k-1]+47*U1_[i][j][k]-13*U1_[i][j][k+1]+2*U1_[i][j][k+2]);
				MR2[i-1][j-1][k-1] = temp*(-3*U2_[i][j][k-2]+27*U2_[i][j][k-1]+47*U2_[i][j][k]-13*U2_[i][j][k+1]+2*U2_[i][j][k+2]);
				MR3[i-1][j-1][k-1] = temp*(-3*U3_[i][j][k-2]+27*U3_[i][j][k-1]+47*U3_[i][j][k]-13*U3_[i][j][k+1]+2*U3_[i][j][k+2]);
				MR4[i-1][j-1][k-1] = temp*(-3*U4_[i][j][k-2]+27*U4_[i][j][k-1]+47*U4_[i][j][k]-13*U4_[i][j][k+1]+2*U4_[i][j][k+2]);
				MR5[i-1][j-1][k-1] = temp*(-3*U5_[i][j][k-2]+27*U5_[i][j][k-1]+47*U5_[i][j][k]-13*U5_[i][j][k+1]+2*U5_[i][j][k+2]);
				
				
				//if (myid == 0 & i == 3 & j == 24)
						//printf("%d\t%f\t%f\t%f\t%f\n",k,ML1[i-1][j-1][k],MR1[i-1][j-1][k-1],ML2[i-1][j-1][k],MR2[i-1][j-1][k-1]);
						
					
			}
		}
	}
	
	for (i = istart; i <= iend; i++) {
		for (j = 2; j <= ny; j++) {
		
			k = 1;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k+1]);
			
			ML1[i-1][j-1][k] = temp*(U1_[i][j][k]+U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(U2_[i][j][k]+U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(U3_[i][j][k]+U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(U4_[i][j][k]+U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(U5_[i][j][k]+U5_[i][j][k+1]);
			
			
				//if (myid == 0 & i == 3 & j == 24)
						//printf("%d\t%f\t%f\t%f\t%f\n",k,ML1[i-1][j-1][k],MR1[i-1][j-1][k],ML2[i-1][j-1][k],MR2[i-1][j-1][k]);
						
			

// ====================================================================================//
		
			/*
			k = 2;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k+1]);
					
			ML1[i-1][j-1][k] = temp*(U1_[i][j][k]+U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(U2_[i][j][k]+U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(U3_[i][j][k]+U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(U4_[i][j][k]+U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(U5_[i][j][k]+U5_[i][j][k+1]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k+1]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j][k+1]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j][k+1]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j][k+1]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j][k+1]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j][k+1]);
			*/
			
			k = 2; 
			
			temp = 1./(-1/J[i][j][k-1]+5/J[i][j][k]+2/J[i][j][k+1]);

			ML1[i-1][j-1][k] = temp*(-U1_[i][j][k-1]+5*U1_[i][j][k]+2*U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(-U2_[i][j][k-1]+5*U2_[i][j][k]+2*U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(-U3_[i][j][k-1]+5*U3_[i][j][k]+2*U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(-U4_[i][j][k-1]+5*U4_[i][j][k]+2*U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(-U5_[i][j][k-1]+5*U5_[i][j][k]+2*U5_[i][j][k+1]);
						
			temp = 1./(2/J[i][j][k-1]+5/J[i][j][k]-1/J[i][j][k+1]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j][k-1]+5*U1_[i][j][k]-U1_[i][j][k+1]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j][k-1]+5*U2_[i][j][k]-U2_[i][j][k+1]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j][k-1]+5*U3_[i][j][k]-U3_[i][j][k+1]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j][k-1]+5*U4_[i][j][k]-U4_[i][j][k+1]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j][k-1]+5*U5_[i][j][k]-U5_[i][j][k+1]);
			
			
			
// ====================================================================================//		
			/*
			k = 3; 
			
			temp = 1./(-1/J[i][j][k-1]+5/J[i][j][k]+2/J[i][j][k+1]);

			ML1[i-1][j-1][k] = temp*(-U1_[i][j][k-1]+5*U1_[i][j][k]+2*U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(-U2_[i][j][k-1]+5*U2_[i][j][k]+2*U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(-U3_[i][j][k-1]+5*U3_[i][j][k]+2*U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(-U4_[i][j][k-1]+5*U4_[i][j][k]+2*U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(-U5_[i][j][k-1]+5*U5_[i][j][k]+2*U5_[i][j][k+1]);
						
			temp = 1./(2/J[i][j][k-1]+5/J[i][j][k]-1/J[i][j][k+1]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j][k-1]+5*U1_[i][j][k]-U1_[i][j][k+1]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j][k-1]+5*U2_[i][j][k]-U2_[i][j][k+1]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j][k-1]+5*U3_[i][j][k]-U3_[i][j][k+1]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j][k-1]+5*U4_[i][j][k]-U4_[i][j][k+1]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j][k-1]+5*U5_[i][j][k]-U5_[i][j][k+1]);
			
			
			k = nz-1;
			
			temp = 1./(-1/J[i][j][k-1]+5/J[i][j][k]+2/J[i][j][k+1]);

			ML1[i-1][j-1][k] = temp*(-U1_[i][j][k-1]+5*U1_[i][j][k]+2*U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(-U2_[i][j][k-1]+5*U2_[i][j][k]+2*U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(-U3_[i][j][k-1]+5*U3_[i][j][k]+2*U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(-U4_[i][j][k-1]+5*U4_[i][j][k]+2*U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(-U5_[i][j][k-1]+5*U5_[i][j][k]+2*U5_[i][j][k+1]);
						
			temp = 1./(2/J[i][j][k-1]+5/J[i][j][k]-1/J[i][j][k+1]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j][k-1]+5*U1_[i][j][k]-U1_[i][j][k+1]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j][k-1]+5*U2_[i][j][k]-U2_[i][j][k+1]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j][k-1]+5*U3_[i][j][k]-U3_[i][j][k+1]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j][k-1]+5*U4_[i][j][k]-U4_[i][j][k+1]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j][k-1]+5*U5_[i][j][k]-U5_[i][j][k+1]);
			*/
			
			
			
// ====================================================================================//


			/*
			k = nz;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k-1]);
					
			ML1[i-1][j-1][k] = temp*(U1_[i][j][k]+U1_[i][j][k-1]);
			ML2[i-1][j-1][k] = temp*(U2_[i][j][k]+U2_[i][j][k-1]);
			ML3[i-1][j-1][k] = temp*(U3_[i][j][k]+U3_[i][j][k-1]);
			ML4[i-1][j-1][k] = temp*(U4_[i][j][k]+U4_[i][j][k-1]);
			ML5[i-1][j-1][k] = temp*(U5_[i][j][k]+U5_[i][j][k-1]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k-1]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j][k-1]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j][k-1]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j][k-1]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j][k-1]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j][k-1]);
			*/
			
			k = nz;
			
			temp = 1./(-1/J[i][j][k-1]+5/J[i][j][k]+2/J[i][j][k+1]);

			ML1[i-1][j-1][k] = temp*(-U1_[i][j][k-1]+5*U1_[i][j][k]+2*U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(-U2_[i][j][k-1]+5*U2_[i][j][k]+2*U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(-U3_[i][j][k-1]+5*U3_[i][j][k]+2*U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(-U4_[i][j][k-1]+5*U4_[i][j][k]+2*U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(-U5_[i][j][k-1]+5*U5_[i][j][k]+2*U5_[i][j][k+1]);
						
			temp = 1./(2/J[i][j][k-1]+5/J[i][j][k]-1/J[i][j][k+1]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j][k-1]+5*U1_[i][j][k]-U1_[i][j][k+1]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j][k-1]+5*U2_[i][j][k]-U2_[i][j][k+1]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j][k-1]+5*U3_[i][j][k]-U3_[i][j][k+1]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j][k-1]+5*U4_[i][j][k]-U4_[i][j][k+1]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j][k-1]+5*U5_[i][j][k]-U5_[i][j][k+1]);
			
			
// ====================================================================================//

			k = nz;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k+1]);
				
			MR1[i-1][j-1][k] = temp*(U1_[i][j][k]+U1_[i][j][k+1]);
			MR2[i-1][j-1][k] = temp*(U2_[i][j][k]+U2_[i][j][k+1]);
			MR3[i-1][j-1][k] = temp*(U3_[i][j][k]+U3_[i][j][k+1]);
			MR4[i-1][j-1][k] = temp*(U4_[i][j][k]+U4_[i][j][k+1]);
			MR5[i-1][j-1][k] = temp*(U5_[i][j][k]+U5_[i][j][k+1]);
			
			
				//if (myid == 0 & i == 3 & j == 24)
						//printf("%d\t%f\t%f\t%f\t%f\n",k,ML1[i-1][j-1][k],MR1[i-1][j-1][k],ML2[i-1][j-1][k],MR2[i-1][j-1][k]);
						
			
				
		}
	}
	
	


	

//// ============================================ ////
			istart = 2;							  ////	
//// ============================================ ////
			iend = gend[myid]-1;	     		  ////
//// ============================================ ////

		/*---Z fluxes---*/
		for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,ZTx,ZTy,ZTz,\
	_rho,_u,_v,_w,iU,_V,_W,__W,_VV,iP,_T,iC,_H,\
	rho_,u_,v_,w_,U_,V_,W_,W__,VV_,P_,T_,C_,H_,\
	rho,u,v,w,U,V,W,_W_,VV,H,C,P,T,\
	dU1,dU2,dU3,dU4,dU5,\
	beta,S,_S_,\
	temp,temp1,temp2,temp3,temp4,temp5,tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2,\
	d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35,\
	d41,d42,d43,d44,d45,d51,d52,d53,d54,d55,\
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	k\
	)

			for (j = 1; j < ny; j++) {
				for (k = 1; k < nzz; k++) {

				xix=xidx_w[i][j][k];
				xiy=xidy_w[i][j][k];
				xiz=xidz_w[i][j][k];
				etx=etdx_w[i][j][k];
				ety=etdy_w[i][j][k];
				etz=etdz_w[i][j][k];          
				ztx=ztdx_w[i][j][k];
				zty=ztdy_w[i][j][k];
				ztz=ztdz_w[i][j][k];
				ZTx=ztx/(sqrt(ztx*ztx+zty*zty+ztz*ztz));
				ZTy=zty/(sqrt(ztx*ztx+zty*zty+ztz*ztz));
				ZTz=ztz/(sqrt(ztx*ztx+zty*zty+ztz*ztz));

				/* lefr parameter */
				_rho = ML1[i][j][k];
				_u = ML2[i][j][k]/_rho;
				_v = ML3[i][j][k]/_rho;
				_w = ML4[i][j][k]/_rho;

				iU = xix*_u+xiy*_v+xiz*_w;
				_V = etx*_u+ety*_v+etz*_w;
				_W = ztx*_u+zty*_v+ztz*_w;

				__W = ZTx*_u+ZTy*_v+ZTz*_w;

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

				W__ = ZTx*u_+ZTy*v_+ZTz*w_;       

				VV_ = u_*u_+v_*v_+w_*w_;
				P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);
				T_ = P_/rho_;
				C_ = K*P_/rho_;
				H_ = 0.5*VV_+C_/(K-1);
				
				/* Low Mach correction */
				beta = sqrt(max(VV_/C_,_VV/iC));
				
				_u = 0.5*(_u+u_)+0.5*beta*(_u-u_);
				_v = 0.5*(_v+v_)+0.5*beta*(_v-v_);
				_w = 0.5*(_w+w_)+0.5*beta*(_w-w_);
				_W = 0.5*(_W+W_)+0.5*beta*(_W-W_);
				__W = 0.5*(__W+W__)+0.5*beta*(__W-W__);
				
				u_ = 0.5*(u_+_u)+0.5*beta*(u_-_u);
				v_ = 0.5*(v_+_v)+0.5*beta*(v_-_v);
				w_ = 0.5*(w_+_w)+0.5*beta*(w_-_w);
				W_ = 0.5*(W_+_W)+0.5*beta*(W_-_W);
				W__ = 0.5*(W__+__W)+0.5*beta*(W__-__W);
				
				/* average flux varibale */
				
				rho = sqrt(_rho*rho_);
				u = (_u*sqrt(_rho)+u_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				v = (_v*sqrt(_rho)+v_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				w = (_w*sqrt(_rho)+w_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				
				W = (_W*sqrt(_rho)+W_*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				
				_W_ = (__W*sqrt(_rho)+W__*sqrt(rho_))/(sqrt(_rho)+sqrt(rho_));
				
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
				Wnp = 0.5*(1.0+thetap)*_W_;
				Cnp = 0.5*sqrt(4.0*C*thetap+(1.0-thetap)*(1.0-thetap)*_W_*_W_);
				
				theta = max(VV/C,e);
				Cnt = sqrt(max(VV/C,e))*sqrt(C);
				
				deltaU = (Cnp-0.5*(1.0-thetap)*_W_*Wnp/Cnt-thetap*fabs(_W_))*(P_-iP)/(rho*theta*C)+Wnp*(W__-__W)/Cnt;
				deltaP = Wnp*(P_-iP)/Cnt+(Cnp-fabs(_W_)+0.5*(1.0-thetap)*_W_*Wnp/Cnt)*rho*(W__-__W);
				
				/* artificial viscosity */
				Fav1 = fabs(_W_)*dU1+deltaU*rho;
				Fav2 = fabs(_W_)*dU2+deltaU*rho*u+deltaP*ZTx;
				Fav3 = fabs(_W_)*dU3+deltaU*rho*v+deltaP*ZTy;
				Fav4 = fabs(_W_)*dU4+deltaU*rho*w+deltaP*ZTz;
				Fav5 = fabs(_W_)*dU5+deltaU*rho*H;

				
				/* inviscid fluxes */
				inFz1[i][j][k] = 0.5*((_rho*_W+rho_*W_-Ep*Fav1*(sqrt(ztx*ztx+zty*zty+ztz*ztz))))/J_w[i][j][k];
				inFz2[i][j][k] = 0.5*((_rho*_u*_W+rho_*u_*W_+ztx*(iP+P_))-Ep*Fav2*(sqrt(ztx*ztx+zty*zty+ztz*ztz)))/J_w[i][j][k];
				inFz3[i][j][k] = 0.5*((_rho*_v*_W+rho_*v_*W_+zty*(iP+P_))-Ep*Fav3*(sqrt(ztx*ztx+zty*zty+ztz*ztz)))/J_w[i][j][k];
				inFz4[i][j][k] = 0.5*((_rho*_w*_W+rho_*w_*W_+ztz*(iP+P_))-Ep*Fav4*(sqrt(ztx*ztx+zty*zty+ztz*ztz)))/J_w[i][j][k];
				inFz5[i][j][k] = 0.5*((_W*(3.5*iP+0.5*_rho*_VV)+W_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5*(sqrt(ztx*ztx+zty*zty+ztz*ztz)))/J_w[i][j][k];

			}
		}
	}

}