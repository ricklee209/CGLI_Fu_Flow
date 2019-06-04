



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
	double temp,temp1,temp2,temp3,temp4,temp5,temp6;
	double deltaU, deltaP, Cdiss;
	double beta,S,_S_,insqr;
    double theda_p, U_p, C_p;

	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
	double XIx,XIy,XIz,ETx,ETy,ETz,ZTx,ZTy,ZTz;
	double _rho,_u,_v,_w,iU,_V,_W,__U,__V,__W,_VV,iP,_T,iC,_H;
	double rho_,u_,v_,w_,U_,V_,W_,U__,V__,W__,VV_,P_,T_,C_,H_;
	double _U_,_V_,_W_;
	double dU1,dU2,dU3,dU4,dU5;
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
	
			
			
// ====================================================================================//


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
	beta,S,\
	temp,temp1,temp2,temp3,temp4,temp5,temp6,\
    deltaU, deltaP, Cdiss,insqr, \
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	k\
	)
    
    

			for (j = 1; j < ny; j++) {
				for (k = 1; k < nzz; k++) {

				ztx=ztdx_w[i][j][k];
				zty=ztdy_w[i][j][k];
				ztz=ztdz_w[i][j][k];

				insqr = 1.0/(sqrt(ztx*ztx+zty*zty+ztz*ztz));

				ZTx=ztx*insqr;
				ZTy=zty*insqr;
				ZTz=ztz*insqr;


				/* lefr parameter */
				_rho = ML1[i][j][k];
				_u = ML2[i][j][k]/_rho;
				_v = ML3[i][j][k]/_rho;
				_w = ML4[i][j][k]/_rho;

				_W = ztx*_u+zty*_v+ztz*_w;

				_VV = _u*_u+_v*_v+_w*_w;
				iP = (ML5[i][j][k]-0.5*_rho*_VV)*(K-1);
				_T = _P/_rho;
				iC = K*_P/_rho;
				_H = 0.5*_VV+iC/(K-1);

				/* right parameter */
				rho_ = MR1[i][j][k];
				u_ = MR2[i][j][k]/rho_;
				v_ = MR3[i][j][k]/rho_;
				w_ = MR4[i][j][k]/rho_;

				W_ = ztx*u_+zty*v_+ztz*w_;

				VV_ = u_*u_+v_*v_+w_*w_;
				P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);
				T_ = P_/rho_;
				C_ = K*P_/rho_;
				H_ = 0.5*VV_+C_/(K-1);
				
            
                /* flux varibale */
                temp5 = sqrt(_rho);
                temp6 = sqrt(rho_);
                
                temp4 = temp5+temp6;

                rho = sqrt(_rho*rho_);
                u = (temp5*_u+temp6*u_)/temp4;
                v = (temp5*_v+temp6*v_)/temp4;;
                w = (temp5*_w+temp6*w_)/temp4;

                W = (temp5*_W+temp6*W_)/temp4;

                VV = u*u+v*v+w*w;
                H = (temp5*_H+temp6*H_)/temp4;
                C = (H-0.5*VV)*(K-1); /**** C = sqrt((H-0.5*VV)*(K-1)); ****/
                P = rho*C/K;
                
                /* jump dU */
            
                dU1 = rho_-_rho;
                dU2 = rho_*u_-_rho*_u;
                dU3 = rho_*v_-_rho*_v;
                dU4 = rho_*w_-_rho*_w;
                dU5 = (P_/(K-1)+0.5*rho_*VV_)-(iP/(K-1)+0.5*_rho*_VV);

                
                beta = max(VV/C,e);    // ---- theda ---- //
                
                temp = 0.5*(1+beta)*W*insqr;    // ---- W' ---- //

                S = 0.5*sqrt(4*beta*C+W*insqr*W*insqr*(1-beta)*(1-beta));   // ---- C' ---- //

                theda_p = VV/C;

                U_p = 0.5*(1+theda_p)*W*insqr;
                C_p = 0.5*sqrt(4*C*theda_p+W*insqr*W*insqr*(1-theda_p)*(1-theda_p));

                temp1 = (P_-iP)/rho/beta/C;
                temp2 = U_p/S*(W_-_W)*insqr;
                temp3 = 0.5*(1-beta)*W*insqr*temp/S;
                temp4 = 0.5*(1-theda_p)*W*insqr*U_p/S;

                deltaU = (S-temp3-beta*fabs(W)*insqr)*temp1+temp2;

                deltaP = U_p/S*(P_-iP)+(C_p-fabs(W)*insqr+temp4)*rho*(W_-_W)*insqr;
                
                
                
                Fav1 = fabs(W*insqr)*dU1+deltaU*rho;
                Fav2 = fabs(W*insqr)*dU2+deltaU*rho*u+deltaP*XIx;
                Fav3 = fabs(W*insqr)*dU3+deltaU*rho*v+deltaP*XIy;
                Fav4 = fabs(W*insqr)*dU4+deltaU*rho*w+deltaP*XIz;
                Fav5 = fabs(W*insqr)*dU5+deltaU*rho*H+deltaP*W*insqr;




				/* inviscid fluxes */
				
                
				inFz1[i][j][k] = 0.5*((_rho*_W+rho_*W_)-Fav1/insqr)/J_w[i][j][k];
				inFz2[i][j][k] = 0.5*((_rho*_u*_W+rho_*u_*W_+ztx*(iP+P_))-Fav2/insqr)/J_w[i][j][k];
				inFz3[i][j][k] = 0.5*((_rho*_v*_W+rho_*v_*W_+zty*(iP+P_))-Fav3/insqr)/J_w[i][j][k];
				inFz4[i][j][k] = 0.5*((_rho*_w*_W+rho_*w_*W_+ztz*(iP+P_))-Fav4/insqr)/J_w[i][j][k];
				inFz5[i][j][k] = 0.5*((_W*(3.5*iP+0.5*_rho*_VV)+W_*(3.5*P_+0.5*rho_*VV_))-Fav5/insqr)/J_w[i][j][k];
				

			}
		}
	}

}