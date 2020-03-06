




#include <stdio.h>
#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <math.h>


#include "Resolution.h"

extern int X_np;

void Statistic_TK
(
// ============================================================================ //
int myid,
int step,
int iteration,
int statistic_step,

double obs1,
double obs2,
double obs3,
double obs4,
double obs5,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*xidy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*MR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ML1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

	
#include "ijk.h"
#include "Viscous_terms.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	istart1 = gstart[myid];

	char LESdata[100];
	
	double rho,U,V,W,VV,P,C,T,h,H;

	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
	// double _rho,_u,_v,_w,_U,_V,_W,__U,__V,__W,_VV,_P,_T,_C,_H;
	double rho_,u_,v_,w_,U_,V_,W_,U__,V__,W__,VV_,P_,T_,C_,H_;
	double _U_,_V_,_W_;
	double dU1,dU2,dU3,dU4,dU5;

	double 	irho,iU,iV,iW,iP,iT,
			jrho,jU,jV,jW,jP,jT,
			krho,kU,kV,kW,kP,kT,

			rhoi,Ui,Vi,Wi,Pi,Ti,
			rhoj,Uj,Vj,Wj,Pj,Tj,
			rhok,Uk,Vk,Wk,Pk,Tk,

			ijrho,ijU,ijV,ijW,ijP,ijT,
			jkrho,jkU,jkV,jkW,jkP,jkT,
			ikrho,ikU,ikV,ikW,ikP,ikT,

			rhoij,Uij,Vij,Wij,Pij,Tij,
			rhojk,Ujk,Vjk,Wjk,Pjk,Tjk,
			rhoik,Uik,Vik,Wik,Pik,Tik,

			irhoj,iUj,iVj,iWj,iPj,iTj,
			irhok,iUk,iVk,iWk,iPk,iTk,

			jrhoi,jUi,jVi,jWi,jPi,jTi,
			jrhok,jUk,jVk,jWk,jPk,jTk,

			krhoj,kUj,kVj,kWj,kPj,kTj,
			krhoi,kUi,kVi,kWi,kPi,kTi,

			du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,
			du2_dx,du2_dy,du2_dz,dv2_dx,dv2_dy,dv2_dz,dw2_dx,dw2_dy,dw2_dz,
			duv_dx,dvw_dx,duw_dx,duv_dy,dvw_dy,duw_dy,duv_dz,dvw_dz,duw_dz,
			dT_dx,dT_dy,dT_dz,

			mu_E,mu_T, Pr_E;

	double a1,a2,a3,a4,a5,a6,a7,
		   b1,b2,b3,b4,b5,b6,b7,
		   c1,c2,c3,c4,c5,c6,c7,
		   d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,
		   e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,
		   f1,f2,f3,f4,f5,f6,f7,f8,f9,f10;

	double Ux,Uy,Uz,
		   Vx,Vy,Vz,
		   Wx,Wy,Wz,
		   Tx,Ty,Tz;

	

	double invXI = 1./(deltaXI);
	double invET = 1./(deltaET);
	double invZT = 1./(deltaZT);

	double inv4XI = 1./(4*deltaXI);
	double inv4ET = 1./(4*deltaET);
	double inv4ZT = 1./(4*deltaZT);
	

// ============================================================================================================= //
	
	static double 
		UXm[X_m][Y_m],VXm[X_m][Y_m],TXm[X_m][Y_m],
		UUXm[X_m][Y_m],VVXm[X_m][Y_m],TTXm[X_m][Y_m],  
		UYm[X_m][Y_m],VYm[X_m][Y_m],TYm[X_m][Y_m],
		UUYm[X_m][Y_m],VVYm[X_m][Y_m],TTYm[X_m][Y_m];

	// ============================ //

	
	int x_np = gcount[myid]+6; 
	

/**** mean profile and turbulenct intensties ****/




//// ============================================ ////
		istart = 3;		             			  ////	
//// ============================================ ////
		iend = gend[myid]+1;					  ////
//// ============================================ ////
		
		for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,\
	a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6,c7,\
	d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,\
	f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,\
	rho,U,V,W,VV,P,T,\
	irho,iU,iV,iW,iP,iT,\
	jrho,jU,jV,jW,jP,jT,\
	krho,kU,kV,kW,kP,kT,\
	rhoi,Ui,Vi,Wi,Pi,Ti,\
	rhoj,Uj,Vj,Wj,Pj,Tj,\
	rhok,Uk,Vk,Wk,Pk,Tk,\
	ijrho,ijU,ijV,ijW,ijP,ijT,\
	jkrho,jkU,jkV,jkW,jkP,jkT,\
	ikrho,ikU,ikV,ikW,ikP,ikT,\
	rhoij,Uij,Vij,Wij,Pij,Tij,\
	rhojk,Ujk,Vjk,Wjk,Pjk,Tjk,\
	rhoik,Uik,Vik,Wik,Pik,Tik,\
	irhoj,iUj,iVj,iWj,iPj,iTj,\
	irhok,iUk,iVk,iWk,iPk,iTk,\
	jrhoi,jUi,jVi,jWi,jPi,jTi,\
	jrhok,jUk,jVk,jWk,jPk,jTk,\
	krhoj,kUj,kVj,kWj,kPj,kTj,\
	krhoi,kUi,kVi,kWi,kPi,kTi,\
	Ux,Uy,Uz,Tx,\
	Vx,Vy,Vz,Ty,\
	Wx,Wy,Wz,Tz,\
	du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,\
	du2_dx,du2_dy,du2_dz,dv2_dx,dv2_dy,dv2_dz,dw2_dx,dw2_dy,dw2_dz,\
	duv_dx,dvw_dx,duw_dx,duv_dy,dvw_dy,duw_dy,duv_dz,dvw_dz,duw_dz,\
	dT_dx,dT_dy,dT_dz,mu_E,Pr_E,\
	_k,__k,k\
	)
	
			for (j = 2; j <= nyy;j++) {
				for (__k = 0,_k = 1, k = 2; k <= nzz; __k++,_k++, k++) {


					// ==== 0 0 0 ==== //
				rho = U1_[i][j][k];
				U = U2_[i][j][k]/rho;
				V = U3_[i][j][k]/rho;
				W = U4_[i][j][k]/rho;     
				VV = U*U+V*V+W*W;
				P = (U5_[i][j][k]-0.5*rho*VV)*(K-1)*J[i][j][k];
				T = P/(rho*R)/J[i][j][k];


				// ==== -1 1 0 ==== //
				irhoj = U1_[i-1][j+1][k];
				iUj = U2_[i-1][j+1][k]/irhoj;
				iVj = U3_[i-1][j+1][k]/irhoj;
				iWj = U4_[i-1][j+1][k]/irhoj;
				iPj = (U5_[i-1][j+1][k]-0.5*irhoj*(iUj*iUj+iVj*iVj+iWj*iWj))*(K-1)*J[i-1][j+1][k];
				iTj = iPj/(irhoj*R)/J[i-1][j+1][k];

				// ==== -1 0 1 ==== //
				irhok = U1_[i-1][j][k+1];
				iUk = U2_[i-1][j][k+1]/irhok;
				iVk = U3_[i-1][j][k+1]/irhok;
				iWk = U4_[i-1][j][k+1]/irhok;
				iPk = (U5_[i-1][j][k+1]-0.5*irhok*(iUk*iUk+iVk*iVk+iWk*iWk))*(K-1)*J[i-1][j][k+1];
				iTk = iPk/(irhok*R)/J[i-1][j][k+1];

				// ==== 1 -1 0 ==== //
				jrhoi = U1_[i+1][j-1][k];
				jUi = U2_[i+1][j-1][k]/jrhoi;
				jVi = U3_[i+1][j-1][k]/jrhoi;
				jWi = U4_[i+1][j-1][k]/jrhoi;
				jPi = (U5_[i+1][j-1][k]-0.5*jrhoi*(jUi*jUi+jVi*jVi+jWi*jWi))*(K-1)*J[i+1][j-1][k];
				jTi = jPi/(jrhoi*R)/J[i+1][j-1][k];

				// ==== 1 0 -1 ==== //
				krhoi = U1_[i+1][j][k-1];
				kUi = U2_[i+1][j][k-1]/krhoi;
				kVi = U3_[i+1][j][k-1]/krhoi;
				kWi = U4_[i+1][j][k-1]/krhoi;
				kPi = (U5_[i+1][j][k-1]-0.5*krhoi*(kUi*kUi+kVi*kVi+kWi*kWi))*(K-1)*J[i+1][j][k-1];
				kTi = kPi/(krhoi*R)/J[i+1][j][k-1];


				// ==== 0 -1 1 ==== //
				jrhok = U1_[i][j-1][k+1];
				jUk = U2_[i][j-1][k+1]/jrhok;
				jVk = U3_[i][j-1][k+1]/jrhok;
				jWk = U4_[i][j-1][k+1]/jrhok;
				jPk = (U5_[i][j-1][k+1]-0.5*jrhok*(jUk*jUk+jVk*jVk+jWk*jWk))*(K-1)*J[i][j-1][k+1];
				jTk = jPk/(jrhok*R)/J[i][j-1][k+1];


				// ==== 0 1 -1 ==== //
				krhoj = U1_[i][j+1][k-1];
				kUj = U2_[i][j+1][k-1]/krhoj;
				kVj = U3_[i][j+1][k-1]/krhoj;
				kWj = U4_[i][j+1][k-1]/krhoj;
				kPj = (U5_[i][j+1][k-1]-0.5*krhoj*(kUj*kUj+kVj*kVj+kWj*kWj))*(K-1)*J[i][j+1][k-1];
				kTj = kPj/(krhoj*R)/J[i][j+1][k-1];

				// ==== -1 0 0 ==== //
				irho = U1_[i-1][j][k];
				iU = U2_[i-1][j][k]/irho;
				iV = U3_[i-1][j][k]/irho;
				iW = U4_[i-1][j][k]/irho;
				iP = (U5_[i-1][j][k]-0.5*irho*(iU*iU+iV*iV+iW*iW))*(K-1)*J[i-1][j][k];
				iT = iP/(irho*R)/J[i-1][j][k];

				// ==== 0 -1 0 ==== //
				jrho = U1_[i][j-1][k];
				jU = U2_[i][j-1][k]/jrho;
				jV = U3_[i][j-1][k]/jrho;
				jW = U4_[i][j-1][k]/jrho;
				jP = (U5_[i][j-1][k]-0.5*jrho*(jU*jU+jV*jV+jW*jW))*(K-1)*J[i][j-1][k];
				jT = jP/(jrho*R)/J[i][j-1][k];

				// ==== 0 0 -1 ==== //
				krho = U1_[i][j][k-1];
				kU = U2_[i][j][k-1]/krho;
				kV = U3_[i][j][k-1]/krho;
				kW = U4_[i][j][k-1]/krho;
				kP = (U5_[i][j][k-1]-0.5*krho*(kU*kU+kV*kV+kW*kW))*(K-1)*J[i][j][k-1];
				kT = kP/(krho*R)/J[i][j][k-1];

				// ==== -1 -1 0 ==== //
				ijrho = U1_[i-1][j-1][k];
				ijU = U2_[i-1][j-1][k]/ijrho;
				ijV = U3_[i-1][j-1][k]/ijrho;
				ijW = U4_[i-1][j-1][k]/ijrho;
				ijP = (U5_[i-1][j-1][k]-0.5*ijrho*(ijU*ijU+ijV*ijV+ijW*ijW))*(K-1)*J[i-1][j-1][k];
				ijT = ijP/(ijrho*R)/J[i-1][j-1][k];


				// ==== 0 -1 -1 ==== //
				jkrho = U1_[i][j-1][k-1];
				jkU = U2_[i][j-1][k-1]/jkrho;
				jkV = U3_[i][j-1][k-1]/jkrho;
				jkW = U4_[i][j-1][k-1]/jkrho;
				jkP = (U5_[i][j-1][k-1]-0.5*jkrho*(jkU*jkU+jkV*jkV+jkW*jkW))*(K-1)*J[i][j-1][k-1];
				jkT = jkP/(jkrho*R)/J[i][j-1][k-1];


				// ==== -1 0 -1 ==== //
				ikrho = U1_[i-1][j][k-1];
				ikU = U2_[i-1][j][k-1]/ikrho;
				ikV = U3_[i-1][j][k-1]/ikrho;
				ikW = U4_[i-1][j][k-1]/ikrho;
				ikP = (U5_[i-1][j][k-1]-0.5*ikrho*(ikU*ikU+ikV*ikV+ikW*ikW))*(K-1)*J[i-1][j][k-1];
				ikT = ikP/(ikrho*R)/J[i-1][j][k-1];
				



				// ==== 1 0 0 ==== //
				rhoi = U1_[i+1][j][k];
				Ui = U2_[i+1][j][k]/rhoi;
				Vi = U3_[i+1][j][k]/rhoi;
				Wi = U4_[i+1][j][k]/rhoi;
				Pi = (U5_[i+1][j][k]-0.5*rhoi*(Ui*Ui+Vi*Vi+Wi*Wi))*(K-1)*J[i+1][j][k];
				Ti = Pi/(rhoi*R)/J[i+1][j][k];

				// ==== 0 1 0 ==== //
				rhoj = U1_[i][j+1][k];
				Uj = U2_[i][j+1][k]/rhoj;
				Vj = U3_[i][j+1][k]/rhoj;
				Wj = U4_[i][j+1][k]/rhoj;
				Pj = (U5_[i][j+1][k]-0.5*rhoj*(Uj*Uj+Vj*Vj+Wj*Wj))*(K-1)*J[i][j+1][k];
				Tj = Pj/(rhoj*R)/J[i][j+1][k];

				// ==== 0 0 1  ==== //
				rhok = U1_[i][j][k+1];
				Uk = U2_[i][j][k+1]/rhok;
				Vk = U3_[i][j][k+1]/rhok;
				Wk = U4_[i][j][k+1]/rhok;
				Pk = (U5_[i][j][k+1]-0.5*rhok*(Uk*Uk+Vk*Vk+Wk*Wk))*(K-1)*J[i][j][k+1];
				Tk = Pk/(rhok*R)/J[i][j][k+1];
        
				xiy=xidy_u[i-1][j-1][k-1];
				ety=etdy_u[i-1][j-1][k-1];
				zty=ztdy_u[i-1][j-1][k-1];

				
				/* derivatives of velocity */
				/* X-direction */
				du_dx = (U-iU)*invXI;
				dv_dx = (V-iV)*invXI;
				dT_dx = (T-iT)*invXI;

				/* Y-direction */
				du_dy = (Uj+iUj-jU-ijU)*inv4ET;
				dv_dy = (Vj+iVj-jV-ijV)*inv4ET;
				dT_dy = (Tj+iTj-jT-ijT)*inv4ET;
        
        /* to computational domain */
        du_dy = du_dy*xiy+du_dy*ety+du_dy*zty;
				dv_dy = dv_dy*xiy+dv_dy*ety+dv_dy*zty;
				dT_dy = dT_dy*xiy+dT_dy*ety+dT_dy*zty;
        
        MR1[i][j][k] = du_dx;
        MR2[i][j][k] = dv_dx;
        MR3[i][j][k] = dT_dx;
        
        ML1[i][j][k] = du_dy;
        ML2[i][j][k] = dv_dy;
        ML3[i][j][k] = dT_dy;
        

				
					}
				}
			}

#pragma omp barrier
    
    
    
    
    
//// =============================== ////
		istart = 3;		             		   ////	
//// =============================== ////
		iend = gend[myid]+1;					   ////
//// =============================== ////
      
	for (i = istart; i <= iend; i++) {
		for (j = 2; j <= ny; j++) {
				for (k = 2; k <= nz; k++) {

					UXm[i][j] = UXm[i][j]+MR1[i][j][k];
					VXm[i][j] = VXm[i][j]+MR2[i][j][k];
					TXm[i][j] = TXm[i][j]+MR3[i][j][k];

					UUXm[i][j] = UUXm[i][j]+MR1[i][j][k]*MR1[i][j][k];
					VVXm[i][j] = VVXm[i][j]+MR2[i][j][k]*MR2[i][j][k];
					TTXm[i][j] = TTXm[i][j]+MR3[i][j][k]*MR3[i][j][k];

          
					UYm[i][j] = UYm[i][j]+ML1[i][j][k];
					VYm[i][j] = VYm[i][j]+ML2[i][j][k];
					TYm[i][j] = TYm[i][j]+ML3[i][j][k];

					UUYm[i][j] = UUYm[i][j]+ML1[i][j][k]*ML1[i][j][k];
					VVYm[i][j] = VVYm[i][j]+ML2[i][j][k]*ML2[i][j][k];
					TTYm[i][j] = TTYm[i][j]+ML3[i][j][k]*ML3[i][j][k];

          
					}
				}
			}




// ======================================================== //
	if ( (step%statistic_step) == 0) {					    //   
// ======================================================== //
		
		double inv = 1./statistic_step/(nz-1);


// =================== //
	istart = 3;        //
	iend = gend[myid]+1; //
// =================== //

		for (i = istart; i <= iend; i++) {
			for (j = 2; j <= ny; j++) {

				UXm[i][j] = UXm[i][j]*inv;
				VXm[i][j] = VXm[i][j]*inv;
				TXm[i][j] = TXm[i][j]*inv;

				UUXm[i][j] = UUXm[i][j]*inv;
				VVXm[i][j] = VVXm[i][j]*inv;
				TTXm[i][j] = TTXm[i][j]*inv;

				UYm[i][j] = UYm[i][j]*inv;
				VYm[i][j] = VYm[i][j]*inv;
				TYm[i][j] = TYm[i][j]*inv;

				UUYm[i][j] = UUYm[i][j]*inv;
				VVYm[i][j] = VVYm[i][j]*inv;
				TTYm[i][j] = TTYm[i][j]*inv;

			}
		}


// =================== //
	istart = 3;        //
	iend = gend[myid]+1; //
// =================== //


	for (i = istart; i <= iend; i++) {

		ii = i+istart1;

		if ((ii-2-nx_inlet)*deltaXI <= obs1 && (ii-1-nx_inlet)*deltaXI >= obs1) {

			FILE *fptr;
			sprintf(LESdata,"TK_obs1_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TXm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTXm[i][j]); }
      
      
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TYm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTYm[i][j]); }

			
			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs2 && (ii-1-nx_inlet)*deltaXI >= obs2) {

			FILE *fptr;
			sprintf(LESdata,"TK_obs2_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TXm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTXm[i][j]); }
      
      
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TYm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTYm[i][j]); }

			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs3 && (ii-1-nx_inlet)*deltaXI >= obs3) {

			FILE *fptr;
			sprintf(LESdata,"TK_obs3_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TXm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTXm[i][j]); }
      
      
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TYm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTYm[i][j]); }

			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs3 && (ii-1)*deltaXI >= obs3) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs4 && (ii-1-nx_inlet)*deltaXI >= obs4) {

			FILE *fptr;
			sprintf(LESdata,"TK_obs4_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TXm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTXm[i][j]); }
      
      
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TYm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTYm[i][j]); }

			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs4 && (ii-1)*deltaXI >= obs4) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs5 && (ii-1-nx_inlet)*deltaXI >= obs5) {

			FILE *fptr;
			sprintf(LESdata,"TK_obs5_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TXm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVXm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTXm[i][j]); }
      
      
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TYm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVYm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTYm[i][j]); }

			fclose(fptr);

		}    // ---- if ((i+1)*deltaXI <= obs5 && (i+2)*deltaXI >= obs5) ---- //
		
		}


// ========================= //
		}					 //
// ========================= //






}
