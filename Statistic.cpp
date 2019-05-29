




#include <stdio.h>
#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <math.h>


#include "Resolution.h"

extern int X_np;

void Statistic
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

double e,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

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

double (*EpY)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
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
	double u,v,w;
	double temp,temp1,temp2,temp3,temp4,temp5;
	double beta,S,_S_;

	double m11,m12,m13,m14,m15,m22,m23,m24,m25,m33,m34,m35,m44,m45,m55;
	double thedac1,thedac2,thedac3,thedad1,thedad2,thedad3;

	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
	double XIx,XIy,XIz,ETx,ETy,ETz,ZTx,ZTy,ZTz;
	double _rho,_u,_v,_w,iU,_V,_W,__U,__V,__W,_VV,iP,_T,iC,_H;
	double rho_,u_,v_,w_,U_,V_,W_,U__,V__,W__,VV_,P_,T_,C_,H_;
	double _U_,_V_,_W_;
	double dU1,dU2,dU3,dU4,dU5;
	double tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2;
	double d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35;
	double d41,d42,d43,d44,d45,d51,d52,d53,d54,d55, Fav1,Fav2,Fav3,Fav4,Fav5;


// ============================================================================================================= //
	
	static double 
		Um[X_m][Y_m],Vm[X_m][Y_m],Wm[X_m][Y_m],Tm[X_m][Y_m],
		UUm[X_m][Y_m],VVm[X_m][Y_m],WWm[X_m][Y_m],TTm[X_m][Y_m],UVm[X_m][Y_m],
		UTm[X_m][Y_m],VTm[X_m][Y_m],UIm[X_m][Y_m],TIm[X_m][Y_m],
		ULm[X_m][Y_m],URm[X_m][Y_m],
		VLm[X_m][Y_m],VRm[X_m][Y_m],
		rhoRm[X_m][Y_m],rhoLm[X_m][Y_m],
		RUVm[X_m][Y_m],SXYm[X_m][Y_m],Roe_Dm[X_m][Y_m],dRLm[X_m][Y_m];

	// ============================ //

	
	int x_np = gcount[myid]+6; 
	

/**** mean profile and turbulenct intensties ****/




// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //

	for (i = istart; i <= iend; i++) {
		for (j = 2; j <= ny; j++) {
				for (k = 2; k <= nz; k++) {

					rho = U1_[i][j][k];
					U = U2_[i][j][k]/rho;
					V = U3_[i][j][k]/rho;
					W = U4_[i][j][k]/rho;     
					VV = U*U+V*V+W*W;
					P = (U5_[i][j][k]-0.5*rho*VV)*(K-1)*J[i][j][k];
					T = P/(rho*R)/J[i][j][k];

					Um[i][j] = Um[i][j]+U;
					Vm[i][j] = Vm[i][j]+V;
					Wm[i][j] = Wm[i][j]+W;
					Tm[i][j] = Tm[i][j]+T;

					
					UUm[i][j] = UUm[i][j]+U*U;
					VVm[i][j] = VVm[i][j]+V*V;
					WWm[i][j] = WWm[i][j]+W*W;
					TTm[i][j] = TTm[i][j]+T*T;

					UVm[i][j] = UVm[i][j]+U*V;

					UTm[i][j] = UTm[i][j]+U*T;
					VTm[i][j] = VTm[i][j]+V*T;
					if (k == 150) {
					UIm[i][j] = U;
					TIm[i][j] = T;
					}
				}
			}

		}



// ======================================================== //
	if ( (step%statistic_step) == 0) {					    //   
// ======================================================== //
		
		double inv = 1./statistic_step/(nz-1);


// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //

		for (i = istart; i <= iend; i++) {
			for (j = 2; j <= ny; j++) {

				Um[i][j] = Um[i][j]*inv;
				Vm[i][j] = Vm[i][j]*inv;
				Wm[i][j] = Wm[i][j]*inv;
				Tm[i][j] = Tm[i][j]*inv;

				UUm[i][j] = UUm[i][j]*inv;
				VVm[i][j] = VVm[i][j]*inv;
				WWm[i][j] = WWm[i][j]*inv;
				TTm[i][j] = TTm[i][j]*inv;

				UVm[i][j] = UVm[i][j]*inv;
				UTm[i][j] = UTm[i][j]*inv;
				VTm[i][j] = VTm[i][j]*inv;

			}
		}


// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //


	for (i = istart; i <= iend; i++) {

		ii = i+istart1;

		if ((ii-2-nx_inlet)*deltaXI <= obs1 && (ii-1-nx_inlet)*deltaXI >= obs1) {

			FILE *fptr;
			sprintf(LESdata,"st_obs1_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Tm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTm[i][j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UIm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TIm[i][j]); }
			
			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs2 && (ii-1-nx_inlet)*deltaXI >= obs2) {

			FILE *fptr;
			sprintf(LESdata,"st_obs2_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Tm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTm[i][j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UIm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TIm[i][j]); }

			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs3 && (ii-1-nx_inlet)*deltaXI >= obs3) {

			FILE *fptr;
			sprintf(LESdata,"st_obs3_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Tm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTm[i][j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UIm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TIm[i][j]); }

			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs3 && (ii-1)*deltaXI >= obs3) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs4 && (ii-1-nx_inlet)*deltaXI >= obs4) {

			FILE *fptr;
			sprintf(LESdata,"st_obs4_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Tm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTm[i][j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UIm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TIm[i][j]); }

			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs4 && (ii-1)*deltaXI >= obs4) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs5 && (ii-1-nx_inlet)*deltaXI >= obs5) {

			FILE *fptr;
			sprintf(LESdata,"st_obs5_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Tm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTm[i][j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UIm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TIm[i][j]); }

			fclose(fptr);

		}    // ---- if ((i+1)*deltaXI <= obs5 && (i+2)*deltaXI >= obs5) ---- //
		
		


// ========================= //
		}					 //
// ========================= //


/**** mean temperature output end ****/



		
	int nx_out = X_out;
	int ny_out = Y_out;

	double (*Tm_out)[Y_out] = new double[X_out+1][Y_out];
	double (*Um_out)[Y_out] = new double[X_out+1][Y_out];

	ii = 2;


//// ========================= ////
	  istart = gstart[myid];   ////	
//// ========================= ////
	  iend = gend0[myid];      ////
//// ========================= ////


	  for (i = istart; i <= iend; i++) {

		  ii = ii+1;

		  for (j = 0; j < ny-1; j++) { 

			  Tm_out[i][j] =Tm[ii][j+2]; 
			  Um_out[i][j] =Um[ii][j+2]; 

		  }
	  }


	  MPI_Comm comm;
	  comm=MPI_COMM_WORLD;
	  MPI_Status istat[8];

	  if (myid > 0) {

		  istart=gstart[myid];
		  icount=gcount[myid]*Y_out;
		  idest=0;

		  itag = 210;
		  MPI_Send((void *)&Tm_out[istart][0], icount, MPI_DOUBLE, idest, itag, comm);
		  itag = 220;
		  MPI_Send((void *)&Um_out[istart][0], icount, MPI_DOUBLE, idest, itag, comm);

	  }

	  else {

		  for ( isrc=1; isrc < nproc; isrc++ ) {

			  istart=gstart[isrc];
			  icount=gcount[isrc]*Y_out;

			  itag = 210;
			  MPI_Recv((void *)&Tm_out[istart][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
			  itag = 220;
			  MPI_Recv((void *)&Um_out[istart][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
			  
		  }

	  }


	  if (myid == 0) {

			char LESdata[100];
			FILE *fptr;
			sprintf(LESdata,"Tave""%0.5d"".bin",step);
			fptr = fopen(LESdata,"wb");

			if (fptr != NULL)
			{

				fwrite(Tm_out,sizeof(double),X_out*Y_out,fptr);
				
				fclose(fptr);

			}

			else printf("File opening Failure\n");
			
			
			sprintf(LESdata,"Uave""%0.5d"".bin",step);
			fptr = fopen(LESdata,"wb");

			if (fptr != NULL)
			{

				fwrite(Um_out,sizeof(double),X_out*Y_out,fptr);
				
				fclose(fptr);

			}

			else printf("File opening Failure\n");

		}    // ---- if (myid == 0) ---- //


	  

	// =================== //
		istart = 3;        //
		iend = gend[myid]; //
	// =================== //

		for (i = istart; i <= iend; i++) {
			for (j = 2; j <= ny; j++) {

				Um[i][j] = 0;
				Vm[i][j] = 0;
				Wm[i][j] = 0;
				Tm[i][j] = 0;

				UUm[i][j] = 0;
				VVm[i][j] = 0;
				WWm[i][j] = 0;
				TTm[i][j] = 0;

				UVm[i][j] = 0;
				UTm[i][j] = 0;
				VTm[i][j] = 0;

			}
		}


		



	delete [] Tm_out;
	delete [] Um_out;



// ======================================================== //
	}  												        //
// ======================================================== //

/**** mean profile and turbulenct intensties - end ****/







}
