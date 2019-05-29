


#define nx_inlet  304
#define nx_outlet 432

#define ny_abs 0

#define X_out 4000  /**** X_out+nx_in+nx_out ****/  
#define Y_out 384   /**** Y_out+ny_abs ****/
#define Z_out 300     

#define X_m 4004   /**** X_out+4 ****/
#define Y_m 388     /**** Y_out+4 ****/
#define Z_m 304      /**** Z_out+4 ****/

#define nx X_out+1    /**** nx+1 ****/
#define ny Y_out+1    /**** ny+1 ****/
#define nz Z_out+1    /**** nz+1 ****/

#define deltaXI 0.00153    
#define deltaET 0.00130208                /**** 1/Y_out because already normalized ****/
#define deltaZT 0.001 

#define nxx nx+1
#define nyy ny+1
#define nzz nz+1

#define nxxx nxx+1
#define nyyy nyy+1
#define nzzz nzz+1 

