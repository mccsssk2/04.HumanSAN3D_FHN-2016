/*
SK June 25, 2015.
FK Mod Version for mouse SAN and atrial cell types.
17/3/2015. The cell types are ready.

This 2 cell system is to find out at what
frequency and strength the SAN can suppress
the paranodal, and what frequencies and 
strength of SAN the paranodal oscillator takes over.

The output will be analysed separately.
The output I am looking for is the input frequency (SAN) vs output frequency (paranodal).

Do simulations at:
a) many SAN frequencies
b) FFT of frequencies
b) GJC values from very low to very high. High is 0.1, which gives synchronisation in 4 or 5 oscillations.
*/

//------------------------------------------------------------------------------
// my standard header and #defines, as of 22 Dec. 2013
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ    */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ    */
#define RTOL  RCONST(1.0e-12)   	  /* scalar relative tolerance            */
#define ATOL  RCONST(1.0e-6)   		  /* scalar absolute tolerance components */
#define MAXSTEPS 5000
#define ZERO  RCONST(0.0)

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

#define doutput 0

/* Problem Constants */

#define num_beats   200
#define NEQ         3                	   /* number of equations  */
#define DELTAT       0.1
#define amp          0.5
#define duration     1.0
#define NUMPARAMS    13
#define NP	     	 1
#define pcl          1000.0

#define VOLTAGE 0

#define aLargeNumber 100000

typedef struct {
	realtype IV;             // stimulus
	realtype lap_val;        // the laplacian.
	realtype p[NUMPARAMS];
	realtype dvdt;
	realtype q; // the factor that is used by yout theory papers (Ostborn JTB 2001 and others) to control the firing rate of SAN.
	realtype dvdt_control;
// the firing follows the theory that peripheral cells produce a SMALLER reaction than the central. This peripheral and central nomenclature
    realtype delta_V; // difference between next cells voltage and my voltage.
// is according to the Manchester lot, I disagree with them - the SAN is not as simple as that.
} *UserData;

/* Functions Called by the Solver */
#include "f.c"
#include "randNum.c"

// these are error weights.
static int ewt(N_Vector y, N_Vector w, void *user_data);

int main(int argc, char *argv[]){

  N_Vector    y;                        // yS is the sensitivity solution with length y and breadth pbar
  void       *cvode_mem;                // memory location pointer.
  UserData    data;                     // a structure that will contain the modelling parameters. This is essentially the same 

/* standard declarations for CVODE only. */
  FILE       *output, *random;			// output file for time straces.
  realtype    t, tout;                  // time variable t, and tout is the output counter.
  int         iout, stimInt;         // iout is an int that tracks the time, stimInt is for multiple stimuli, plist is array for parameters.
  char       *str;                      // str is generally used in file i/o
  int         NOUT;                     // 
  double      sstime = 0.0;             // "time" that I can write to file.
  tout      = DELTAT;                   // initial tout value.
  iout      = 0;
  NOUT = (int)(((num_beats)*pcl)/DELTAT); // how often do you output.

// cell # 0 is SAN, cell # 1 is paranodal.
double usr_u0_loc[2], usr_u1_loc[2], usr_u2_loc[2], usr_dvdt[2];
double old_v0, old_v1;
double gjc, q_value, dvdt_control_value;

  booleantype err_con;
  err_con    = TRUE;

realtype dvdtmax0 = -10000000000.0;
realtype dvdtmax1 = -10000000000.0;

int file_counter;

/* Initialising the random number generator is important. */
random = fopen("/dev/urandom","rb");
fread(&somerandomnumber, sizeof(unsigned int), 1, random);
fclose(random);
init_genrand(somerandomnumber);

for(file_counter = 0; file_counter < 1000000; file_counter++){

err_con = TRUE;
dvdtmax0 = -1000000000.0    ;
dvdtmax1 = -1000000000.0    ;
sstime  = 0.0              ;
iout    = 0                ;
stimInt = (int)(pcl/DELTAT);
usr_u0_loc[0] = usr_u0_loc[1] = 0; usr_u1_loc[0] = usr_u1_loc[1] = 0; usr_u2_loc[0] = usr_u2_loc[1] = 0;

// genrand_real1()
gjc                  = 0.015      *genrand_real1();   // input parameter 1. range is 0 to 0.3.
q_value              = 10.0 + 2000*genrand_real1();   // input parameter 2. range is 10 to 2010
dvdt_control_value   = 1.0;                           // input parameter 3. do not do anything in this pass.

// cell 0.
double t00, v00, t01, v01, t02, v02; // first min., max, second min.
int mes0_done = -1;

// cell 1
double t10, v10, t11, v11, t12, v12; // first min., max, second min.
int mes1_done = -1;

double cl0[aLargeNumber] ,  cl1[aLargeNumber];
double std0, std1, mean0, mean1;
double max_cl0, min_cl0;
double max_cl1, min_cl1;
int n0, n1;
int number_of_cls0, number_of_cls1;

int count0, count1;
count0 = count1 = 0;

for(n0 = 0; n0 < aLargeNumber; n0++){
	cl0[n0] = cl1[n0] = -1000.0;
	mean0 = mean1 = 0;
	std0  = std1  = 0;
	number_of_cls0 = number_of_cls1 = 0;
	max_cl0 = min_cl0 = -1000.0;
	max_cl1 = min_cl1 = -1000.0;
}

//printf("%d\n",file_counter);

n0 = n1 = 0;

  for(iout = 0; iout <= NOUT; iout++) { // start of time loop.

// san cell.**************************************************************************************
  cvode_mem = NULL;
  data      = NULL;
  y         = NULL;
/* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ);
  data = (UserData) malloc(sizeof *data); // now it is created.
   Ith(y,1) = usr_u0_loc[0];
   Ith(y,2) = usr_u1_loc[0];  
   Ith(y,3) = usr_u2_loc[0];  
   old_v0 = usr_u0_loc[0];
  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
    // solver
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    CVodeInit(cvode_mem, f, 0.0, y);
    /* Error weights. This is used during run time. */
    CVodeWFtolerances(cvode_mem, ewt);
    CVodeSetUserData(cvode_mem, data); // you must tell cvode_mem about data.
    CVodeSStolerances(cvode_mem, RTOL, ATOL);
    CVodeSetSensErrCon(cvode_mem, err_con);              // always have err_con switched on.
    CVDense(cvode_mem, NEQ);
    CVodeSetMaxStep(cvode_mem, DELTAT);
    tout = DELTAT;
	/* set up time independent user data. */
	data->lap_val      = 0.0;
// 0.05 130.0 50.0 12.50 10.0 18.2 18.2 1.5 600.0 0.1 0.5 0.01 5.0
// 0.05 130.0 50.0 12.50 10.0 18.2 18.2 1.5 600.0 0.1 0.5 0.01 5.0 2
// SAN parameters.

data->dvdt_control = dvdt_control_value; 
data->q = q_value;

// SAN cell.
			data->p[0] /* tau_d   */  =  0.05 ; // upstroke
			data->p[1] /* tau_r   */  =  130.0;
			data->p[2] /* tau_si  */  =  50.0 ;
			data->p[3] /* tau_0   */  =  12.5 ;
			data->p[4] /* tau_vp  */  =  10.0 ;
			data->p[5] /* tau_vm1 */  =  18.2 ;
			data->p[6] /* tau_vm2 */  =  18.2 ;
			data->p[7] /* tau_wp  */  =  1.5  ;
			data->p[8] /* tau_wm  */  =  data->q; // 600.0; // pacing frequency.
			data->p[9] /* u_c     */  =  0.1  ;
			data->p[10] /* u_v     */ =  0.5  ;
			data->p[11] /* usi_c   */ =  0.01 ;
			data->p[12] /* k       */ =  5.0  ;
			data->IV = 0.0;
			data->delta_V = gjc*(old_v1 - old_v0);
		    CVode(cvode_mem, tout, y, &t, CV_NORMAL);
usr_dvdt[0] = data->dvdt;
   usr_u0_loc[0] = Ith(y,1);
   usr_u1_loc[0] = Ith(y,2);  
   usr_u2_loc[0] = Ith(y,3);  

// do your measurement of cell # 0 here.
// dvdtmax during this oscillation.
if(dvdtmax0 < usr_dvdt[0]) dvdtmax0 = usr_dvdt[0];
// minimum.
if( (usr_dvdt[0]>0.0) && (fabs(usr_dvdt[0])<0.001) && (old_v0 < usr_u0_loc[0]) && (mes0_done==-1) && (usr_u0_loc[0]<0.2)){

cl0[n0] = sstime - t00; n0++; /*********************/

t00 = sstime; v00 = usr_u0_loc[0];
mes0_done = 1;
}

// this is close to the minimum, but 1 time step before.
if( (usr_dvdt[0]<0.0) && (fabs(usr_dvdt[0])<0.001) && (old_v0 > usr_u0_loc[0]) && (mes0_done==1) && (usr_u0_loc[0]>0.2) ){

if(count0>10){

str = malloc(32*sizeof(char));
sprintf(str,"cell0data.%d",atoi(argv[1]));
output = fopen(str,"a+");
free(str);

//fprintf(output,"%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\n",file_counter, 0, gjc, q_value, data->dvdt_control, cl0[n0-1], dvdtmax0,n0);
fclose(output);
}

dvdtmax0 = -100000000.0;
mes0_done = -1;
count0++;
}

  N_VDestroy_Serial(y);
  CVodeFree(&cvode_mem);
  free(data);
//***************************************************************************************
// paranodal cell.***********************************************************************
  cvode_mem = NULL;
  data      = NULL;
  y         = NULL;
/* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ);
  data = (UserData) malloc(sizeof *data); // now it is created.
   Ith(y,1) = usr_u0_loc[1];
   Ith(y,2) = usr_u1_loc[1];  
   Ith(y,3) = usr_u2_loc[1];  
   old_v1 = usr_u0_loc[1];
  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
    // solver
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    CVodeInit(cvode_mem, f, 0.0, y);
    /* Error weights. This is used during run time. */
    CVodeWFtolerances(cvode_mem, ewt);
    CVodeSetUserData(cvode_mem, data); // you must tell cvode_mem about data.
    CVodeSStolerances(cvode_mem, RTOL, ATOL);
    CVodeSetSensErrCon(cvode_mem, err_con);              // always have err_con switched on.
    CVDense(cvode_mem, NEQ);
    CVodeSetMaxStep(cvode_mem, DELTAT);
    tout = DELTAT;
	/* set up time independent user data. */
	data->lap_val      = 0.0;
// 0.05 130.0 50.0 12.50 10.0 18.2 18.2 1.5 800.0 0.1 0.5 0.01 5.0
// paranodal parameters.

data->dvdt_control = -1000.0; 
data->q            = -1000.0;

			data->p[0] /* tau_d   */  =  0.175 ; // upstroke modulator.
			data->p[1] /* tau_r   */  =  130.0;
			data->p[2] /* tau_si  */  =  50.0 ;
			data->p[3] /* tau_0   */  =  12.5 ;
			data->p[4] /* tau_vp  */  =  10.0 ;
			data->p[5] /* tau_vm1 */  =  18.2 ;
			data->p[6] /* tau_vm2 */  =  18.2 ;
			data->p[7] /* tau_wp  */  =  1.5  ;
			data->p[8] /* tau_wm  */  =  1400.0; // frequency modulator.
			data->p[9] /* u_c     */  =  0.1  ;
			data->p[10] /* u_v     */ =  0.5  ;
			data->p[11] /* usi_c   */ =  0.01 ;
			data->p[12] /* k       */ =  5.0  ;
			data->IV = 0.0;
			data->delta_V = 0.0*gjc*(old_v0 - old_v1);
		    CVode(cvode_mem, tout, y, &t, CV_NORMAL);
usr_dvdt[1] = data->dvdt;
   usr_u0_loc[1] = Ith(y,1);
   usr_u1_loc[1] = Ith(y,2);  
   usr_u2_loc[1] = Ith(y,3);  

// do your measurement of cell # 0 here.
// dvdtmax during this oscillation.
if(dvdtmax1 < usr_dvdt[1]) dvdtmax1 = usr_dvdt[1];
// minimum.
if( (usr_dvdt[1]>0.0) && (fabs(usr_dvdt[1])<0.001) && (old_v1 < usr_u0_loc[1]) && (mes1_done==-1) && (usr_u0_loc[1]<0.2)){

cl1[n1] = sstime - t10; n1++; /******************************/

t10 = sstime; v10 = usr_u0_loc[1];
mes1_done = 1;
}
if( (usr_dvdt[1]<0.0) && (fabs(usr_dvdt[1])<0.001) && (old_v1 > usr_u0_loc[1]) && (mes1_done==1) && (usr_u0_loc[1]>0.2) ){

if(count1>10){

str = malloc(32*sizeof(char));
sprintf(str,"cell1data.%d",atoi(argv[1]));
output = fopen(str,"a+");
free(str);

// q_value and dvdt_control are for the SAN cell, the para cells q is not altered.
//fprintf(output,"%d\t%d\t%f\t%f\t%f\t%f\t%f\n",file_counter, 1, gjc, q_value, data->dvdt_control, cl1[n1-1], dvdtmax1);
fclose(output);
}

dvdtmax1 = -100000000.0;
mes1_done = -1;
count1++;
}

  N_VDestroy_Serial(y);
  CVodeFree(&cvode_mem);
free(data);
//***************************************************************************************
//***************************************************************************************
    /* output */
    if(iout%10==0&&sstime>(NOUT*DELTAT*0.8) ){ // last 100 s only.
	/* Line diagrams. */
	if(doutput==1){
	str = malloc(32*sizeof(char));
	sprintf(str,"output.dat.%d",file_counter);
	output = fopen(str,"a+");
	free(str);
	fprintf(output,"%10.10f\t",sstime);
	fprintf(output,"%10.10f\t",gjc);
	fprintf(output,"%10.10f\t%10.10f\t",usr_u0_loc[0],usr_u0_loc[1]);
	fprintf(output,"\n");
	fclose(output);
	}
    } // end of output.

      sstime = (double)iout*DELTAT;
      tout += DELTAT;

if(count0>110&&count1>100) break;

} // end of time loop.

// cell 0
mean0 = 0.0;
max_cl0 = -10000000.0;
min_cl0 =  10000000.0;
for(iout = 10; iout < n0; iout++){
mean0 += cl0[iout];
if(max_cl0 <= cl0[iout]) max_cl0 = cl0[iout];
if(min_cl0 >= cl0[iout]) min_cl0 = cl0[iout];
}

if((n0-10.0)>0)
mean0 = mean0/((double)(n0-10.0));
else
mean0 = -1000.0;

std0 = 0.0;
for(iout = 10; iout < n0; iout++){
std0 = std0 + (cl0[iout] - mean0)*(cl0[iout] - mean0);
}
if((n0-10.0)>0)
std0 = sqrt(std0)/((double)n0-10.0);
else
std0 = -1000.0;

// bubble sort.
double temp;
for(count0 = 10; count0 < n0-1; count0++){
	for(count1 = 10; count1 < n0 - 1 - count0; count1++){
		if(cl0[count1]>cl0[count1+1])
		{
			temp = cl0[count1+1];
			cl0[count1+1] = cl0[count1];
			cl0[count1] = temp;
		}
	}
} // end of sort.

if((n0-10)>0.0)
number_of_cls0 = 1;
else
number_of_cls0 = -1000.0;

for(count0 = 10; count0 < n0-1; count0++){
		if( fabs(cl0[count0+1]-cl0[count0]) > 1.0 ) number_of_cls0++;
}

// cell 1
mean1 = 0.0;
max_cl1 = -10000000.0;
min_cl1 =  10000000.0;
for(iout = 10; iout < n1; iout++){
mean1 += cl1[iout];
if(max_cl1 <= cl1[iout]) max_cl1 = cl1[iout];
if(min_cl1 >= cl1[iout]) min_cl1 = cl1[iout];
}

if((n1-10.0)>0)
mean1 = mean1/((double)(n1-10.0));
else
mean1 = -1000.0;

std1 = 0.0;
for(iout = 10; iout < n1; iout++){
std1 = std1 + (cl1[iout] - mean1)*(cl1[iout] - mean1);
}
if((n1-10)>0)
std1 = sqrt(std1)/((double)n1-10.0);
else
std1 = -1000.0;

// bubble sort.
for(count0 = 10; count0 < n1-1; count0++){
	for(count1 = 10; count1 < n1 - 1 - count0; count1++){
		if(cl1[count1]>cl1[count1+1])
		{
			temp = cl1[count1+1];
			cl1[count1+1] = cl1[count1];
			cl1[count1] = temp;
		}
	}
} // end of sort.

if((n1-10)>0.0)
number_of_cls1 = 1;
else
number_of_cls1 = -1000.0;

for(count0 = 10; count0 < n1-1; count0++){
		if( fabs(cl1[count0+1]-cl1[count0]) > 1.0 ) number_of_cls1++;
}

// now put it all into a file.
str = malloc(32*sizeof(char));
sprintf(str,"mydata.data.%d",atoi(argv[1]));
output = fopen(str,"a+");
free(str);
fprintf(output,"%d\t%f\t%f",   file_counter, gjc,     q_value               );
fprintf(output,"\t%f\t%f\t%d",               mean0,   std0,   number_of_cls0);
fprintf(output,"\t%f\t%f",                   min_cl0, max_cl0               );
fprintf(output,"\t%f\t%f\t%d",               mean1,   std1,   number_of_cls1);
fprintf(output,"\t%f\t%f\n",                 min_cl1, max_cl1);
fclose(output);

} // end of file_counter

  return 0;
} // end of cell main.

/*
 * error weights. this does not help much with eps_3_0
 */
static int ewt(N_Vector u, N_Vector w, void *user_data)
{
  int i;
  for (i=1; i<=NEQ; i++) Ith(w,i) = 1.0/(RTOL * ABS(Ith(u,i)) + ATOL);  
  return(0);
}

