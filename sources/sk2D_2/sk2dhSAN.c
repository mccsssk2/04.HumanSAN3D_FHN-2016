/*
26 August. Morph of mouse sk2d.c for human SAN 2D calculations.

The macro-reentry becomes very stable when interference from the SEPs is removed.
This is done by reducing diffusion at the SAN-atrial border (see RHS function).

In the macro-rentry case we do:
1. Reduce the SAN-atrial diffusion (RHS function).
2. Do 24 runs with various SAN atrial AP properties and overall diffusion properties.
*/

static char help[] = "Morph of sk2d.c from the Physoc 2015 2D SAN.\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscvec.h>
#include <petscsys.h>
#include <petscviewer.h> // this is for the VTK binary output in 3D.

/* my standard sundials headers, constants, and macros that will fit Petsc.
   Petsc cannot handle cvodes, at least not my installations. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>       /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>            /* prototype for CVDense                */
#include <sundials/sundials_dense.h>      /* definitions DlsMat DENSE_ELEM        */
#include <sundials/sundials_types.h>      /* definition of type realtype          */
#include <sundials/sundials_math.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ        */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ        */

#define RTOL        RCONST(1.0e-6)   	  /* scalar relative tolerance            */
#define ATOL        RCONST(1.0e-3)        /* scalar absolute tolerance components */

#define MAXSTEPS    50000
#define ZERO        RCONST(0.0)

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ        */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ        */

#define NEQ         3                     /* number of reaction variables. This is FK3 */

/* Problem Constants */
#define DELTAT       0.075                /* time step                              */
#define amp          0.500                /* stimulation amplitude                  */
#define duration     1.000                /* stimulus duration                      */

#define NUMPARAMS    13                   /* ODE parameters. Often or all the time I want to specify these as input */

#define pcl          1000.0                /* period of pacing                           */

#define DX           0.23496                 /* X internode spacing                       */
#define DY           0.23496                 /* Y internode spacing. This is the Z direction of the 3D hSAN     */

#define usr_MX 132  // this needs to be read from the geometry data.
#define usr_MY 166  // this needs to be read from the geometry data.
#define NBS    5    // 2D arrays have 5 neighbours, itself, and 4 surrounding it in a standard 1st order FD stencil.

#define DATA   3 // This model has 3 bits of data.

#define SEPS_low 0 // 0 for low, 1 for not low. This is in the diffusion function. case 0 gives ALL re-entry stuck to SAN heterogeniety.

/* User-defined data structures and routines           */
/* AppCtx: used by FormIFunction() and FormIJacobian() */
typedef struct {
  DM        da; // DM instance
/* All model data such as geomtery, fibers, cell type data, extracellular distributions must go through the context. */
  PetscInt  ***geometry; // 2D models have 3D array geometry. 0 = itself, 1 is y+1, 2 is x+1, 3 is y-1, 4 is x-1.
  PetscReal ***data;  // 2D models have 3D data. rightmost is however many data you have.
  PetscReal DISTMAX;
} AppCtx;

/* sundials contexts. */
// atrial
typedef struct {
	realtype IV;             // stimulus
	realtype lap_val;        // the laplacian.
	realtype p[NUMPARAMS];
	realtype dvdt;
	realtype factor;
	realtype q;
} *UserData;

// SAN
typedef struct {
	realtype dvdt;
	realtype q;
} *UserData_SAN;


// petsc functions.
extern PetscInt FormIFunction(TS,PetscReal,Vec,Vec,Vec,void*);

/* SUNDIALS functions called by the cvode solver                                      */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{

 realtype Y[NEQ];
 realtype dY[NEQ];

 realtype IV;
 UserData data;
 data = (UserData) user_data;

 int i;

 IV = data->IV;

  for(i=0;i<NEQ;i++)
    Y[i] = Ith(y,i+1);

// FK basal parameters.
/* data->p[0] */ realtype tau_d      = data->p[0];  //   0.25;
/* data->p[1] */ realtype tau_r      = data->p[1];  //   33  ;
/* data->p[2] */ realtype tau_si     = data->p[2];  //   30  ;
/* data->p[3] */ realtype tau_0      = data->p[3];  //   12.5;
/* data->p[4] */ realtype tau_vp     = data->p[4];  //   3.33;
/* data->p[5] */ realtype tau_vm1    = data->p[5];  //   1250;
/* data->p[6] */ realtype tau_vm2    = data->p[6];  //   19.6;
/* data->p[7] */ realtype tau_wp     = data->p[7];  //   870 ;
/* data->p[8] */ realtype tau_wm     = data->p[8];  //   41  ;
/* data->p[9] */ realtype u_c        = data->p[9];  //   0.13;
/* data->p[10] */ realtype u_v       = data->p[10]; //   0.04;
/* data->p[11] */ realtype usi_c     = data->p[11]; //   0.85;
/* data->p[12] */ realtype k         = data->p[12]; //   1.0 ;

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

  realtype Jfi, Jso, Jsi, tau_vm;

  if (Y[0]<u_c) {
    Jfi    =  0.0;
    Jso    =  Y[0]/(tau_0);
    tau_vm = (Y[0] > u_v) ? (tau_vm1) : (tau_vm2);
    dY[1]  = (1.0 - Y[1])/tau_vm;
    dY[2]  = (1.0 - Y[2])/tau_wm;
  } else {
    Jfi    = - Y[1]/(tau_d)*(1.-Y[0])*(Y[0]-u_c);
    Jso    =   1.0 /(tau_r);
    dY[1]  = - Y[1]/tau_vp;
    dY[2]  = - Y[2]/tau_wp;
  }
    Jsi    = - Y[2]/(2.0*tau_si)*(1.0 + tanh(k*(Y[0]-usi_c)));
  dY[0]    = - Jfi - Jso - Jsi + IV;

// printf("%f %f %f\n",Y[0],Y[1],Y[2]);

 for(i=0;i<NEQ;i++)
  Ith(ydot,i+1) = dY[i];

data->dvdt = dY[0];

 return(0);
}

// these are error weights in the Sundials solvers.
static int ewt(N_Vector y, N_Vector w, void *user_data);

/* For Random Numbers */
/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

int solution_broke = 0, somerandomnumber = 23334556;

/*
The random number generator.
*/

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/************************************************************************/

int main(int argc,char **argv)
{
/* Program reduced to use of FD with colouring. */
  TS             ts;                   /* nonlinear solver          */
  Vec            u, r, noscills;                  /* solution, residual vectors. You need the r for SNES */
  Mat            J,Jmf = NULL;         /* Jacobian matrices          */
  PetscInt       maxsteps = 500;      /* iterations for convergence */
  DM             da;
  AppCtx         user;                  /* user-defined work context  */
  SNES           snes; // scalable nonlinear equations solvers (snes)
//  int *usr_mybase_x, *usr_mybase_y, *usr_mysize_x, *usr_mysize_y, *usr_myblocksize; // I will use the same language as Petsc and Sundials.
  PetscInt       mybase_x, mybase_y, mysize_x, mysize_y, usr_i, usr_j, time_int, total_time_int;
  PetscInt       file_Counter = 0;
  PetscScalar    **u_localptr;

  double sanDiff_value = 0.0;

// MPI declarations.
  PetscInt    size = 0;                   /* Petsc/MPI                               */
  PetscInt    rank = 0;

/* Standard declarations */
/* Sundials out of time loop declarations. All my variables start with usr_···      */
  FILE      *geometry, *spiralFile, *random;
  realtype  usr_stimtime, usr_time, usr_t, usr_tout;
  double    **usr_u0_loc, **usr_u1_loc, **usr_u2_loc, **usr_u0_loc_old;
  double    **usr_num_Oscillations;
  int       int1, int2, int3, int4, int5, int6, int7;
  realtype  dist, q_value, distmax;
  PetscReal distx;
  double    radius;
  double    tmp_u0, tmp_u1, tmp_u2;
  char      str[1000]; // so you do not keep creating and destroying this.

  PetscInitialize(&argc,&argv,(char*)0,help);
  /* Initialize user application context */
  user.da           = NULL;

/*
cases: These depend on the job arrary number from 1 to 24.
Here are the base values.
*/
int    control_short_ap      = 0;       // 0 is for control (long) AP in atrium, 1 is for short AP in atrium.
int    san_control_fast_slow = 0;       // 0 = Control pacing (q=600), 1 = fast pacing (q=200), 2 = slow pacing rate (q=1200).
int    if_no_border          = 0;       // 0 = with border, 1 = remove border in this program.
// diffusions of all tissue types possible.
double diffusion_atrial      = 0.325;   // tissue type 1. This is validated from FIND A HSAN paper.

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

/******************************************************************************/
/*******************write the 24  cases****************************************/
/******************************************************************************/

if( atoi(argv[1])==1 ){ // Control
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 0    ;
}

if( atoi(argv[1])==2 ){ // fast SAN
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 1    ;
}

if( atoi(argv[1])==3 ){ // slow SAN
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 2    ;
}

if( atoi(argv[1])==4 ){ // Control SAN, short atrial AP
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 0    ;
}

if( atoi(argv[1])==5 ){ // fast SAN, short atrial AP
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 1    ;
}

if( atoi(argv[1])==6 ){ // slow SAN, short atrial AP
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 2    ;
}


if( atoi(argv[1])==7 ){ // low GJC,  Control
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 0    ;
}

if( atoi(argv[1])==8 ){ // low GJC,  fast SAN
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 1    ;
}

if( atoi(argv[1])==9 ){ // low GJC,  slow SAN
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 2    ;
}

if( atoi(argv[1])==10 ){ // low GJC,  Control SAN, short atrial AP
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 0    ;
}

if( atoi(argv[1])==11 ){ // low GJC,  fast SAN, short atrial AP
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 1    ;
}

if( atoi(argv[1])==12 ){ // low GJC,  slow SAN, short atrial AP
	// Border, no border
	if_no_border          = 0    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 2    ;
}


if( atoi(argv[1])==13 ){ // no border, Control
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 0    ;
}

if( atoi(argv[1])==14 ){ // no border, fast SAN
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 1    ;
}

if( atoi(argv[1])==15 ){ // no border, slow SAN
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 2    ;
}

if( atoi(argv[1])==16 ){ // no border, Control SAN, short atrial AP
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 0    ;
}

if( atoi(argv[1])==17 ){ // no border, fast SAN, short atrial AP
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 1    ;
}

if( atoi(argv[1])==18 ){ // no border, slow SAN, short atrial AP
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 2    ;
}

if( atoi(argv[1])==19 ){ // no border, low GJC,  Control
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 0    ;
}

if( atoi(argv[1])==20 ){ // no border, low GJC,  fast SAN
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 1    ;
}

if( atoi(argv[1])==21 ){ // no border, low GJC,  slow SAN
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 0    ;
	san_control_fast_slow = 2    ;
}

if( atoi(argv[1])==22 ){ // no border, low GJC,  Control SAN, short atrial AP
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 0    ;
}

if( atoi(argv[1])==23 ){ // no border, low GJC,  fast SAN, short atrial AP
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 1    ;
}

if( atoi(argv[1])==24 ){ // no border, low GJC,  slow SAN, short atrial AP
	// Border, no border
	if_no_border          = 1    ;
	// diffusion control (0.325) or low (0.1*0.325).
	diffusion_atrial      = 0.1*0.325;
	// EP
	control_short_ap      = 1    ;
	san_control_fast_slow = 2    ;
}

/******************************************************************************/
/******************************************************************************/

// start time.
time_int = 0; usr_time = 0.0; usr_stimtime = 0.0;
total_time_int = (int)(5000.0/DELTAT); // 5 seconds of activity.
// total_time_int = (int)(50/DELTAT);
for(time_int = 0; time_int < total_time_int; time_int++){ // I decided to do 2 seconds, which takes 1.6 hours on 4 procs.
   DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,usr_MX,usr_MY,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da); 
   user.da = da;
   DMCreateGlobalVector(da,&u); 
   VecDuplicate(u,&r);
   VecDuplicate(u,&noscills);

		/* size_y is the width of the array. Do this once. */
		if(time_int < 1){
			   DMDAGetCorners(da,&mybase_x,&mybase_y,NULL,&mysize_x,&mysize_y,NULL); 

				usr_u0_loc       = (double **) calloc(mysize_y, sizeof(double*));
				usr_u0_loc_old   = (double **) calloc(mysize_y, sizeof(double*));
				usr_u1_loc       = (double **) calloc(mysize_y, sizeof(double*));
				usr_u2_loc       = (double **) calloc(mysize_y, sizeof(double*));
			usr_num_Oscillations = (double **) calloc(mysize_y, sizeof(double*));
/* the usr_mysize_ and usr_mybase_ arrays are not the same on all processors. */
					for (usr_j = 0; usr_j < mysize_y; usr_j++){ 
						usr_u0_loc[usr_j]     = (double*) calloc(mysize_x,sizeof(double));
						usr_u0_loc_old[usr_j] = (double*) calloc(mysize_x,sizeof(double));
						usr_u1_loc[usr_j]     = (double*) calloc(mysize_x,sizeof(double));
						usr_u2_loc[usr_j]     = (double*) calloc(mysize_x,sizeof(double));
			            usr_num_Oscillations[usr_j] = (double *) calloc(mysize_x, sizeof(double));
					}

			/* Declare the geometry memory. This declares memory on each proc. */
			user.geometry = (PetscInt ***) calloc(mysize_y, sizeof(PetscInt**));
			for (usr_j = 0; usr_j < mysize_y; usr_j++){
				user.geometry[usr_j] = (PetscInt**) calloc(mysize_x,sizeof(PetscInt*));
				for (usr_i = 0; usr_i < mysize_x; usr_i++)
					user.geometry[usr_j][usr_i] = (PetscInt*) calloc(NBS, sizeof(PetscInt));				
			}

			/* Declare the data memory. This declares memory on each proc. */
			user.data = (PetscReal ***) calloc(mysize_y, sizeof(PetscReal**));
			for (usr_j = 0; usr_j < mysize_y; usr_j++){
				user.data[usr_j] = (PetscReal**) calloc(mysize_x,sizeof(PetscReal*));
				for (usr_i = 0; usr_i < mysize_x; usr_i++)
					user.data[usr_j][usr_i] = (PetscReal*) calloc(DATA, sizeof(PetscReal));				
			}


/* Initialising the random number generator is important. So this will have a different value on each processor. */
random = fopen("/dev/urandom","rb");
fread(&somerandomnumber, sizeof(unsigned int), 1, random);
fclose(random);
init_genrand(somerandomnumber);

				/* read in your geometry here.  		*/
				/* This reads the ASCII geometry on each proc. 
				The file is opened for reading on each proc.    */
				user.DISTMAX = -10.0;
					geometry = fopen("hSAN2015_2D.geom","r");
			while(fscanf(geometry,"%d %d %d %d %d %d %d %lf %lf %lf",&int1,&int2,&int3,&int4,&int5,&int6,&int7,&dist,&q_value, &distmax)!=EOF){

				usr_j = (PetscInt)int1 - mybase_y; usr_i = (PetscInt)int2  - mybase_x;
				if(usr_j>=0&&usr_j<mysize_y&&usr_i>=0&&usr_i<mysize_x){
				user.geometry[usr_j][usr_i][0] = (PetscInt)int3; // 5 types of neighbours: me, top, right, bottom, left.

				user.geometry[usr_j][usr_i][1] = (PetscInt)int4;
				user.geometry[usr_j][usr_i][2] = (PetscInt)int5;
				user.geometry[usr_j][usr_i][3] = (PetscInt)int6;
				user.geometry[usr_j][usr_i][4] = (PetscInt)int7;

				// removed border here.
				if(if_no_border==1&&user.geometry[usr_j][usr_i][0]==4) // without border
				user.geometry[usr_j][usr_i][0] = (PetscInt)(1);
				if(if_no_border==1&&user.geometry[usr_j][usr_i][1]==4)
				user.geometry[usr_j][usr_i][1] = (PetscInt)(1);
				if(if_no_border==1&&user.geometry[usr_j][usr_i][2]==4)
				user.geometry[usr_j][usr_i][2] = (PetscInt)(1);
				if(if_no_border==1&&user.geometry[usr_j][usr_i][3]==4)
				user.geometry[usr_j][usr_i][3] = (PetscInt)(1);
				if(if_no_border==1&&user.geometry[usr_j][usr_i][4]==4)
				user.geometry[usr_j][usr_i][4] = (PetscInt)(1);
// model data.
				user.data[usr_j][usr_i][0]     = (PetscReal)dist; // output this once.

				user.data[usr_j][usr_i][1]     = 1.0; // this is redundant.
				if(user.geometry[usr_j][usr_i][0]==2){
				if(san_control_fast_slow==0) user.data[usr_j][usr_i][1] = (7.00 + (1.0 - 0.5*genrand_real1())*2.0); // Control.
				if(san_control_fast_slow==1) user.data[usr_j][usr_i][1] = (5.00 + (1.0 - 0.5*genrand_real1())*2.0); // fast.
				if(san_control_fast_slow==2) user.data[usr_j][usr_i][1] = (20.0 + (1.0 - 0.5*genrand_real1())*2.0); // slow.
				}

				user.DISTMAX                   = (PetscReal)distmax;

// diffusion because there are first order terms in the heterogeneous data. Also, this opens up the random number for me here.
// me
		distx = user.data[usr_j][usr_i][0];
		sanDiff_value = (0.0000325 + diffusion_atrial * distx*distx );
		
		if(user.geometry[usr_j][usr_i][0]==0) // nothing 
		       user.data[usr_j][usr_i][2] = (PetscReal)(0.0);
		if(user.geometry[usr_j][usr_i][0]==1) // atrial uniform
		       user.data[usr_j][usr_i][2] = (PetscReal)(diffusion_atrial);
		if(user.geometry[usr_j][usr_i][0]==2) // SAN
		       user.data[usr_j][usr_i][2] = (PetscReal) sanDiff_value;
		if(user.geometry[usr_j][usr_i][0]==3) // apoptosis
		       user.data[usr_j][usr_i][2] = ((PetscReal) sanDiff_value); // a parameter.
		if(user.geometry[usr_j][usr_i][0]==4) // insulating border.
		       user.data[usr_j][usr_i][2] = (PetscReal)(0.0); // this is not going to matter any more.
				} // end of if statement
			} // end of while statement


// write the q. q is a real, usr_u0_loc is a real. use that for now to write the q.
		for(usr_j = 0; usr_j < mysize_y; usr_j++)
		for(usr_i = 0; usr_i < mysize_x; usr_i++)
			usr_u0_loc[usr_j][usr_i] = user.data[usr_j][usr_i][1];
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_j = 0; usr_j < mysize_y; usr_j++)
	for(usr_i = 0; usr_i < mysize_x; usr_i++)
		u_localptr[mybase_y + usr_j][mybase_x + usr_i] = usr_u0_loc[usr_j][usr_i];
    DMDAVecRestoreArray(da,u,&u_localptr);
			sprintf(str,"q_valueMouse2D.vts");
			PetscViewer viewer4;
			PetscViewerCreate(PETSC_COMM_WORLD, &viewer4);
			PetscViewerSetType(viewer4, PETSCVIEWERVTK);
			PetscViewerFileSetName(viewer4, str);
			VecView(u, viewer4);
			PetscViewerDestroy(&viewer4);

	/* Sundials initial conditions for resting state.                                                             */
		for(usr_j = 0; usr_j < mysize_y; usr_j++)
		for(usr_i = 0; usr_i < mysize_x; usr_i++){
			usr_num_Oscillations[usr_j][usr_i] = 0.0;
			if(user.geometry[usr_j][usr_i][0]>0&&user.geometry[usr_j][usr_i][0]<4)
			{ usr_u0_loc[usr_j][usr_i] = 0.0; usr_u1_loc[usr_j][usr_i] = 0.0; usr_u2_loc[usr_j][usr_i] = 0.0; }
		else    { usr_u0_loc[usr_j][usr_i] = usr_u1_loc[usr_j][usr_i] = usr_u2_loc[usr_j][usr_i] = -0.01;         }
		}

// random initial conditions.
/*
		for(usr_j = 0; usr_j < mysize_y; usr_j++)
		for(usr_i = 0; usr_i < mysize_x; usr_i++){
			if(user.geometry[usr_j][usr_i][0]==2) // SAN only.
			{ usr_u0_loc[usr_j][usr_i] = genrand_real1(); usr_u1_loc[usr_j][usr_i] = genrand_real1(); usr_u2_loc[usr_j][usr_i] = genrand_real1(); }
		else    { usr_u0_loc[usr_j][usr_i] = usr_u1_loc[usr_j][usr_i] = usr_u2_loc[usr_j][usr_i] = -0.01;         }
		}
*/

// spiral wave initial conditions. This file has initial conditions for all the FD points in this geometry.

spiralFile = fopen("spiralInit2D.data","r");
while(fscanf(spiralFile,"%d %d %lf %lf %lf",&int1, &int2, &tmp_u0, &tmp_u1, &tmp_u2)!=EOF){
   usr_j = (PetscInt)int1 - mybase_y; usr_i = (PetscInt)int2  - mybase_x;
				if(usr_j>=0&&usr_j<mysize_y&&usr_i>=0&&usr_i<mysize_x){
					if(user.geometry[usr_j][usr_i][0]==2&&user.geometry[usr_j][usr_i][1]==2&&user.geometry[usr_j][usr_i][2]==2&&user.geometry[usr_j][usr_i][3]==2&&user.geometry[usr_j][usr_i][4]==2)   { // for macrorenetry, it is 1, for micro, it is 2. for micro, you need fibrosis as well. Simple fibrosis here, patches from the MATLAB programs.
					 	usr_u0_loc[usr_j][usr_i] = (double)tmp_u0;
						usr_u1_loc[usr_j][usr_i] = (double)tmp_u1;
						usr_u2_loc[usr_j][usr_i] = (double)tmp_u2;
					}
				}
}
fclose(spiralFile);

// write initial conditions.
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_j = 0; usr_j < mysize_y; usr_j++)
	for(usr_i = 0; usr_i < mysize_x; usr_i++)
		u_localptr[mybase_y + usr_j][mybase_x + usr_i] = usr_u0_loc[usr_j][usr_i];
    DMDAVecRestoreArray(da,u,&u_localptr);
			sprintf(str,"InitialConditions_HSAN2D.vts");
			PetscViewer viewer5;
			PetscViewerCreate(PETSC_COMM_WORLD, &viewer5);
			PetscViewerSetType(viewer5, PETSCVIEWERVTK);
			PetscViewerFileSetName(viewer5, str);
			VecView(u, viewer5);
			PetscViewerDestroy(&viewer5);

		} // end of time_int < 1.

	for(usr_j = 0; usr_j < mysize_y; usr_j++) // these usr_j and usr_i are for local locations.
	for(usr_i = 0; usr_i < mysize_x; usr_i++){
		usr_u0_loc_old[usr_j][usr_i] = usr_u0_loc[usr_j][usr_i];
		if(user.geometry[usr_j][usr_i][0]>0&&user.geometry[usr_j][usr_i][0]<4){
	/*******************************************************************************/
		    N_Vector y;                    /* The state vector                 */
		    void *cvode_mem;               /* Sundials memory pointer          */
		    UserData data;                 /* Sundials data structure instance */ 
		    /* Create serial vector of length NEQ for I.C. and abstol          */
		    cvode_mem = NULL;
		    data      = NULL;
		    y         = NULL;
		    y         = N_VNew_Serial(NEQ);
		    data      = (UserData) malloc(sizeof *data);
		    /* user specified Sundials data.                */
			data->IV = 0.0;
		    data->factor = 1.0; /* user specified Sundials data.                */
			data->q = user.data[usr_j][usr_i][1]; // q changes only in the SAN part.


// 0.025 100.0 45 16.5 16.5 12500 75 1.5 50 0.2 0.5 0.01 2.5 6.0
// 0.025 100.0 45.0 16.5 16.50 12500.0 75.0 1.5 50.0 0.2 0.5 0.01 2.5 6.0
if(user.geometry[usr_j][usr_i][0]==2){ // SAN part


//  0.05 130.0 50.0 12.50 10.0 18.2 18.2 1.5 600.0 0.1 0.5 0.01 5.0 2
		data->p[0] /* double tau_d     */  =    0.20; // 0.050 ;
		data->p[1] /* double tau_r     */  =    130.0 ;
		data->p[2] /* double tau_si    */  =    50.00 ;
		data->p[3] /* double tau_0     */  =    12.5  ;
		data->p[4] /* double tau_vp    */  =    10.0  ;
		data->p[5] /* double tau_vm1   */  =    18.2  ;
		data->p[6] /* double tau_vm2   */  =    18.2  ;
		data->p[7] /* double tau_wp    */  =    1.5   ;

// this needs to be revised carefully for pacing at 850. I got it, do a bit more checking. revision 1 done.
		data->p[8] /* double tau_wm    */  =    100.0*data->q; // control is 600, fast is 200, slow is 1200 (max. vales).

		data->p[9] /* double u_c       */  =    0.1   ;
		data->p[10] /* double u_v       */ =    0.5   ; // does not matter because tau_vm1 and tau_vm2 are equal.
		data->p[11] /* double usi_c     */ =    0.01  ;
		data->p[12] /* double k         */ =    5.0   ;
}
else if(user.geometry[usr_j][usr_i][0]==1){ // atrial part.
	// Control
	if(control_short_ap==0){
			data->p[0] /* double tau_d     */  =    0.2000;
			data->p[1] /* double tau_r     */  =    130.0 ;
			data->p[2] /* double tau_si    */  =    127.0 ;
			data->p[3] /* double tau_0     */  =    12.5  ;
			data->p[4] /* double tau_vp    */  =    10.0  ;
			data->p[5] /* double tau_vm1   */  =    18.2  ;
			data->p[6] /* double tau_vm2   */  =    18.2  ;
			data->p[7] /* double tau_wp    */  =    1020.0;
			data->p[8] /* double tau_wm    */  =    80.0  ;
			data->p[9] /* double u_c       */  =    0.3   ;
			data->p[10] /* double u_v       */ =    0.5   ; // does not matter because tau_vm1 and tau_vm2 are equal.
			data->p[11] /* double usi_c     */ =    0.85  ;
			data->p[12] /* double k         */ =    10.0  ;
	}
	if(control_short_ap==1){
	// Short AP. AF.
	//   1     2     3     4    5    6    7     8     9   10  11   12  13   14  15
	// 0.0200 130.0 310.0 12.5 10.0 10.0 30.0 1020.0 352. 0.3 0.5 0.85 10.0 1.0 1
			data->p[0] /* double tau_d     */  =    0.2000;
			data->p[1] /* double tau_r     */  =    40.0  ;
			data->p[2] /* double tau_si    */  =    310.0 ;
			data->p[3] /* double tau_0     */  =    12.5  ;
			data->p[4] /* double tau_vp    */  =    10.0  ;
			data->p[5] /* double tau_vm1   */  =    10.0  ;
			data->p[6] /* double tau_vm2   */  =    30.0  ;
			data->p[7] /* double tau_wp    */  =    10.0; // AF
			data->p[8] /* double tau_wm    */  =    352.0 ;
			data->p[9] /* double u_c       */  =    0.3   ;
			data->p[10] /* double u_v       */ =    0.5   ; // does not matter because tau_vm1 and tau_vm2 are equal.
			data->p[11] /* double usi_c     */ =    0.85  ;
			data->p[12] /* double k         */ =    10.0  ;
	}
} // end of atrial if
// test this apoptosis cell before using it. It must be inexcitable, marginal AP producing, and must not produce propagations.
else if(user.geometry[usr_j][usr_i][0]==3){ // Apoptosis: Use this to make the SAN-atrial junctions inactive.
		data->p[0] /* double tau_d     */  =    10.000; // low upstroke.
		data->p[1] /* double tau_r     */  =    130.0 ;
		data->p[2] /* double tau_si    */  =    127.0 ;
		data->p[3] /* double tau_0     */  =    12.5  ;
		data->p[4] /* double tau_vp    */  =    10.0  ;
		data->p[5] /* double tau_vm1   */  =    18.2  ;
		data->p[6] /* double tau_vm2   */  =    18.2  ;
		data->p[7] /* double tau_wp    */  =    100.00;
		data->p[8] /* double tau_wm    */  =    80.0  ;
		data->p[9] /* double u_c       */  =    0.05  ;
		data->p[10] /* double u_v       */ =    1.0   ; // does not matter because tau_vm1 and tau_vm2 are equal.
		data->p[11] /* double usi_c     */ =    1.00  ;
		data->p[12] /* double k         */ =    10.0  ;
}

		    usr_t     = 0; usr_tout = DELTAT;
		    Ith(y,1)  = usr_u0_loc[usr_j][usr_i];
		    Ith(y,2)  = usr_u1_loc[usr_j][usr_i];
		    Ith(y,3)  = usr_u2_loc[usr_j][usr_i];
		    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
		    CVodeInit(cvode_mem, f, 0.0, y);
		    CVodeWFtolerances(cvode_mem, ewt);
		    CVodeSetUserData(cvode_mem, data);
		    CVodeSStolerances(cvode_mem, RTOL, ATOL);
		    CVDense(cvode_mem, NEQ);
//		    CVodeSetMaxStep(cvode_mem, DELTAT); // dont set the max step unless necessary
		    CVode(cvode_mem, usr_tout, y, &usr_t, CV_NORMAL);

		    usr_u0_loc[usr_j][usr_i] = Ith(y,1);
		    usr_u1_loc[usr_j][usr_i] = Ith(y,2);
		    usr_u2_loc[usr_j][usr_i] = Ith(y,3);

		    N_VDestroy_Serial(y);
		    CVodeFree(&cvode_mem);
		    free(data);
	/*******************************************************************************/
		} /* end of geometry if.                                                        */
} // end of Sundials loop.

// petsc part, part II
   TSCreate(PETSC_COMM_WORLD,&ts); 
   TSSetProblemType(ts,TS_NONLINEAR); 
   TSSetType(ts,TSBEULER); 
   TSSetDM(ts,da); 
   TSSetIFunction(ts,r,FormIFunction,&user); 
   TSSetDuration(ts,maxsteps,DELTAT); // duration of simulation.
/********** put your usr_u0_loc into u here. **************************************/
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_j = 0; usr_j < mysize_y; usr_j++)
	for(usr_i = 0; usr_i < mysize_x; usr_i++)
		u_localptr[mybase_y + usr_j][mybase_x + usr_i] = usr_u0_loc[usr_j][usr_i];
    DMDAVecRestoreArray(da,u,&u_localptr);

/**********************************************************************************/
   TSSetSolution(ts,u); 
   TSSetInitialTimeStep(ts,0.0,DELTAT); 
   DMSetMatType(da,MATAIJ); 
   DMCreateMatrix(da,&J); 
   TSGetSNES(ts,&snes); 
   MatCreateSNESMF(snes,&Jmf); 
   SNESSetJacobian(snes,Jmf,J,SNESComputeJacobianDefaultColor,0); 
   TSSetFromOptions(ts); 
   TSSolve(ts,u); 
/********** put your u back into usr_u0_loc here. ********************************/
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_j = 0; usr_j < mysize_y; usr_j++)
	for(usr_i = 0; usr_i < mysize_x; usr_i++){
		 usr_u0_loc[usr_j][usr_i] = u_localptr[mybase_y + usr_j][mybase_x + usr_i]; /* LHS is my 2D array, RHS is PetSc's */
if(usr_u0_loc[usr_j][usr_i]<-0.15||usr_u0_loc[usr_j][usr_i]>1.25){ printf("It failed.\n"); exit(1);}
}
    DMDAVecRestoreArray(da,u,&u_localptr);

	for(usr_j = 0; usr_j < mysize_y; usr_j++)
	for(usr_i = 0; usr_i < mysize_x; usr_i++){
		if(usr_u0_loc_old[usr_j][usr_i]<0.5&&usr_u0_loc[usr_j][usr_i]>=0.5)
		usr_num_Oscillations[usr_j][usr_i]++;
	}

    DMDAVecGetArray(da,noscills,&u_localptr);
	for(usr_j = 0; usr_j < mysize_y; usr_j++)
	for(usr_i = 0; usr_i < mysize_x; usr_i++){
		 u_localptr[mybase_y + usr_j][mybase_x + usr_i] = usr_num_Oscillations[usr_j][usr_i]; 
	}
    DMDAVecRestoreArray(da,noscills,&u_localptr);

/********************************************************************************/
// VTK and binary data output.
		if( (time_int%( (int)(1.0/DELTAT) )==0) &&(time_int>0.0*total_time_int) ){ // every ms, the DELTAT is now 0.05 because it kept failing.

			sprintf(str,"my_2d%d.vts",file_Counter);
			PetscViewer viewer;
			PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
			PetscViewerSetType(viewer, PETSCVIEWERVTK);
			PetscViewerFileSetName(viewer, str);
			VecView(u, viewer);
			PetscViewerDestroy(&viewer);

			// put out a binary every milli-second so that you can do your calculations as post processing.
			PetscViewer viewer2;
			sprintf(str,"my_2d%d.bin",file_Counter++);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_WRITE,&viewer2);
			VecView(u,viewer2);
			PetscViewerDestroy(&viewer2);

			sprintf(str,"nOSC.vts");
			PetscViewer viewer3;
			PetscViewerCreate(PETSC_COMM_WORLD, &viewer3);
			PetscViewerSetType(viewer3, PETSCVIEWERVTK);
			PetscViewerFileSetName(viewer3, str);
			VecView(noscills, viewer3);
			PetscViewerDestroy(&viewer3);

		}

// this is in your time loop for now.
   VecDestroy(&u); 
   VecDestroy(&r); 
   VecDestroy(&noscills);
   MatDestroy(&J); 
   MatDestroy(&Jmf); 
   TSDestroy(&ts); 
   DMDestroy(&da); 

usr_time = (PetscReal)time_int*DELTAT; usr_stimtime = (PetscReal)time_int*DELTAT;
if(usr_stimtime>pcl) usr_stimtime = 0.0;

} // end of time.

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   free(usr_u0_loc);
   free(usr_u1_loc);
   free(usr_u2_loc);
free(usr_u0_loc_old);
free(usr_num_Oscillations);

   PetscFinalize();
   PetscFunctionReturn(0);
}

/* FormIFunction = Udot - RHSFunction */
PetscInt FormIFunction(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx)
{
  AppCtx         *user=(AppCtx*)ctx;
  DM             da   = (DM)user->da;
  PetscInt       i,j,xs,ys,xm,ym, gj, gi;
  PetscReal      hx,hy,sx,sy;
  PetscReal      U0, U1, U2, U3, U4;
  PetscScalar    uxx,uyy, **uarray, **f, **udot;
  Vec            localU;

   PetscFunctionBeginUser;
   DMGetLocalVector(da,&localU);

  /*
     Scatter ghost points to local vector,using the 2-step process
        DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
   DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU); 
   DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU); 

  /* Get pointers to vector data */
   DMDAVecGetArray(da,localU,&uarray); 

   DMDAVecGetArray(da,F,&f); 
   DMDAVecGetArray(da,Udot,&udot); 


  /* Get local grid boundaries */
   DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); 

/* wrong. I started with ex15.c in the ts examples. fix it with Preconditioner friendly method. */
/* Compute function over the locally owned part of the grid                                     */
  for (j=ys; j<ys+ym; j++)
    for (i=xs; i<xs+xm; i++)
	f[j][i] = 0.0; // so that i do not have to deal with empty space.

    for (j=ys; j<ys+ym; j++) 
    for (i=xs; i<xs+xm; i++) {

gj = j - ys;
gi = i - xs;

hx = DX; sx = (user->data[gj][gi][2])/(hx*hx); // sx1 = user->data[gj][gi][7]/(2.0*hx);
hy = DY; sy = (user->data[gj][gi][2])/(hy*hy); // sy1 = user->data[gj][gi][8]/(2.0*hy);

//if(sep==1) x0 = 58; y0 = 18;  end;
//if(sep==2) x0 = 71; y0 = 53;  end;
//if(sep==3) x0 = 56; y0 = 104; end;
//if(sep==4) x0 = 31; y0 = 104; end;

//for i=1:1:maxxN
//	for j=1:1:maxyN
//	if( img3(j,i)==4)
//		radiusTemp = (i - x0)*(i-x0) + (j-y0)*(j-y0);
//		if(radiusTemp < sizeSquared) img5(j,i) = 1; end;
//	end;
//	end;
//end;

// SEP conduction is damaged. If there is no SEP, then it is just a border.
if(SEPS_low==0)
if( (user->geometry[gj][gi][0]==1)&&
    (user->geometry[gj][gi][1]==2||user->geometry[gj][gi][2]==2||user->geometry[gj][gi][3]==2||user->geometry[gj][gi][4]==2||
     user->geometry[gj][gi][1]==3||user->geometry[gj][gi][2]==3||user->geometry[gj][gi][3]==3||user->geometry[gj][gi][4]==3  ) )
{
	hx = DX; sx = 0.001*(user->data[gj][gi][2])/(hx*hx); // sx1 = user->data[gj][gi][7]/(2.0*hx);
	hy = DY; sy = 0.001*(user->data[gj][gi][2])/(hy*hy); // sy1 = user->data[gj][gi][8]/(2.0*hy);
}
		/* This is still hacky. I am not happy with the way boundaries are broken down into separate components. */
		if(user->geometry[gj][gi][0]>0&&user->geometry[gj][gi][0]<4){
			if(i>0&&j>0&&i<usr_MX-1&&j<usr_MY-1){

			                     U1 = uarray[j+1][i]; 
			U4 = uarray[j][i-1]; U0 = uarray[j  ][i]; U2 = uarray[j][i+1]; 
	                                     U3 = uarray[j-1][i]; 

		/* this is not the best no flux b.c., improve later. */
if(user->geometry[gj][gi][0]==1){
// 1 is close to 0, 2, 3 and 4. 2 and 3 are ok, but 0 and 4 are nothing.
		if(user->geometry[gj][gi][1]  ==0) 				     U1 = U3;
		if(user->geometry[gj][gi][3]  ==0) 				     U3 = U1;
		if( (user->geometry[gj][gi][1]==0)&&(user->geometry[gj][gi][3]==0) ) U3 = U1 = U0;
		if(user->geometry[gj][gi][2]  ==0) 				     U2 = U4;
		if(user->geometry[gj][gi][4]  ==0) 				     U4 = U2;
		if( (user->geometry[gj][gi][2]==0)&&(user->geometry[gj][gi][4]==0) ) U4 = U2 = U0;

		if(user->geometry[gj][gi][1]  ==4) 				     U1 = U3;
		if(user->geometry[gj][gi][3]  ==4) 				     U3 = U1;
		if( (user->geometry[gj][gi][1]==4)&&(user->geometry[gj][gi][3]==4) ) U3 = U1 = U0;
		if(user->geometry[gj][gi][2]  ==4) 				     U2 = U4;
		if(user->geometry[gj][gi][4]  ==4) 				     U4 = U2;
		if( (user->geometry[gj][gi][2]==4)&&(user->geometry[gj][gi][4]==4) ) U4 = U2 = U0;

		if( (user->geometry[gj][gi][1]==0)&&(user->geometry[gj][gi][3]==4) ) U3 = U1 = U0;
		if( (user->geometry[gj][gi][1]==4)&&(user->geometry[gj][gi][3]==0) ) U3 = U1 = U0;
		if( (user->geometry[gj][gi][2]==0)&&(user->geometry[gj][gi][4]==4) ) U4 = U2 = U0;
		if( (user->geometry[gj][gi][2]==4)&&(user->geometry[gj][gi][4]==0) ) U4 = U2 = U0;
}

if(user->geometry[gj][gi][0]==2||user->geometry[gj][gi][0]==3){
// 2 or 3 are close to 1, 3 and 4. 4 is non-conducting.

		if(user->geometry[gj][gi][1]  ==4) 				     U1 = U3;
		if(user->geometry[gj][gi][3]  ==4) 				     U3 = U1;
		if( (user->geometry[gj][gi][1]==4)&&(user->geometry[gj][gi][3]==4) ) U3 = U1 = U0;
		if(user->geometry[gj][gi][2]  ==4) 				     U2 = U4;
		if(user->geometry[gj][gi][4]  ==4) 				     U4 = U2;
		if( (user->geometry[gj][gi][2]==4)&&(user->geometry[gj][gi][4]==4) ) U4 = U2 = U0;

}

// if(user->geometry[gj][gi][0]==3||user->geometry[gj][gi][0]==4) U1 = U2 = U3 = U4 = U0 = -0.1; // lets see if this sticks.

			uxx = (-2.0*U0 + U4 + U2);
			uyy = (-2.0*U0 + U3 + U1);
// the convection term causes an instability. Maybe the accuracy of the derivative is not good.
			f[j][i] = udot[j][i] - (uxx*sx + uyy*sy ); // the udot is not zero.	
			}
			else
			{
			/* Boundary conditions. */
				if(j==usr_MY-1) f[j][i] = udot[j][i] + sy*(uarray[j][i] - uarray[j-1][i]); // these are first order approximations.
				if(i==usr_MX-1) f[j][i] = udot[j][i] + sx*(uarray[j][i] - uarray[j][i-1]);
				if(j==0)        f[j][i] = udot[j][i] + sy*( uarray[j][i] - uarray[j+1][i]);
				if(i==0)        f[j][i] = udot[j][i] + sx*( uarray[j][i] - uarray[j][i+1]);
			if((j==0)&&(i==0))               f[j][i] = udot[j][i] + 2.0*(sx + sy)*( uarray[j][i] - uarray[j+1][i+1]);
			if((j==0)&&(i==usr_MX-1))        f[j][i] = udot[j][i] + 2.0*(sx + sy)*( uarray[j][i] - uarray[j-1][i+1]);
			if((j==usr_MY-1)&&(i==0))        f[j][i] = udot[j][i] + 2.0*(sx + sy)*( uarray[j][i] - uarray[j-1][i+1]);
			if((j==usr_MY-1)&&(i==usr_MX-1)) f[j][i] = udot[j][i] + 2.0*(sx + sy)*( uarray[j][i] - uarray[j-1][i-1]);
		} // this else above is a hack, do properly.
	}
    } // end of for loops.

  /* Restore vectors */
   DMDAVecRestoreArray(da,localU,&uarray); 

   DMDAVecRestoreArray(da,F,&f);
   DMDAVecRestoreArray(da,Udot,&udot); 
   DMRestoreLocalVector(da,&localU); 
  PetscFunctionReturn(0);
}



/*
 * error weights. this does not help much with eps_3_0 Sundials function I like.
 */
static int ewt(N_Vector u, N_Vector w, void *user_data)
{
  int i;
  for (i=1; i<=NEQ; i++) Ith(w,i) = 1.0/(RTOL * ABS(Ith(u,i)) + ATOL);  
  return(0);
}

