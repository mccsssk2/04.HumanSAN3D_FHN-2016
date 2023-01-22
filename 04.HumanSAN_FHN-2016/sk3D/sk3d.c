/*
SK. 16 April 2015.

Each simulation does 4 cases of anatomy.

To make a 3D model of the human SAN based on HD geometry.
It is going to have an insulating layer, as well as a heterogeneous SAN.

June 14.
Now you got paranodal volume in the model. Test this out.

0 = no tissue
1 = SAN
2 = SAN
3 = paranodal oscillator
4 = atrial tissue
5 = adipose tissue
6 = insulating border

7 = inactive cells. This is to produce the intermmitent atrial pacing following Ostborne.
*/

static char help[] = "Human SAN 3D geometry.\n";

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

/*
cases:
control_short_ap      takes 0 or 1
san_control_fast_slow takes 0 1 and 2
overall_diffusion_fraction is 1 or 0.1 This alters GJC throughout the model.
san_upstroke is 0.05 (control) or 0.2 (low) in the SAN only.
*/
#define control_short_ap           0    // Int.  0 is control (long) atrium AP, 1 is short atrium AP (both ISO and ACH reduced APD in Fedorov's papers.
#define san_control_fast_slow      0    // Int.  0 = Control SAN pacing rate, 1 = fast pacing rate (ISO), 2 = slow pacing rate (ACH).
#define overall_diffusion_fraction 1.00 // Real. this takes values 1 or 0.1
#define san_upstroke               0.05 // Real. this is 0.05 for Control SAN upstroke, 0.20 for less SAN upstroke.

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ        */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ        */
#define RTOL        RCONST(1.0e-12)   	  /* scalar relative tolerance            */
#define ATOL        RCONST(1.0e-6)        /* scalar absolute tolerance components */
#define MAXSTEPS    5000
#define ZERO        RCONST(0.0)

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ        */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ        */

#define NEQ          3                    /* number of reaction variables. This is FK3                              */
#define NEQ_SAN	     2                    // FitzHugh Nagumo.
#define NUMPARAMS    13                   /* ODE parameters. Often or all the time I want to specify these as input */

/* Problem Constants */
#define DELTAT       0.100                /* time step. The 2D wants 0.075 ms.                              */
#define amp          0.600                /* stimulation amplitude                  */
#define duration     1.000                /* stimulus duration                      */
// 1200 is less than paranodal pacing, 750 to 840 is between paranodal and SAN, less than 750 is less than SAN.
#define pcl          200.0                /* period of pacing                       */

// #define diffusionDX  overall_diffusion_fraction*0.3548 // 1.064*5.0*0.0650               /* internode diffusion                    */
// #define diffusionDY  overall_diffusion_fraction*0.3548 // 1.064*5.0*0.0650               /* internode diffusion                    */
// #define diffusionDZ  overall_diffusion_fraction*0.6916 // 1.064*5.0*0.1300               /* internode diffusion                    */

#define diffusionDX  overall_diffusion_fraction*0.0650 // 1.064*5.0*0.0650               /* internode diffusion                    */
#define diffusionDY  overall_diffusion_fraction*0.0650 // 1.064*5.0*0.0650               /* internode diffusion                    */
#define diffusionDZ  overall_diffusion_fraction*0.1300 // 1.064*5.0*0.1300               /* internode diffusion                    */

#define DX           0.250                /* X internode spacing                    */
#define DY           0.250                /* Y internode spacing                    */
#define DZ           0.500                /* Z internode spacing: This may not be 0.25 mm in the Chandler model. It may be 0.5 mm. See suppl */
#define usr_MX 128
#define usr_MY 128 // this needs to be read from the geometry data.
#define usr_MZ 60   // the z dimension is small so that my laptop pc can handle it.
/* The stencil is u0 = x,y,z; , u1 = y+1, u2 = x+1, u3 = y-1, u4 = x-1, u5 = z+1, u6 = z - 1 */
#define NBS    7    // 3D arrays have 7 units. Itself, and 6 surrounding it in a standard 1st order FD stencil.
#define DATA   3 // this model may have just 1 bit of data (q value), but I will keep it uniform with hSAN2d.c

/* User-defined data structures and routines           */
/* AppCtx: used by FormIFunction() and FormIJacobian() */
typedef struct {
  DM        da; // DM instance
/* All model data such as geomtery, fibers, cell type data, extracellular distributions must go through the context. */
  PetscInt  ****geometry; // 3D models have 4D array geometry. 0 = itself, 1 is y+1, 2 is x+1, 3 is y-1, 4 is x-1, 5 is z+1, 6 is z-1.
  PetscReal    ***dist; // distance measure.
  PetscReal   ****data;
} AppCtx;

/* sundials context. */
// atrial
typedef struct {
	realtype IV;             // stimulus
	realtype lap_val;        // the laplacian.
	realtype p[NUMPARAMS];
	realtype dvdt;
	realtype factor;
	realtype q; // need this to make SAN heterogeneous.
} *UserData;

// SAN
typedef struct {
	realtype dvdt;
} *UserData_SAN;

// petsc functions.
extern PetscInt FormIFunction(TS,PetscReal,Vec,Vec,Vec,void*);

int getid(int col, int y, int x){
	return (col*y + x);
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

/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
int main(int argc,char **argv)
{
/* Program reduced to use of FD with colouring. */
  TS             ts;                    /* nonlinear solver                                      */
  Vec            u, r;                  /* solution, residual vectors. You need the r for SNES   */
  Mat            J,Jmf = NULL;          /* Jacobian matrices                                     */
  PetscInt       maxsteps = 100;        /* iterations for convergence                            */
  DM             da;                    /* DM object                                             */
  AppCtx         user;                  /* user-defined work context                             */
  SNES           snes;                  /* scalable nonlinear equations solvers (snes)           */
  PetscInt       mybase_x, mybase_y, mybase_z, mysize_x, mysize_y, mysize_z, usr_i, usr_j, usr_k, usr_nbs, time_int;
  PetscInt       file_Counter = 0;
  PetscScalar    ***u_localptr, ***u_localptr2;

// MPI declarations.
  PetscInt    size = 0;                   /* Petsc/MPI                                           */
  PetscInt    rank = 0;

double distance;

/* Standard declarations                                                                          */
/* Sundials out of time loop declarations. All my variables start with usr_···                    */
FILE       *geometry, *spiralFile, *random;
realtype    usr_stimtime, usr_time, usr_t, usr_tout                          ;
double   ***usr_u0_loc, ***usr_u1_loc, ***usr_u2_loc                         ; // the SAN model uses u0 and u1 only.
int         intx, inty, intz, intu0, intu1, intu2, intu3, intu4, intu5, intu6, int_count;
double radius;

double tmp_u0, tmp_u1, tmp_u2;
int int1, int2, int3;

char str[1000];  /* Have several chars, so that we do not keep creating and destroying the str array */

  PetscInitialize(&argc,&argv,(char*)0,help);
  /* Initialize user application context                                                             */
  user.da           = NULL;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

/* start time loop                                                                                    */
time_int = 0; usr_time = 0.0; usr_stimtime = 0.0;
int total_time_int = 200000; // I need 20 s for the paranodal to show its trick.
for(time_int = 0; time_int < total_time_int; time_int++){
   DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,usr_MX,usr_MY,usr_MZ, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,NULL,&da);  /* this decides that u is non-ghosted 1st order star stencil vector. */
   user.da = da;

   DMCreateGlobalVector(da,&u); 
   VecDuplicate(u,&r);

		/* size_y is the width of the array. Do this once. */
		if(time_int < 1){
			   DMDAGetCorners(da,&mybase_x,&mybase_y,&mybase_z,&mysize_x,&mysize_y,&mysize_z); 

		usr_u0_loc = (double ***) calloc(mysize_z, sizeof(double**));
		usr_u1_loc = (double ***) calloc(mysize_z, sizeof(double**));
		usr_u2_loc = (double ***) calloc(mysize_z, sizeof(double**));
		user.dist     = (PetscReal ***)  calloc(usr_MZ, sizeof(PetscReal** ));
	for (usr_k = 0; usr_k < mysize_z; usr_k++){ 
				usr_u0_loc[usr_k] = (double **) calloc(mysize_y, sizeof(double*));
				usr_u1_loc[usr_k] = (double **) calloc(mysize_y, sizeof(double*));
				usr_u2_loc[usr_k] = (double **) calloc(mysize_y, sizeof(double*));
			user.dist[usr_k]     = (PetscReal**)  calloc(usr_MY,sizeof(PetscReal* ));
				/* the usr_mysize_ and usr_mybase_ arrays are not the same on all processors. */
					for (usr_j = 0; usr_j < mysize_y; usr_j++){ 
						usr_u0_loc[usr_k][usr_j] = (double*) calloc(mysize_x,sizeof(double));
						usr_u1_loc[usr_k][usr_j] = (double*) calloc(mysize_x,sizeof(double));
						usr_u2_loc[usr_k][usr_j] = (double*) calloc(mysize_x,sizeof(double));
				user.dist[usr_k][usr_j]     = (PetscReal* ) calloc(usr_MX,sizeof(PetscReal));
					}
	}

			/* Declare the geometry memory. This declares memory on each proc. */

	user.geometry = (PetscInt ****) calloc(mysize_z, sizeof(PetscInt***));
	for (usr_k = 0; usr_k < mysize_z; usr_k++){
			user.geometry[usr_k] = (PetscInt ***) calloc(mysize_y, sizeof(PetscInt**));
			for (usr_j = 0; usr_j < mysize_y; usr_j++){
				user.geometry[usr_k][usr_j] = (PetscInt**) calloc(mysize_x,sizeof(PetscInt*));
				for (usr_i = 0; usr_i < mysize_x; usr_i++)
					user.geometry[usr_k][usr_j][usr_i] = (PetscInt*) calloc(NBS,sizeof(PetscInt));				
			}
	}

			/* Declare the data memory. This declares memory on each proc. */
			user.data = (PetscReal ****) calloc(mysize_z, sizeof(PetscReal***));
			for (usr_k = 0; usr_k < mysize_z; usr_k++){
			user.data[usr_k] = (PetscReal ***) calloc(mysize_y, sizeof(PetscReal**));
			for (usr_j = 0; usr_j < mysize_y; usr_j++){
				user.data[usr_k][usr_j] = (PetscReal**) calloc(mysize_x,sizeof(PetscReal*));
				for (usr_i = 0; usr_i < mysize_x; usr_i++)
					user.data[usr_k][usr_j][usr_i] = (PetscReal*) calloc(DATA, sizeof(PetscReal));				
			}
		}

/* Initialising the random number generator is important. */
random = fopen("/dev/urandom","rb");
fread(&somerandomnumber, sizeof(unsigned int), 1, random);
fclose(random);
init_genrand(somerandomnumber);

// SAN is 1 and 2. It has q != 1 values between 1 and 5.
/* This reads the ASCII geometry on each proc. The file is opened for reading on each proc.                     */
			geometry = fopen("humanSAN.geom","r");
			while(fscanf(geometry,"%d %d %d %d %d %d %d %d %d %d %lf",&intz, &inty, &intx, &intu0, &intu1, &intu2, &intu3, &intu4, &intu5, &intu6, &distance)!=EOF){
				usr_i = (PetscInt)intx - mybase_x; usr_j = (PetscInt)inty - mybase_y;
				usr_k = (PetscInt)intz - mybase_z;
				if(usr_k>=0&&usr_k<mysize_z&&usr_j>=0&&usr_j<mysize_y&&usr_i>=0&&usr_i<mysize_x){
				user.geometry[usr_k][usr_j][usr_i][0] = (PetscInt)intu0;
				user.geometry[usr_k][usr_j][usr_i][1] = (PetscInt)intu1;
				user.geometry[usr_k][usr_j][usr_i][2] = (PetscInt)intu2;
				user.geometry[usr_k][usr_j][usr_i][3] = (PetscInt)intu3;
				user.geometry[usr_k][usr_j][usr_i][4] = (PetscInt)intu4;
				user.geometry[usr_k][usr_j][usr_i][5] = (PetscInt)intu5;
				user.geometry[usr_k][usr_j][usr_i][6] = (PetscInt)intu6;
				user.dist[usr_k][usr_j][usr_i]        = (PetscReal)distance;

user.data[usr_k][usr_j][usr_i][0] = -1.0;
user.data[usr_k][usr_j][usr_i][1] = -1.0; // just make sure that your data is all there.
user.data[usr_k][usr_j][usr_i][2] = -1.0;



						user.data[usr_k][usr_j][usr_i][1]     = 1.0;

if(user.geometry[usr_k][usr_j][usr_i][0]==1||user.geometry[usr_k][usr_j][usr_i][0]==2){
	if(san_control_fast_slow==0) user.data[usr_k][usr_j][usr_i][1] = (7.00 + (1.0 - 0.5*genrand_real1())*2.0); // Control.
	if(san_control_fast_slow==1) user.data[usr_k][usr_j][usr_i][1] = (5.00 + (1.0 - 0.5*genrand_real1())*2.0); // fast.
	if(san_control_fast_slow==2) user.data[usr_k][usr_j][usr_i][1] = (20.0 + (1.0 - 0.5*genrand_real1())*2.0); // slow.
}

if(user.geometry[usr_k][usr_j][usr_i][0]<0||user.geometry[usr_k][usr_j][usr_i][0]>6){ printf("blah\n"); exit(1);}

// fix the geometry to one of your 5 conditions.
// case 1. SAN   only         . In this case, set paranodal to atrial, border to atrial.
// case 2. SAN + border       . In this case, set paranodal to atrial.
// case 3. SAN + para         . In this case set border to atrial.
// case 4. para only
// case 5. all together.
if(atoi(argv[1])==1){ // SAN only.
// rid of paranodal
	for(int_count = 0; int_count < 7; int_count++){
	if(user.geometry[usr_k][usr_j][usr_i][int_count]==3) user.geometry[usr_k][usr_j][usr_i][int_count] = 4;
	}
// rid of border.
	for(int_count = 0; int_count < 7; int_count++){
	if(user.geometry[usr_k][usr_j][usr_i][int_count]==6) user.geometry[usr_k][usr_j][usr_i][int_count] = 4;
	}
}
if(atoi(argv[1])==2){ // SAN+border
// rid of paranodal
	for(int_count = 0; int_count < 7; int_count++){
	if(user.geometry[usr_k][usr_j][usr_i][int_count]==3) user.geometry[usr_k][usr_j][usr_i][int_count] = 4;
	}
}
if(atoi(argv[1])==3){ // SAN + para
// rid of border.
	for(int_count = 0; int_count < 7; int_count++){
	if(user.geometry[usr_k][usr_j][usr_i][int_count]==6) user.geometry[usr_k][usr_j][usr_i][int_count] = 4;
	}
}
if(atoi(argv[1])==4){ // para only
// rid of border and SAN.
//printf("I have decided not to do the para only case. Exiting.\n");
//exit(1);
	for(int_count = 0; int_count < 7; int_count++){
	if(user.geometry[usr_k][usr_j][usr_i][int_count]==1) user.geometry[usr_k][usr_j][usr_i][int_count] = 4;
	if(user.geometry[usr_k][usr_j][usr_i][int_count]==2) user.geometry[usr_k][usr_j][usr_i][int_count] = 4;
	if(user.geometry[usr_k][usr_j][usr_i][int_count]==6) user.geometry[usr_k][usr_j][usr_i][int_count] = 4;
	}
}
if(atoi(argv[1])==5){} // SAN + para + border. do nothing.

if(atoi(argv[1])<1||atoi(argv[1])>5){ printf("Do not know this case. Exiting at line 507.\n"); exit(1);}
				} 
			}

// output the q values
		for(usr_k = 0; usr_k < mysize_z; usr_k++)
		for(usr_j = 0; usr_j < mysize_y; usr_j++)
		for(usr_i = 0; usr_i < mysize_x; usr_i++)
			usr_u0_loc[usr_k][usr_j][usr_i] = user.data[usr_k][usr_j][usr_i][1];
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_k = 0; usr_k < mysize_z; usr_k++)
	for(usr_j = 0; usr_j < mysize_y; usr_j++)
	for(usr_i = 0; usr_i < mysize_x; usr_i++)
		u_localptr[mybase_z + usr_k][mybase_y + usr_j][mybase_x + usr_i] = usr_u0_loc[usr_k][usr_j][usr_i];
    DMDAVecRestoreArray(da,u,&u_localptr);
            sprintf(str,"q_values3D.bin");
			PetscViewer viewer4;
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_WRITE,&viewer4);
			VecView(u,viewer4);
			PetscViewerDestroy(&viewer4);

	/* Sundials initial conditions                                                             */
		for(usr_k = 0; usr_k < mysize_z; usr_k++)
		for(usr_j = 0; usr_j < mysize_y; usr_j++)
		for(usr_i = 0; usr_i < mysize_x; usr_i++){
			if(user.geometry[usr_k][usr_j][usr_i][0]>0){ // all 5 cell types.
				usr_u0_loc[usr_k][usr_j][usr_i] = 0.0; 
				usr_u1_loc[usr_k][usr_j][usr_i] =  0.5; usr_u2_loc[usr_k][usr_j][usr_i] = 0.0;
			}
			else
			{	
				/* impossible value. For computations and paraview visualisation. */
				usr_u0_loc[usr_k][usr_j][usr_i] = usr_u1_loc[usr_k][usr_j][usr_i] = usr_u2_loc[usr_k][usr_j][usr_i] = -0.1;
			}
		}

// random initial conditions. bad idea.
/*
		for(usr_k = 0; usr_k < mysize_z; usr_k++)
		for(usr_j = 0; usr_j < mysize_y; usr_j++)
		for(usr_i = 0; usr_i < mysize_x; usr_i++){
			if(user.geometry[usr_k][usr_j][usr_i][0]>0){ // all 5 cell types.
				usr_u0_loc[usr_k][usr_j][usr_i] = genrand_real1(); 
				usr_u1_loc[usr_k][usr_j][usr_i] = genrand_real1(); usr_u2_loc[usr_k][usr_j][usr_i] = genrand_real1();
			}
			else
			{	
				usr_u0_loc[usr_k][usr_j][usr_i] = usr_u1_loc[usr_k][usr_j][usr_i] = usr_u2_loc[usr_k][usr_j][usr_i] = -0.1;
			}
		}
*/

// put in your spiral wave here.
/*
spiralFile = fopen("spiralInit3D.data","r"); // format: x z u0 u1 u2. So this is for every y.
while(fscanf(spiralFile,"%d %d %lf %lf %lf",&int1, &int2, &tmp_u0, &tmp_u1, &tmp_u2)!=EOF){
   usr_k = (PetscInt)int1 - mybase_z; usr_i = (PetscInt)int2  - mybase_x;
				if(usr_k>=0&&usr_k<mysize_z&&usr_i>=0&&usr_i<mysize_x){
					for(usr_j = 0; usr_j < mysize_y; usr_j++){
int3 = usr_j + mybase_y; // actual Y coordinate of the location. I want the spiral between plus minus 5 of y = 48
// if(int3>(52-3)&&int3<(52+3)){
//						if(user.geometry[usr_k][usr_j][usr_i][0]>3&&user.geometry[usr_k][usr_j][usr_i][0]<6){ // macro-rentry
						if(user.geometry[usr_k][usr_j][usr_i][0]>0&&user.geometry[usr_k][usr_j][usr_i][0]<3){ // micro-rentry
					 	usr_u0_loc[usr_k][usr_j][usr_i] = (double)tmp_u0;
						usr_u1_loc[usr_k][usr_j][usr_i] = (double)tmp_u1;
						usr_u2_loc[usr_k][usr_j][usr_i] = (double)tmp_u2;
						}
// 					   }
					}
				}
}
fclose(spiralFile);
*/

// output the initial conditions.
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_k = 0; usr_k < mysize_z; usr_k++)
	for(usr_j = 0; usr_j < mysize_y; usr_j++)
	for(usr_i = 0; usr_i < mysize_x; usr_i++)
		u_localptr[mybase_z + usr_k][mybase_y + usr_j][mybase_x + usr_i] = usr_u0_loc[usr_k][usr_j][usr_i];
    DMDAVecRestoreArray(da,u,&u_localptr);
            sprintf(str,"initialConditions3D.bin");
			PetscViewer viewer5;
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_WRITE,&viewer5);
			VecView(u,viewer5);
			PetscViewerDestroy(&viewer5);

} // end of time_int < 1.


// part I. Reaction solution.
	for(usr_k = 0; usr_k < mysize_z; usr_k++) // these usr_k, usr_j and usr_i are for local locations.
	for(usr_j = 0; usr_j < mysize_y; usr_j++) 
	for(usr_i = 0; usr_i < mysize_x; usr_i++){


		if(user.geometry[usr_k][usr_j][usr_i][0]>0&&user.geometry[usr_k][usr_j][usr_i][0]<6){
	/*******************************************************************************/

		    N_Vector y          ;                    /* The state vector                 */
		    void    *cvode_mem  ;               /* Sundials memory pointer          */
		    UserData data       ;                 /* Sundials data structure instance */ 
		    /* Create serial vector of length NEQ for I.C. and abstol          */
		    cvode_mem     = NULL;
		    data          = NULL;
		    y             = NULL;
		    y             = N_VNew_Serial(NEQ);
		    data          = (UserData) malloc(sizeof *data);
		    data->factor = 1.0; /* user specified Sundials data.                */
		    data->IV     = 0.0;
			data->q = user.data[usr_k][usr_j][usr_i][1]; // the q is needed only for tissues 1 and 2.

// 0.01 120 45 16.5 16.5 1000 75 1.5 100 0.2 0.5 0.01 2.5 0.5 2.0
if(user.geometry[usr_k][usr_j][usr_i][0]==1||user.geometry[usr_k][usr_j][usr_i][0]==2){ // SAN. SAN in the geometry is tissue type 1 and 2.

//  0.05 130.0 50.0 12.50 10.0 18.2 18.2 1.5 600.0 0.1 0.5 0.01 5.0 2
		data->p[0] /* double tau_d     */  =    san_upstroke; // 0.050 ; // one of my parameters.
		data->p[1] /* double tau_r     */  =    130.0 ;
		data->p[2] /* double tau_si    */  =    50.00 ;
		data->p[3] /* double tau_0     */  =    12.5  ;
		data->p[4] /* double tau_vp    */  =    10.0  ;
		data->p[5] /* double tau_vm1   */  =    18.2  ;
		data->p[6] /* double tau_vm2   */  =    18.2  ;
		data->p[7] /* double tau_wp    */  =    1.5   ;

// this needs to be revised carefully for pacing at 850. I got it, do a bit more checking.
		data->p[8] /* double tau_wm    */  =    100.0*data->q; // control is 600, fast is 200, slow is 1200 (max. vales).

		data->p[9] /* double u_c       */  =    0.1   ;
		data->p[10] /* double u_v       */ =    0.5   ; // does not matter because tau_vm1 and tau_vm2 are equal.
		data->p[11] /* double usi_c     */ =    0.01  ;
		data->p[12] /* double k         */ =    5.0   ;

}else if(user.geometry[usr_k][usr_j][usr_i][0]==3){ // paranodal volume.

		data->p[0] /* double tau_d     */  =    0.05; // 0.175 ; // lower upstroke than SAN.
		data->p[1] /* double tau_r     */  =    130.0 ;
		data->p[2] /* double tau_si    */  =    50.00 ;
		data->p[3] /* double tau_0     */  =    12.5  ;
		data->p[4] /* double tau_vp    */  =    10.0  ;
		data->p[5] /* double tau_vm1   */  =    18.2  ;
		data->p[6] /* double tau_vm2   */  =    18.2  ;
		data->p[7] /* double tau_wp    */  =    1.5   ;
		data->p[8] /* double tau_wm    */  =    800.0; // 1400.0; // longer cycle length than SAN, but shorter than the ACH case.
		data->p[9] /* double u_c       */  =    0.1   ;
		data->p[10] /* double u_v       */ =    0.5   ; // does not matter because tau_vm1 and tau_vm2 are equal.
		data->p[11] /* double usi_c     */ =    0.01  ;
		data->p[12] /* double k         */ =    5.0   ;

} else if(user.geometry[usr_k][usr_j][usr_i][0]>3){ // atrial part (adipose in Chandler paper is treated as atrial tissue here).
	// Control atrial AP
	if(control_short_ap==0){
			data->p[0] /* double tau_d     */  =    0.005; // 0.2000;
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
}
else if(user.geometry[usr_k][usr_j][usr_i][0]==7){ // Apoptosis: Use this to make the SAN-atrial junctions inactive.
// 20.0  40.0 310.0 12.5 10.0 10.0 30.0  1.0 552.0 0.95 0.5 0.1 10.0 1
		data->p[0] /* double tau_d     */  =    20.000; // low upstroke.
		data->p[1] /* double tau_r     */  =    40.0  ;
		data->p[2] /* double tau_si    */  =    310.0 ;
		data->p[3] /* double tau_0     */  =    12.5  ;
		data->p[4] /* double tau_vp    */  =    10.0  ;
		data->p[5] /* double tau_vm1   */  =    10.0  ;
		data->p[6] /* double tau_vm2   */  =    30.0  ;
		data->p[7] /* double tau_wp    */  =    1.00  ;
		data->p[8] /* double tau_wm    */  =    352.0 ;
		data->p[9] /* double u_c       */  =    0.5   ;
		data->p[10] /* double u_v       */ =    0.85  ; // does not matter because tau_vm1 and tau_vm2 are equal.
		data->p[11] /* double usi_c     */ =    1.00  ;
		data->p[12] /* double k         */ =    10.0  ;
}
		    usr_t = 0; usr_tout = DELTAT;
		    Ith(y,1) = usr_u0_loc[usr_k][usr_j][usr_i];
		    Ith(y,2) = usr_u1_loc[usr_k][usr_j][usr_i];
		    Ith(y,3) = usr_u2_loc[usr_k][usr_j][usr_i];
		    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
		    CVodeInit(cvode_mem, f, 0.0, y);
		    CVodeWFtolerances(cvode_mem, ewt);
		    CVodeSetUserData(cvode_mem, data);
		    CVodeSStolerances(cvode_mem, RTOL, ATOL);
		    CVDense(cvode_mem, NEQ);
		    CVodeSetMaxStep(cvode_mem, DELTAT);
		    CVode(cvode_mem, usr_tout, y, &usr_t, CV_NORMAL);

		    usr_u0_loc[usr_k][usr_j][usr_i] = Ith(y,1);
		    usr_u1_loc[usr_k][usr_j][usr_i] = Ith(y,2);
		    usr_u2_loc[usr_k][usr_j][usr_i] = Ith(y,3);

		    N_VDestroy_Serial(y);
		    CVodeFree(&cvode_mem);
		    free(data);

	/*******************************************************************************/
		} /* end of geometry if.                                                        */
} // end of Sundials loop.

// petsc part II
   TSCreate(PETSC_COMM_WORLD,&ts); 
   TSSetProblemType(ts,TS_NONLINEAR); 
   TSSetType(ts,TSBEULER); 
   TSSetDM(ts,da); 
   TSSetIFunction(ts,r,FormIFunction,&user); 
   TSSetDuration(ts,maxsteps,DELTAT); // duration of simulation.
/********** put your usr_u0_loc into u here. ********************************/
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_k = 0; usr_k < mysize_z; usr_k++)
	for(usr_j = 0; usr_j < mysize_y; usr_j++)
	for(usr_i = 0; usr_i < mysize_x; usr_i++)
		u_localptr[mybase_z + usr_k][mybase_y + usr_j][mybase_x + usr_i] = usr_u0_loc[usr_k][usr_j][usr_i];
    DMDAVecRestoreArray(da,u,&u_localptr);

/****************************************************************************/
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
	for(usr_k = 0; usr_k < mysize_z; usr_k++)
	for(usr_j = 0; usr_j < mysize_y; usr_j++)
	for(usr_i = 0; usr_i < mysize_x; usr_i++)
		 usr_u0_loc[usr_k][usr_j][usr_i] = u_localptr[mybase_z + usr_k][mybase_y + usr_j][mybase_x + usr_i];
    DMDAVecRestoreArray(da,u,&u_localptr);
/****************************************************************************/
// do some VTK/binary output.
// example in snes/ex61view.c. This example is serial, I need parallel file i/o.
if(time_int%10==0 &&(time_int>(int)(0.0*total_time_int) ) ){ // do not change this output rate: you made decisions based on this output (file number)
            sprintf(str,"my_3d%d.bin",file_Counter++);
			PetscViewer viewer2;
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_WRITE,&viewer2);
			VecView(u,viewer2);
			PetscViewerDestroy(&viewer2);
}

// this is in your time loop for now.
   VecDestroy(&u); 
   VecDestroy(&r); 
   MatDestroy(&J); 
   MatDestroy(&Jmf); 
   TSDestroy(&ts); 
   DMDestroy(&da); 

usr_time = (PetscReal)time_int*DELTAT; // actual time.

} // end of time.

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   free(usr_u0_loc);
   free(usr_u1_loc);
   free(usr_u2_loc);

   PetscFinalize();
   PetscFunctionReturn(0);
}

/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/

/* FormIFunction = Udot - RHSFunction */
PetscInt FormIFunction(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx)
{
  AppCtx         *user=(AppCtx*)ctx;
  DM             da   = (DM)user->da;
  PetscInt       i,j,k,xs,ys,zs,xm,ym,zm,gj, gi, gk;
  PetscReal      hx,hy,hz,sx,sy,sz;
  PetscReal      U0, U1, U2, U3, U4, U5, U6;
  PetscScalar    uxx,uyy,uzz, ***uarray, ***f, ***udot;
  PetscScalar    ux, uy, uz;
  Vec            localU;

PetscReal dd;

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
   DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 
/*
This calculation of the RHS, or udot - RHS as they call it is
wrong. I started with ex15.c in the ts examples. So fix this.
Numerically, it is stable. Still some work to be done.
*/
  /* Compute function over the locally owned part of the grid */
  for (k=zs; k<zs+zm; k++)
  for (j=ys; j<ys+ym; j++)
    for (i=xs; i<xs+xm; i++)
	f[k][j][i] = 0.0; // so that i do not have to deal with empty space.

    for (k=zs; k<zs+zm; k++)
    for (j=ys; j<ys+ym; j++) 
    for (i=xs; i<xs+xm; i++) {

// at these values if i, j, k, geometry is definitely 0. Also, calculating first order terms here is not ok.
//	if( (k==0)||(k==usr_MZ-1)||(j==0)||(j==usr_MY-1)||(i==0)||(i==usr_MX-1) ) continue;
		

// the gk go from start of box to end of box in each direction.
gk = k - zs;
gj = j - ys;
gi = i - xs;

if(user->geometry[gk][gj][gi][0]==1||user->geometry[gk][gj][gi][0]==2){ // SAN
  dd = user->dist[gk][gj][gi];
  hx = DX; sx = (0.0000325 + dd*dd*diffusionDX)/(hx*hx); // same as the 2D now.
  hy = DY; sy = (0.0000325 + dd*dd*diffusionDY)/(hy*hy);
  hz = DZ; sz = (0.0000325 + dd*dd*diffusionDZ)/(hz*hz);
}
else if(user->geometry[gk][gj][gi][0]==3){ // paranodal oscillator.
  hx = DX; sx = 0.5*diffusionDX/(hx*hx);
  hy = DY; sy = 0.5*diffusionDY/(hy*hy);
  hz = DZ; sz = 0.5*diffusionDZ/(hz*hz);
}
else if(user->geometry[gk][gj][gi][0]==6){ // insulating border.
  hx = DX; sx = 0.0*diffusionDX/(hx*hx);
  hy = DY; sy = 0.0*diffusionDY/(hy*hy);
  hz = DZ; sz = 0.0*diffusionDZ/(hz*hz);
}else{
  hx = DX; sx = diffusionDX/(hx*hx); // atrium and adipose tissue.
  hy = DY; sy = diffusionDY/(hy*hy);
  hz = DZ; sz = diffusionDZ/(hz*hz);
}
		/* This is still hacky. I am not happy with the way boundaries are broken down into separate components. */
		if(user->geometry[gk][gj][gi][0]>0){
			if(i>0&&j>0&&k>0&&i<usr_MX-1&&j<usr_MY-1&&k<usr_MZ-1){
			                         U1 = uarray[k][j+1][i]; 
		     U4 = uarray[k][j][i-1]; U0 = uarray[k][j  ][i]; U2 = uarray[k][j][i+1]; 
		                             U3 = uarray[k][j-1][i];

								     U5 = uarray[k+1][j][i]; 
								     U6 = uarray[k-1][j][i]; 
// you have not considered the cases 1 & 2 == 0, 2 & 3 == 0, 3 & 4 == 0 , 4 & 1 ==0
// the 3 neighbours being 0 cases are already in the code.
		if(user->geometry[gk][gj][gi][1]  ==0) 				             U1 = U3;
		if(user->geometry[gk][gj][gi][3]  ==0) 				             U3 = U1;
		if((user->geometry[gk][gj][gi][1] ==0)&&(user->geometry[gk][gj][gi][3]==0) )  U3 = U1 = U0;
		if(user->geometry[gk][gj][gi][2]  ==0) 				             U2 = U4;
		if(user->geometry[gk][gj][gi][4]  ==0) 				             U4 = U2;
		if((user->geometry[gk][gj][gi][2] ==0)&&(user->geometry[gk][gj][gi][4]==0) )  U4 = U2 = U0;
		if(user->geometry[gk][gj][gi][5]  ==0) 				             U5 = U6;
		if(user->geometry[gk][gj][gi][6]  ==0) 				             U6 = U5;
		if((user->geometry[gk][gj][gi][5] ==0)&&(user->geometry[gk][gj][gi][5]==0) ) U6 = U5 = U0;

		if(user->geometry[gk][gj][gi][1]  ==6) 				             U1 = U3;
		if(user->geometry[gk][gj][gi][3]  ==6) 				             U3 = U1;
		if((user->geometry[gk][gj][gi][1] ==6)&&(user->geometry[gk][gj][gi][3]==6) )  U3 = U1 = U0;
		if(user->geometry[gk][gj][gi][2]  ==6) 				             U2 = U4;
		if(user->geometry[gk][gj][gi][4]  ==6) 				             U4 = U2;
		if((user->geometry[gk][gj][gi][2] ==6)&&(user->geometry[gk][gj][gi][4]==6) )  U4 = U2 = U0;
		if(user->geometry[gk][gj][gi][5]  ==6) 				             U5 = U6;
		if(user->geometry[gk][gj][gi][6]  ==6) 				             U6 = U5;
		if((user->geometry[gk][gj][gi][5] ==6)&&(user->geometry[gk][gj][gi][5]==6) ) U6 = U5 = U0;

		if((user->geometry[gk][gj][gi][1]==0)&&(user->geometry[gk][gj][gi][3]==6) )  U3 = U1 = U0;
		if((user->geometry[gk][gj][gi][1]==6)&&(user->geometry[gk][gj][gi][3]==0) )  U3 = U1 = U0;
		if((user->geometry[gk][gj][gi][2]==0)&&(user->geometry[gk][gj][gi][4]==6) )  U4 = U2 = U0;
		if((user->geometry[gk][gj][gi][2]==6)&&(user->geometry[gk][gj][gi][4]==0) )  U4 = U2 = U0;
		if((user->geometry[gk][gj][gi][5]==0)&&(user->geometry[gk][gj][gi][5]==6) )  U6 = U5 = U0;
		if((user->geometry[gk][gj][gi][5]==6)&&(user->geometry[gk][gj][gi][5]==0) )  U6 = U5 = U0;

		if(user->geometry[gk][gj][gi][0]  ==0) U1 = U2 = U3 = U4 = U5 = U6 = -0.1;
		if(user->geometry[gk][gj][gi][0]  ==6) U1 = U2 = U3 = U4 = U5 = U6 = -0.1;

			uxx = (-2.0*U0 + U4 + U2);
			uyy = (-2.0*U0 + U3 + U1);
			uzz = (-2.0*U0 + U5 + U6);
			ux = (U2 - U4)/(2.0*DX);
			uy = (U1 - U3)/(2.0*DY);
			uz = (U5 - U6)/(2.0*DX);
		if((user->geometry[gk][gj][gi][1]==0)||(user->geometry[gk][gj][gi][2]==0)||(user->geometry[gk][gj][gi][3]==0)||
		    (user->geometry[gk][gj][gi][4]==0)||(user->geometry[gk][gj][gi][5]==0)||(user->geometry[gk][gj][gi][6]==0)){
			ux = 0.0; uy = 0.0; uz = 0.0;
		}
		f[k][j][i] = udot[k][j][i] - (uxx*sx + uyy*sy + uzz*sz /* + 0.01*ux + 0.01*uy + 0.01*uz */ ); // the udot is not zero
            } // end of if statement to keep it inside the box.
	} // end of if geometry statement.
    } // end of for loops.

  /* Restore vectors */
   DMDAVecRestoreArray(da,localU,&uarray); 

   DMDAVecRestoreArray(da,F,&f);
   DMDAVecRestoreArray(da,Udot,&udot); 
   DMRestoreLocalVector(da,&localU); 
  PetscFunctionReturn(0);
}

