/*
Sanjay Kharche.
4 July 2015.
Figure 3:
Col 2: ECG
Col 3: AP profiles from SAN centre, SAN periphery, atrium, and paranodal area.
*/
static char help[] = "Mouse 2015 3D SAN measurements, post-processing.\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscvec.h>
#include <petscsys.h>

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
#define RTOL        RCONST(1.0e-12)   	  /* scalar relative tolerance            */
#define ATOL        RCONST(1.0e-6)        /* scalar absolute tolerance components */
#define MAXSTEPS    5000
#define ZERO        RCONST(0.0)

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ        */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ        */

#define NEQ         3                     /* number of reaction variables. This is FK3 */

/* Problem Constants */
#define DELTAT       0.250                /* time step                              */
#define amp          0.500                /* stimulation amplitude                  */
#define duration     1.000                /* stimulus duration                      */

#define NUMPARAMS    14                   /* ODE parameters. Often or all the time I want to specify these as input */

#define pcl          100.0                /* period of pacing                           */

// #define diffusion0    0.00000                /* internode diffusion                  */
// #define diffusion1    0.01000                /* internode diffusion                  */
// #define diffusion2    0.00010                /* internode diffusion                  */

#define DX           0.250                 /* X internode spacing                       */
#define DY           0.250                 /* Y internode spacing                       */
#define DZ	     0.500

#define usr_MX 128
#define usr_MY 128 // this needs to be read from the geometry data.
#define usr_MZ 60   // the z dimension is small so that my laptop pc can handle it.
/* The stencil is u0 = x,y,z; , u1 = y+1, u2 = x+1, u3 = y-1, u4 = x-1, u5 = z+1, u6 = z - 1 */
#define NBS    7    // 3D arrays have 7 units. Itself, and 6 surrounding it in a standard 1st order FD stencil.

/* User-defined data structures and routines           */
/* AppCtx: used by FormIFunction() and FormIJacobian() */
typedef struct {
  DM        da; // DM instance
/* All model data such as geomtery, fibers, cell type data, extracellular distributions must go through the context. */
  PetscInt  ****geometry; // 3D models have 4D array geometry. 0 = itself, 1 is y+1, 2 is x+1, 3 is y-1, 4 is x-1, 5 is z+1, 6 is z-1.
} AppCtx;

		//! Byte swap int
		int32_t swap_int32( int32_t val )
		{
		    val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF ); 
		    return (val << 16) | ((val >> 16) & 0xFFFF);
		}

		float ReverseFloat( const float inFloat )
		{
		   float retVal;
		   char *floatToConvert = ( char* ) & inFloat;
		   char *returnFloat = ( char* ) & retVal;

		   // swap the bytes into a temporary buffer
		   returnFloat[0] = floatToConvert[3];
		   returnFloat[1] = floatToConvert[2];
		   returnFloat[2] = floatToConvert[1];
		   returnFloat[3] = floatToConvert[0];

		   return retVal;
		}


int main(int argc,char **argv)
{

  PetscInt   file_Counter = 0;
  FILE      *geometry, *output;
  char       str[1000]; // so you do not keep creating and destroying this.
float   ***usr_u0_loc, ***usr_u0_loc_old, ***dvdt, ***dvdt2;
float ***activationTime;
int         intx, inty, intz, intu0,  intu1 , intu2 , intu3 , intu4 , intu5 , intu6 , intu7 , intu8 ;
int                           intu9,  intu10, intu11, intu12, intu13, intu14, intu15, intu16, intu17;
int                           intu18, intu19, intu20, intu21, intu22, intu23, intu24, intu25, intu26;
double distance;
float ecg;
int rx, ry, rz, xl, yl, zl;
float r, r3;

/* Program reduced to use of FD with colouring. */
  Vec            u, u_old;              /* solution, residual vectors. You need the r for SNES */
  AppCtx         user;                  /* user-defined work context  */
  DM             da;
  PetscInt       usr_i, usr_j, usr_k, time_int;
  PetscScalar    ***u_localptr;
PetscInt start, end, leadingPacemakerLocation, lp_x, lp_y, lp_z, file_Counter_start;

char fileNames[100];

	double cl_value = 0;

  PetscInitialize(&argc,&argv,(char*)0,help);
  // allocate your measurement vectors.
		usr_u0_loc     = (float ***) calloc(usr_MZ, sizeof(float**));
		usr_u0_loc_old = (float ***) calloc(usr_MZ, sizeof(float**));
		dvdt           = (float ***) calloc(usr_MZ, sizeof(float**));
		dvdt2          = (float ***) calloc(usr_MZ, sizeof(float**));
		activationTime = (float ***) calloc(usr_MZ, sizeof(float**));
	for (usr_k = 0; usr_k < usr_MZ; usr_k++){ 
				usr_u0_loc[usr_k]     = (float **) calloc(usr_MY, sizeof(float*));
				usr_u0_loc_old[usr_k] = (float **) calloc(usr_MY, sizeof(float*));
				dvdt[usr_k]           = (float **) calloc(usr_MY, sizeof(float*));
				dvdt2[usr_k]          = (float **) calloc(usr_MY, sizeof(float*));
				activationTime[usr_k] = (float **) calloc(usr_MY, sizeof(float*));
				/* the usr_mysize_ and usr_mybase_ arrays are not the same on all processors. */
					for (usr_j = 0; usr_j < usr_MY; usr_j++){ 
						usr_u0_loc[usr_k][usr_j]     = (float*) calloc(usr_MX,sizeof(float));
						usr_u0_loc_old[usr_k][usr_j] = (float*) calloc(usr_MX,sizeof(float));
						dvdt[usr_k][usr_j]           = (float*) calloc(usr_MX,sizeof(float));
						dvdt2[usr_k][usr_j]          = (float*) calloc(usr_MX,sizeof(float));
						activationTime[usr_k][usr_j]          = (float*) calloc(usr_MX,sizeof(float));
					}
	}

			/* Declare the geometry memory. This declares memory on each proc. */
	user.geometry = (PetscInt ****) calloc(usr_MZ, sizeof(PetscInt***));
	for (usr_k = 0; usr_k < usr_MZ; usr_k++){
			user.geometry[usr_k] = (PetscInt ***) calloc(usr_MY, sizeof(PetscInt**));
			for (usr_j = 0; usr_j < usr_MY; usr_j++){
				user.geometry[usr_k][usr_j] = (PetscInt**) calloc(usr_MX,sizeof(PetscInt*));
				for (usr_i = 0; usr_i < usr_MX; usr_i++)
					user.geometry[usr_k][usr_j][usr_i] = (PetscInt*) calloc(NBS,sizeof(PetscInt));				
			}
	}

				/* read in your geometry here.  		*/
				/* This reads the ASCII geometry on each proc. 
				The file is opened for reading on each proc.    */
					geometry = fopen("humanSAN.geom","r");
			while(fscanf(geometry,"%d %d %d %d %d %d %d %d %d %d %lf",&intz, &inty, &intx, &intu0, &intu1, &intu2, &intu3, &intu4, &intu5, &intu6, &distance)!=EOF){
				usr_k = (PetscInt)intz; usr_j = (PetscInt)inty; usr_i = (PetscInt)intx;
				if(usr_j>=0&&usr_j<usr_MY&&usr_i>=0&&usr_i<usr_MX){
				user.geometry[usr_k][usr_j][usr_i][0]  = (PetscInt)intu0 ;
				user.geometry[usr_k][usr_j][usr_i][1]  = (PetscInt)intu1 ;
				user.geometry[usr_k][usr_j][usr_i][2]  = (PetscInt)intu2 ;
				user.geometry[usr_k][usr_j][usr_i][3]  = (PetscInt)intu3 ;
				user.geometry[usr_k][usr_j][usr_i][4]  = (PetscInt)intu4 ;
				user.geometry[usr_k][usr_j][usr_i][5]  = (PetscInt)intu5 ;
				user.geometry[usr_k][usr_j][usr_i][6]  = (PetscInt)intu6 ;
/*
				user.geometry[usr_k][usr_j][usr_i][7]  = (PetscInt)intu7 ;
				user.geometry[usr_k][usr_j][usr_i][8]  = (PetscInt)intu8 ;
				user.geometry[usr_k][usr_j][usr_i][9]  = (PetscInt)intu9 ;
				user.geometry[usr_k][usr_j][usr_i][10] = (PetscInt)intu10;
				user.geometry[usr_k][usr_j][usr_i][11] = (PetscInt)intu11;
				user.geometry[usr_k][usr_j][usr_i][12] = (PetscInt)intu12;
				user.geometry[usr_k][usr_j][usr_i][13] = (PetscInt)intu13;
				user.geometry[usr_k][usr_j][usr_i][14] = (PetscInt)intu14;
				user.geometry[usr_k][usr_j][usr_i][15] = (PetscInt)intu15;
				user.geometry[usr_k][usr_j][usr_i][16] = (PetscInt)intu16;
				user.geometry[usr_k][usr_j][usr_i][17] = (PetscInt)intu17;
				user.geometry[usr_k][usr_j][usr_i][18] = (PetscInt)intu18;
				user.geometry[usr_k][usr_j][usr_i][19] = (PetscInt)intu19;
				user.geometry[usr_k][usr_j][usr_i][20] = (PetscInt)intu20;
				user.geometry[usr_k][usr_j][usr_i][21] = (PetscInt)intu21;
				user.geometry[usr_k][usr_j][usr_i][22] = (PetscInt)intu22;
				user.geometry[usr_k][usr_j][usr_i][23] = (PetscInt)intu23;
				user.geometry[usr_k][usr_j][usr_i][24] = (PetscInt)intu24;
				user.geometry[usr_k][usr_j][usr_i][25] = (PetscInt)intu25;
				user.geometry[usr_k][usr_j][usr_i][26] = (PetscInt)intu26;
*/
				} 
			}

// initialise the activationTime
	for(usr_k = 0; usr_k < usr_MZ; usr_k++)
	for(usr_j = 0; usr_j < usr_MY; usr_j++)
	for(usr_i = 0; usr_i < usr_MX; usr_i++){
		activationTime[usr_k][usr_j][usr_i] = -1.0; // impossible time.
	}

leadingPacemakerLocation = -1; file_Counter_start = -1; ecg = 0.0;

// For AP and ECG, I will use all the data I have
start = 2;
end   = atoi(argv[1]); // count in bash, and pass it to this program.
for(file_Counter = start; file_Counter < end; file_Counter++){ // this is where I saw 1 activation of the SAN model.

  /* Initialize user application context */
  user.da           = NULL;

   DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,usr_MX,usr_MY,usr_MZ, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,NULL,&da);  /* this decides that u is non-ghosted 1st order star/box stencil vector. */
   user.da = da;
   DMDASetUniformCoordinates(da,0,usr_MX, 0, usr_MY, 0, usr_MZ); // this is to set up the dimensions so that we get a good vts file.
   DMCreateGlobalVector(da,&u);      // at this time.

sprintf(fileNames,"gunzip my_3d%d.bin.gz",file_Counter);
system(fileNames);
memset(&fileNames[0],0,sizeof(fileNames));
sprintf(fileNames,"gunzip my_3d%d.bin.gz",file_Counter-1);
system(fileNames);
memset(&fileNames[0],0,sizeof(fileNames));

   // inputs are Petsc binary files.
			PetscViewer viewer_in;
			sprintf(str,"my_3d%d.bin",file_Counter);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_READ,&viewer_in);
			VecLoad(u,viewer_in);
			PetscViewerDestroy(&viewer_in);

   DMCreateGlobalVector(da,&u_old);      // at this time.
   // inputs are Petsc binary files.
			PetscViewer viewer_in_old;
			sprintf(str,"my_3d%d.bin",file_Counter-1);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_READ,&viewer_in_old);
			VecLoad(u_old,viewer_in_old);
			PetscViewerDestroy(&viewer_in_old);

sprintf(fileNames,"gzip my_3d%d.bin",file_Counter);
system(fileNames);
memset(&fileNames[0],0,sizeof(fileNames));
sprintf(fileNames,"gzip my_3d%d.bin",file_Counter-1);
system(fileNames);
memset(&fileNames[0],0,sizeof(fileNames));


//**********************************************************************************************
// get your vectors (i.e. 3D arrays)  for measurements here.
// now.
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_k = 0; usr_k < usr_MZ; usr_k++)
	for(usr_j = 0; usr_j < usr_MY; usr_j++)
	for(usr_i = 0; usr_i < usr_MX; usr_i++)
	 usr_u0_loc[usr_k][usr_j][usr_i] = (float)u_localptr[usr_k][usr_j][usr_i]; /* LHS is my 3D array, RHS is PetSc's */
    DMDAVecRestoreArray(da,u,&u_localptr);

    DMDAVecGetArray(da,u_old,&u_localptr);
	for(usr_k = 0; usr_k < usr_MZ; usr_k++)
	for(usr_j = 0; usr_j < usr_MY; usr_j++)
	for(usr_i = 0; usr_i < usr_MX; usr_i++)
	 usr_u0_loc_old[usr_k][usr_j][usr_i] = (float)u_localptr[usr_k][usr_j][usr_i]; /* LHS is my 3D array, RHS is PetSc's */
    DMDAVecRestoreArray(da,u,&u_localptr);


/**********************************************************************************************/
/**********************************************************************************************/
// do some measurements




// do the ecg and write it.
ecg = 0.0;
xl = 64; yl = -50; zl = 30;
	for(usr_k = 0; usr_k < usr_MZ; usr_k++)
	for(usr_j = 0; usr_j < usr_MY; usr_j++)
	for(usr_i = 0; usr_i < usr_MX; usr_i++){

if(user.geometry[usr_k][usr_j][usr_i][0]>0&&user.geometry[usr_k][usr_j][usr_i+1][0]>0&&user.geometry[usr_k][usr_j][usr_i-1][0]>0&&user.geometry[usr_k][usr_j+1][usr_i][0]>0&&user.geometry[usr_k][usr_j-1][usr_i][0]>0&&user.geometry[usr_k+1][usr_j][usr_i][0]>0&&user.geometry[usr_k-1][usr_j][usr_i][0]>0){

rx = (float)(usr_i - xl)*DX; ry = (float)(usr_j - yl)*DY; rz = (float)(usr_k - zl)*DZ;
r = sqrt(rx*rx + ry*ry + rz*rz);
r3 = r*r*r;

ecg = ecg - (usr_u0_loc[usr_k][usr_j][usr_i+1] - usr_u0_loc[usr_k][usr_j][usr_i-1])*DX/r3 \
          - (usr_u0_loc[usr_k][usr_j+1][usr_i] - usr_u0_loc[usr_k][usr_j-1][usr_i])*DY/r3 \
          - (usr_u0_loc[usr_k+1][usr_j][usr_i] - usr_u0_loc[usr_k-1][usr_j][usr_i])*DZ/r3;
}

}

// now write it.
output = fopen("lineDiagrams.dat","a+");
// time, ecg, san centre AP, san periphery AP, paranodal AP, atrial AP.
fprintf(output,"%d\t%f\t%f\t%f\t%f\t%f\n",file_Counter,ecg, usr_u0_loc[41][50][50], usr_u0_loc[22][56][37], usr_u0_loc[35][70][56], usr_u0_loc[35][70][70]);
fclose(output);

/**********************************************************************************************/
/**********************************************************************************************/

// this is in your time loop for now.
   VecDestroy(&u); 
   VecDestroy(&u_old); 
   DMDestroy(&da); 

 } // fileCounter loop end


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   free(usr_u0_loc);
   free(usr_u0_loc_old);
   free(dvdt);
  free(dvdt2);

   PetscFinalize();
   PetscFunctionReturn(0);
}

