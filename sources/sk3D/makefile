
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

MYSUNDIALS_INCLUDE = -I/home/srk209/sundials/include
MYSUNDIALS_LIB     = -L/home/srk209/sundials/lib -lsundials_cvode -lsundials_nvecserial -lsundials_nvecparallel
MYBLAS_INCLUDE	   = -I/home/BLAS
MYBLAS_LIB	   = -L/home/srk209/BLAS -llibfblas
MYLAPACK_INCLUDE   = -I/home/srk209/lapack-3.5.0/lapacke/include
MYLAPACK_LIB	   = -L/home/srk209/lapack-3.5.0/ -lliblapack.a

CPPFLAGS	   = ${MYSUNDIALS_INCLUDE}

sk3d: 	sk3d.o  chkopts
	-${CLINKER} -o sk3d sk3d.o ${PETSC_TS_LIB} ${MYSUNDIALS_LIB} -lm -inline-calloc
#	mpirun -n 4 ./sk3d 5

# make simple vtks to see the time series 2D/3D data.
vt: 	vtk.o  chkopts
	-${CLINKER} -o vt vtk.o  ${PETSC_TS_LIB} -lm -inline-calloc

# do some activation time.
at: 	activationtime.o  chkopts
	-${CLINKER} -o at activationtime.o  ${PETSC_TS_LIB} -lm -inline-calloc

ecg:	ecg_ap.o  chkopts
	-${CLINKER} -o ecg ecg_ap.o  ${PETSC_TS_LIB} -lm -inline-calloc

run:
	mpiexec -n 3 ./sk3d &

# -draw_pause 0.01 -ts_monitor_draw_solution
	
edit:
	gedit sk3d.c makefile &

veryclean:
	rm -rf *.dat *.vtk *.vts ex15 *.o sk2d sk3d vt *~ mes *.bin*
