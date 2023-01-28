
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

MYSUNDIALS_INCLUDE = -I/home/srk209/sundials/include
MYSUNDIALS_LIB     = -L/home/srk209/sundials/lib -lsundials_cvode -lsundials_nvecserial -lsundials_nvecparallel
MYBLAS_INCLUDE	   = -I/home//BLAS
MYBLAS_LIB	   = -L/home/srk209/BLAS -llibfblas
MYLAPACK_INCLUDE   = -I/home/srk209/lapack-3.5.0/lapacke/include
MYLAPACK_LIB	   = -L/home/srk209/lapack-3.5.0/ -lliblapack.a

CPPFLAGS	   = ${MYSUNDIALS_INCLUDE}


geom:
	gcc -o geom myGeometry.c -lm

sk2d: 	sk2dhSAN.o  chkopts
	-${CLINKER} -o sk2d sk2dhSAN.o  ${PETSC_TS_LIB} ${MYSUNDIALS_LIB} -lm -inline-calloc
	${RM} sk2dhSAN.o
#	time mpiexec -n 4 ./sk2d 17

mes: 	measurement.o  chkopts
	-${CLINKER} -o mes measurement.o  ${PETSC_TS_LIB} -lm -inline-calloc
#	${RM} sk2d.o
#	./mes

edit:
	gedit sk2d.c f.c makefile measurement.c &

veryclean:
	rm -rf *.dat *.vtk *.vts ex15 *.o sk2d *~ *.bin *.info mes
