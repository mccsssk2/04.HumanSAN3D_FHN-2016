# 2 cells
# SHELL = /bin/sh

CC       = mpicc
CFLAGS   = -Wall -g -fno-inline
INCLUDE  = /home/srk209/sundials/include
MY_APP	 = hsan
LIB	 = -L/home/srk209/sundials/lib

hsan:	fkmodsk_2cell.c
	make veryclean
	${CC} ${CFLAGS} -I${INCLUDE} -c fkmodsk_2cell.c -fpic -o fkmodsk_2cell.o
	${CC} ${CFLAGS}  fkmodsk_2cell.o -I${INCLUDE} -g -lm ${LIB} -lsundials_cvodes -lsundials_nvecserial -mcmodel=large -o ${MY_APP}
#	make san

# the 600.0 is basal SAN. Make it 800.0 for paranodal AP. You can use the 5 as well. reduce 5 and increase CL.
san:
	./${MY_APP}
clean:
	rm -f *.o ex*.txt ${MY_APP} ${MY_APP1} *.dat

veryclean:
	rm -f *.o *.dat*
	rm -f ${MY_APP} ${MY_APP1} *~

edit:
	gedit *.c Makefile &

