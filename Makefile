# Makefile compilation parameters
CC = gfortran
CFLAGS = -O1 -fbounds-check -fbacktrace -fcheck=all -g -Wall -Wextra -Wrealloc-lhs-all
TARGET = ExecuteOptimization
EXT = -llapack -lblas -L/opt/homebrew/Cellar/metis/5.1.0/lib -lmetis -fopenmp
FILES = Modules/MA87Routines/sdeps90.f90 \
		Modules/MA87Routines/ddeps90.f90 \
		Modules/MA87Routines/common90.f90 \
		Modules/MA87Routines/hsl_ma87s.f90 \
		Modules/MA87Routines/hsl_ma87d.f90 \
		Modules/MMARoutines/MMA_Routines.f90 \
		Modules/MMARoutines/MMA_Interface.f90 \
		Modules/Base_Module.f90 \
		Modules/Solver_MA87Module.f90 \
		Modules/FEA_Module.f90 \
		Modules/Optimization_Module.f90 \
		Modules/Paraview_Module.f90 \
		MainOptimization.f90
DIR = Paraview DataResults DataResults/.InternalData
# Object files
OBJ1 = ${FILES:.f90=.o}

all: create_dirs $(TARGET)

# dir
create_dirs:
	mkdir -p $(DIR)

# compilation and cleanup command
%.o : %.f90
	${CC} ${CFLAGS} -o $@ -c $<

${TARGET} : ${OBJ1}
	${CC} ${CFLAGS} -o $@ ${OBJ1} $(EXT)

.PHONY : clean
clean :
	@rm -f *.o *.mod ${TARGET} ${OBJ1} DataResults/*.txt DataResults/.InternalData/*.txt \
	Paraview/*.geom Paraview/*.esca *.case