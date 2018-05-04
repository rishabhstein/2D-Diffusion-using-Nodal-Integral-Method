#the compiler: g95 for FORTRAN program, define as g++ for C++ and gcc for C
g95 = gfortran

  # compiler flags:
  #  -g    adds debugging information to the executable file
  #  -ffree-line-length-none is to compile FORTRAN-90 code
FLAGS  = -g -ffree-line-length-none

  # the build target executable:
TARGET = Diffusion
EXEC=	Diffu
FILES=	T_P1 T_P2 Ts Tt Xs Ys Xt Yt
LOGFILE= log.txt
all: $(TARGET)

$(TARGET): $(TARGET).f90
	$(g95) $(FLAGS) -o $(EXEC) $(TARGET).f90

clean:
	$(RM) $(EXEC) $(FILES) $(LOGFILE)

run:
	./$(EXEC)
