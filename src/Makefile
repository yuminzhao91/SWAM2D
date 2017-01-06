FC=gfortran
FFLAGS=-W -Wall -O3 -ffast-math -fno-protect-parens -fstack-arrays -march=native
LDFLAGS=
BIN=../bin/
EXEC=fdpsv
#SRCF=$(wildcard *.f08)
SRC=param.f08 acquisition.f08 types.f08 boundaries.f08 deriv.f08 marching.f08 source.f08 model.f08 main.f08
OBJ=$(SRC:.f08=.o)

all: $(BIN)/$(EXEC)

$(BIN)/$(EXEC): $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.f08
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean mrproper

clean:
	@rm -rf *.o *.mod

mrproper: clean
	@rm -rf $(BIN)/$(EXEC)
