# >> FORTRAN compiler
FC=gfortran

# >> FORTRAN compiler flags for release
FFLAGS=-W -Wall -O3 -ffast-math -fno-protect-parens -fstack-arrays -march=native

# >> FORTRAN compiler flags for debug
FFLAGS_DEBUG=-W -Wall -C -g

# >> FORTRAN library flags
LDFLAGS=

# >> Binary folder
BIN=../bin/

# >> Executable name
EXEC=swam2d
EXECD=swam2d

#SRC=mod_types.f08 $(wildcard mod_*.f08) main.f08
SRC=mod_param.f08 mod_acquisition.f08 mod_types.f08 mod_boundaries.f08 mod_deriv.f08 mod_marching.f08 mod_source.f08 mod_model.f08 main.f08
OBJ=$(SRC:.f08=.o)

all: $(BIN)/$(EXEC)

debug: FFLAGS=$(FFLAGS_DEBUG)
debug: $(BIN)/$(EXEC)

$(BIN)/$(EXEC): $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.f08
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean mrproper

clean:
	@rm -rf *.o *.mod

mrproper: clean
	@rm -rf $(BIN)/$(EXEC)
