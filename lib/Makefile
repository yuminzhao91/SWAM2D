FC=gfortran
FFLAGS= -fpic -W -Wall -O2 -g
LDFLAGS= -shared
TARGET_LIB= libfdpsv.so
SRC=../src/param.f08 ../src/acquisition.f08 ../src/types.f08 ../src/boundaries.f08 ../src/deriv.f08 ../src/marching.f08 ../src/source.f08 ../src/model.f08
OBJ=$(SRC:.f08=.o)

all: $(TARGET_LIB)

$(TARGET_LIB): $(OBJ)
	$(FC) $(LDFLAGS) -o $@ $^

%.o: %.f08
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean mrproper

clean:
	@rm -rf *.o *.mod

mrproper: clean
	@rm -rf $(TARGET_LIB)
