
PROG=COSMO-d2-2hybrid

all: $(PROG)


$(PROG): $(PROG).f90
	gfortran $(PROG).f90 -I /usr/include -o $(PROG) -leccodes  -leccodes_f90 #&& ./$(PROG)

