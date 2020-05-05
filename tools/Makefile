
PROG=COSMO-d2-2hybrid

all: $(PROG)


$(PROG): $(PROG).f90
	gfortran $(PROG).f90  -o $(PROG) -I$(HOME)/.local/include  -L$(HOME)/.local/lib -leccodes  -leccodes_f90 #&& ./$(PROG)

