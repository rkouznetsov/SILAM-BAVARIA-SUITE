
PROG=ICON-d2-2hybrid

all: $(PROG)


$(PROG): $(PROG).f90
	gfortran $(PROG).f90  -o $(PROG) -I$(HOME)/.local/eccodes/include  -L$(HOME)/.local/eccodes/lib -leccodes  -leccodes_f90 #&& ./$(PROG)

