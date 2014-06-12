#cpl = ifort
#lflags_precise = -fp-model precise -fp-model source -fimf-arch-consistency=true
#lflags_debug = -debug -save-temps -fpic -heap-arrays -O2 -g -traceback -check all -fp-stack-check -vec-report5 $(lflags_precise)
#lflags_fast = -fast
#lflags_profiling = $(lflags_fast) -profile-functions -profile-loops=all -profile-loops-report=2

cpl	= gfortran
lflags_precise = -g -Wall -fcheck=all -ffpe-trap=zero,invalid,overflow,underflow -pedantic-errors
lflags_debug = $(lflags_precise)
lflags_fast = -O3
lflags_profiling = $(lflags_fast) -pg

OBJS = common_MC_MO.o common_MC_MO_TP.o subroutines_share_triv.o main.o \
       subroutines_share_proc.o subroutines_mo.o subroutines_proc_reac.o opkd*.o

LIBS = -lm -lcfitsio

exe_name = isotopize
debug_name = isotopize_debug
prof_name = isotopize_prof

all: lflags = $(lflags_fast)
all: $(exe_name)

precise: lflags = $(lflags_fast) $(lflags_precise)
precise: $(exe_name)

debug: lflags = $(lflags_debug)
debug: $(debug_name)

prof: lflags = $(lflags_profiling)
prof: $(prof_name)

cflags = $(lflags) -c

$(exe_name):  $(OBJS)
	$(cpl) $(lflags) -o $(exe_name) $(OBJS) $(LIBS)
	echo -e '\a'

$(debug_name):  $(OBJS)
	$(cpl) $(lflags) -o $(debug_name) $(OBJS) $(LIBS)
	echo -e '\a'

$(prof_name):  $(OBJS)
	$(cpl) $(lflags) -o $(prof_name) $(OBJS) $(LIBS)
	echo -e '\a'

$(exe_alt):  $(OBJS)
	$(cpl) $(lflags) -o $(exe_alt) $(OBJS) $(LIBS)
	echo -e '\a'

main.o: main.f90 configure.o disk.o sub_trivials.o montecarlo.o ray_tracing.o hitran.o lamda.o cdms.o ray_propagating.o
	$(cpl) $(cflags) main.f90

disk.o: disk.f90 chemistry.o heating_cooling.o load_Visser_CO_selfshielding.o data_struct.o grid.o montecarlo.o data_dump.o load_Draine_dusts.o ray_propagating.o vertical_structure.o
	$(cpl) $(cflags) disk.f90

vertical_structure.o: vertical_structure.f90 data_struct.o sub_trivials.o grid.o
	$(cpl) $(cflags) vertical_structure.f90

sync:
	#rsync -acvz *.f90 fdu@shire.astro.lsa.umich.edu:~/now/src/

syncall:
	#rsync -acvz *.f90 fdu@shire.astro.lsa.umich.edu:~/now/src/
	#rsync -acvz makefile fdu@shire.astro.lsa.umich.edu:~/now/src/
clean:
	#rm ./*.mod ./*.o
