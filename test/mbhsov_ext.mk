PROJECT=MBHSOVEXT   # for historical reasons we always like to 
                     # call the executable int2, but you could set it to 
                     # something more descriptive

###HOST=windows
#HOST=linux
#HOST=macosx
HOST=linux-gfortran
#HOST=linux-gfortran-openmp


ifeq ($(HOST),linux)

OBJSUF=o
MODSUF=mod
FC=f77 -c 
FFLAGS=-fast
FLINK=f77 -o $(PROJECT)
WITH_SECOND=1

else

ifeq ($(HOST),macosx)

OBJSUF=o
MODSUF=mod
FC=gfortran -c 
FFLAGS=-O2
FLINK=gfortran -o $(PROJECT)

else

ifeq ($(HOST),linux-gfortran)

# buggy compiler, do not use -O2
OBJSUF=o
MODSUF=mod
FC=gfortran -c 
FFLAGS=-g
FLINK=gfortran -g -o $(PROJECT)

else

ifeq ($(HOST),linux-gfortran-openmp)

# buggy compiler, proceed anyway, use gfortran > 4.4.0
OBJSUF=o
MODSUF=mod
FC=gfortran -c 
FFLAGS=-O3 --openmp
FLINK=gfortran -o $(PROJECT) --openmp
### export OMP_NUM_THREADS=4
### export OMP_STACKSIZE=1024M

else

ifeq ($(HOST),linux-intel-openmp)

OBJSUF=o
MODSUF=mod
FC=ief77 -c
FFLAGS= -xT -O3 -ip -openmp 
FLINK=ief77 -o $(PROJECT) -xT -O3 -ip -openmp -static
WITH_SECOND=1

else 

ifeq ($(HOST),linux-intel)

OBJSUF=o
MODSUF=mod
FC=ief77 -c
FFLAGS=-xT -O3 -ip
FLINK=ief77 -o $(PROJECT) -xT -O3 -ip
WITH_SECOND=1

else 

ifeq ($(HOST),linux-intel-profile)

OBJSUF=o
MODSUF=mod
FC=ief77 -c
FFLAGS=-xT -O3 -ip -p
FLINK=ief77 -o $(PROJECT) -xT -O3 -ip -p 
WITH_SECOND=1

else 

ifeq ($(HOST),linux-intel64)

OBJSUF=o
MODSUF=mod
FC=ief77 -c
FFLAGS=-xT -O3 -ip -mcmodel medium -i-dynamic 
FLINK=ief77 -o $(PROJECT) -xT -O3 -ip -mcmodel medium -i-dynamic 
###export LD_LIBRARY_PATH=/opt/intel/fce/10.1.008/lib
WITH_SECOND=1

else 

ifeq ($(HOST),linux-lahey64)

OBJSUF=o
MODSUF=mod
FC=lf77_x64 -c
FFLAGS=-fast --model medium
FLINK=lf77_x64 -o $(PROJECT) -fast --model medium /opt/lf6481/lib64/libelf.so.0
###FLINK=lf77_x64 -o $(PROJECT) -fast --model medium 
WITH_SECOND=1

else 
     
ifeq ($(HOST),linux-fort77)

OBJSUF=o
MODSUF=mod
FC=fort77 -c
FFLAGS=-fast 
FLINK=fort77 -o $(PROJECT) -static -g
WITH_SECOND=1

else 

ifeq ($(HOST),linux-profile)

OBJSUF=o
MODSUF=mod
FC=fort77-i386-m32 -c
FFLAGS=-fast -profile -pg
FLINK=fort77-i386-m32 -o $(PROJECT) -static -pg -g
WITH_SECOND=1

else 

OBJSUF=obj
MODSUF=mod
FC=lf95 -c
FFLAGS=-O1
FLINK=lf95 -out $(PROJECT)

endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif


.PHONY: $(PROJECT) clean list

.f.$(OBJSUF):
	$(FC) $(FFLAGS) $<

.f90.$(OBJSUF):
	$(FC) $(FFLAGS) $<

.f.$(MODSUF):
	$(FC) $(FFLAGS) $<

.SUFFIXES: $(MODSUF) .$(OBJSUF) .f .c .f90

# SOURCE FILE LIST
#
vpath %.f .:../src

vpath %.f90 .

FMODS = 

F90SRCS = mbhsov_ext_dr.f90

FSRCS = modbhrouts.f helmrouts2d_quad.f helmrouts2d.f \
	laprouts2d.f laprouts2d_quad.f hank103.f prini.f \
	cdjseval2d.f dfft.f hkrand.f dlaran.f  svd.f \
	mbhrouts2dsov.f yukrouts2d.f mbhrouts2d.f

ifeq ($(WITH_SECOND),1) 
FSRCS += second-r8.f
endif

#
# object files list
MODS    =  $(FMODS:.f=.$(MODSUF)) 
OBJS    =  $(FMODS:.f=.$(OBJSUF)) $(FSRCS:.f=.$(OBJSUF)) \
	$(F90SRCS:.f90=.$(OBJSUF)) 
#
$(PROJECT):   $(MODS)   $(OBJS)
	rm -f $(PROJECT)
	$(FLINK) $(OBJS)
	./$(PROJECT)
#
clean: 
	rm -f $(OBJS)
# 
list: $(FSRCS)
	echo $^
#
pack: $(FSRCS)
	cat $^ > _tmp_.pack.f
#
