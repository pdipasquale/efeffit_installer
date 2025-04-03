##
## Fortran Makefile
##
## 1: replace with the contents of Config.mak from
##    /usr/local/share/ifeffit/config

#==  Ifeffit build configuration:
#    LIB_IFF = ifeffit library  
#    LIB_PLT = PGPLOT  libraries
#    LIB_F77 = Fortran libraries
#    LIB_X11 = X Libaries
#    INC_IFF = location of ifeffit.h
LIB_IFF  = -L/usr/local/lib -lifeffit
LIB_PLT  = -L/usr/local/pgplot -lpgplot -lpng -lz -lX11
LIB_F77  = -L/usr/lib/gcc/x86_64-linux-gnu/3.4.6 -L/usr/lib/gcc/x86_64-linux-gnu/3.4.6/../../../../lib -L/lib/../lib -L/usr/lib/../lib -lfrtbegin -lg2c -lm -lgcc_s 
INC_IFF  = -I/usr/local/share/ifeffit/config
#==  

#F77       = f77 -O1
#LOAD      = $(F77)

LIB_IFF  = -L/cygdrive/c/AnalysisTauhid/IFEFIT/ifeffit-1.2.5_lucas/src/lib/ -lifeffit 
LIB_PLT  = -L/usr/local/pgplot -lpgplot -L/lib -lpng
LIB_F77  = 
LIB_X11  = -L/lib -lX11

INC_IFF  = 

F77       = g77
LOAD      = $(F77)


TARGET = mu2chi
SRC    = mu2chi.f
OBJ    = mu2chi.o

default: $(TARGET)


.f.o:
	$(F77) -c $*.f

#ALL_LIBS = $(LIB_IFF) $(LIB_PLT) $(LIB_X11) $(LIB_F77) 

ALL_LIBS = $(LIB_IFF) $(LIB_PLT) $(LIB_F77) $(LIB_X11)




$(TARGET): $(OBJ)
	$(LOAD) $(OBJ) -o $(TARGET) $(ALL_LIBS)

clean :
	rm -f $(OBJ) $(TARGET)

