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
LIB_IFF  = -L/home/lucasfs/ifeffit/lib -lifeffit 
LIB_PLT  = -L/home/lucasfs/pglib -lpgplot 
LIB_F77  = -lg2c -lm -L/usr/lib/gcc-lib/i386-linux/3.3.1 -L/usr/lib/gcc-lib/i386-linux/3.3.1/../../.. -lm -lfrtbegin -lgcc_s 
LIB_X11  = -L/usr/X11R6/lib -lX11
INC_IFF  = -I/home/lucasfs/ifeffit/share/ifeffit/congif 
#==  

F77       = f77 -O1
LOAD      = $(F77)

TARGET = mu2chi
SRC    = mu2chi.f
OBJ    = mu2chi.o

default: $(TARGET)


.f.o:
	$(F77) -c $*.f

ALL_LIBS = $(LIB_IFF) $(LIB_PLT) $(LIB_X11) $(LIB_F77) 

$(TARGET): $(OBJ)
	$(LOAD) $(OBJ) -o $(TARGET) $(ALL_LIBS)

clean :
	rm -f $(OBJ) $(TARGET)

