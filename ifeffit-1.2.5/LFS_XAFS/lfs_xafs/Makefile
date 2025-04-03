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
LIB_IFF  = -L/cygdrive/c/AnalysisTauhid/IFEFIT/ifeffit-1.2.5_lucas/src/lib/ -lifeffit 
LIB_PLT  = -L/usr/local/pgplot -lpgplot -L/lib -lpng
LIB_F77  = 
LIB_X11  = -L/lib -lX11

INC_IFF  = 
#==  


F77       = g77
LOAD      = $(F77)

TARGET = lfs_xafs
SRC    = lfs_xafs.f
OBJ    = lfs_xafs.o

default: $(TARGET)


.f.o:
	$(F77) -c $*.f

ALL_LIBS = $(LIB_IFF) $(LIB_PLT) $(LIB_F77) $(LIB_X11)

$(OBJ): /cygdrive/c/AnalysisTauhid/IFEFIT/ifeffit-1.2.5_lucas/src/lib/libifeffit.a

$(TARGET): $(OBJ)
	$(LOAD) $(OBJ) -o $(TARGET) $(ALL_LIBS)

clean :
	rm -f $(OBJ) $(TARGET)

install: default
	cp lfs_xafs ~/bin/
