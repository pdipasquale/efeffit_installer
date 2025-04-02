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
