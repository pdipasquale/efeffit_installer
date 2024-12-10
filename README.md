# I wrote an installer file for this
If you want to use the installer file run the following command in a terminal:

gcc installer/full_install.c -o install.o || sudo ./install.o

# ifeffit_lucas

# Need to install a VM with a particular Ubuntu version
# For this we'll install 16.04, once the efeffit install is done you can upgrade your ubuntu freely.

# First we need to install a g77 compiler, so we need to change our source lists up a bit:

# run this to open the list of sources:
sudo gedit /etc/apt/sources.list

# Scroll to the bottom of the text file that just opened and paste these to the very bottom of the document:

deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe

deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe

deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe

deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe

# Now in a terminal type:
sudo apt update

sudo apt install g77

#That should update the list of sources, and then the second command will let us install the g77 compiler

#if you want to test that the compiler works, just type 'g77' into a terminal and you should see a compilation error (rather than a "command not found" error)

# open a terminal in the directory of this readme and copy the ifeffit lucas folder into your lib directory:
sudo cp ifeffit-1.2.5_lucas /usr/lib

# and open another terminal in /usr/lib/ifeffit-1.2.5_lucas
cd /usr/lib/ifeffit-1.2.5_lucas

#now in this terminal we can start installing ifeffit
#First by installing pgplot, run the following command, which will fail

sudo ./PGPLOT_install
# we will first see an error about png.h, xos.h, etc.
#need to install the renderer for png files
sudo apt-get install libpng-dev

sudo apt-get install libx11-dev

# Now to resolve an error about a missing png.h rule, run the following:
sudo mkdir /usr/local/pgplot
sudo ln -s /usr/include/png.h /usr/local/pgplot/png.h

sudo ln -s /usr/include/pngconf.h /usr/local/pgplot/pngconf.h

sudo ln -s /usr/include/zlib.h /usr/local/pgplot/zlib.h

sudo ln -s /usr/include/zconf.h /usr/local/pgplot/zconf.h

#note that for you zconf might be in: sudo ln /usr/include/x86_64-linux-gnu/zconf.h             Check to make sure
# should throw an error that crt1.o, crti.o and crtn.o are missing, also that libgcc_s.so is missing
#run the following commands and it should resolve those errors:

sudo ln -s /usr/lib/x86_64-linux-gnu/crti.o /usr/lib/

sudo ln -s /usr/lib/x86_64-linux-gnu/crt1.o /usr/lib/ 

sudo ln -s /usr/lib/x86_64-linux-gnu/crtn.o /usr/lib/

sudo ln -s /usr/lib/gcc/x86_64-linux-gnu/4.8/libgcc_s.so /usr/lib/x86_64-linux-gnu/

# Now after all that ./PGPLOT_install should work fine, and we can install ifeffit!
# in the terminal in the ifeffit-lucas directory, ("cd /usr/lib/ifeffit-1.2.5_lucas" if you can't find it)
sudo ./configure

sudo make

# should throw an error about some missing files, we need to install:
sudo apt-get install libtinfo-dev

sudo ln -s /usr/lib/x86_64-linux-gnu/libtermcap.a /lib/termcap/

# now try make again
sudo make

sudo make install

# Now you're (probabaly) done!

To test that the intallation worked, open up a terminal anywhere and just type 'ifeffit' and hit enter. it should load into a weird environment with a '>>'
We can now use this to feed commands into ifeffit.

ifeffit itself behaves like a scripting language, you'll find that you can write scripts that ifeffit can read in and run, for example:
    if you make a fit_script.iff
    to run it you would just write into your terminal:
    "ifeffit fit_script.iff"
    similar to python

# ENJOY :)
# Note that i plan on updating this repo with an installation executable, that resolves all of these weird errors, but it will only work on x86 processors.
