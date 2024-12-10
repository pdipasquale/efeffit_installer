#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
    // giving system command and storing return value
    int returnCode = system("echo I am running");
    
    FILE *sourcesFile = fopen("/etc/apt/sources.list", "a");
    // for testing use the file below
    // Need a means of intelligently avoiding directly overwriting the sources file 
    //in case this breaks and needs to be run from the start

    // I'm using the work intelligently very loosely here
    //FILE *sourcesFile = fopen(".\\test_files\\sources.list", "a");
        fputs("#we are editing the file here#\n", sourcesFile);
        fputs("##################\n", sourcesFile);
        fputs("deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe\n", sourcesFile);
        fputs("deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe\n", sourcesFile);
        fputs("deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe\n", sourcesFile);
        fputs("deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe\n", sourcesFile);
        fputs("##################\n", sourcesFile);

    fclose(sourcesFile);

    // installing the g77 compiler
    const char *g77install = "sudo apt install -y g77";
    printf("installing g77 compiler");
    int result = system(g77install);

    if (result == 0) {
        printf("Dependencies installed successfully.\n");
    } else {
        printf("Failed to install dependencies. Error code: %d\n", result);
    }

    // move the ifeffit folder into the /usr/lib folder
    const char *copyifeffitcommand = "sudo cp ifeffit-1.2.5_lucas /usr/lib";
    system(copyifeffitcommand);

    if (result == 0) {
        printf("ifeffit directory copied to /usr/lib/.\n");
    } else {
        printf("Failed to install dependencies. Error code: %d\n", result);
    }
    system("cd /usr/lib/ifeffit-1.2.5_lucas");

    system("sudo apt-get install -y libpng-dev");
    system("sudo apt-get install -y libx11-dev");
    printf("installed libpng-dev & libx11-dev \n");
    printf("linking files to correct directories \n");

    system("sudo ln -s /usr/include/png.h /usr/local/pgplot");
    system("sudo ln -s /usr/include/pngconf.h /usr/local/pgplot");
    system("sudo ln -s /usr/include/zlib.h /usr/local/pgplot");
    system("sudo ln -s /usr/include/zconf.h /usr/local/pgplot");

    system("sudo ln -s /usr/lib/x86_64-linux-gnu/crti.o /usr/lib/");
    system("sudo ln -s /usr/lib/x86_64-linux-gnu/crt1.o /usr/lib/ ");
    system("sudo ln -s /usr/lib/x86_64-linux-gnu/crtn.o /usr/lib/");
    system("sudo ln -s /usr/lib/gcc/x86_64-linux-gnu/5/libgcc_s.so /usr/lib/x86_64-linux-gnu/");

    system("sudo ./PGPLOT_install");
    printf("PGPLOT is now installed!");
    
    system("sudo apt-get install -y libtinfo-dev");

    system("sudo mkdir /lib/termcap/");
    system("sudo ln -s /usr/lib/x86_64-linux-gnu/libtermcap.a /lib/termcap/");
    // the dependancies are now sorted, with that we can install ifeffit
    system("sudo ./configure");
    system("sudo make");
    system("sudo make install");

    printf("ifeffit-lucas is now installed!\n Test out by typing 'ifeffit' into a command window and type feffit2 to check that the custom version is working.\n You should see this error: 'feffit2: no chi(k) data array?'");

    return 0;
}