#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main()
{
    // giving system command and storing return value.....
    int returnCode = system("echo I am running");
    FILE *deps_file = fopen("installer/deps_check.txt", "r");
    char c[100];  // Buffer to hold the line
    fscanf(deps_file, "%[^\n]", c);  // Read a line up to the newline character
    printf("reading deps_file, value: %s\n", c);

    FILE *sourcesFile = fopen("/etc/apt/sources.list", "a+");
    if (strcmp(c, "0") == 0) {  // Compare the string to "0"
    // Need a means of intelligently avoiding directly overwriting the sources file 
    //in case this breaks and needs to be run from the start
    // I'm using the work intelligently very loosely here
    //FILE *sourcesFile = fopen(".\\test_files\\sources.list", "a");
    // Basically if the installation fucks up we need to go into the sources file and delete this section
        fclose(deps_file);
        fputs("#we are editing the file here#\n", sourcesFile);
        fputs("##################\n", sourcesFile);
        fputs("deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe\n", sourcesFile);
        fputs("deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe\n", sourcesFile);
        fputs("deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe\n", sourcesFile);
        fputs("deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe\n", sourcesFile);
        fputs("##################\n", sourcesFile);
        FILE *deps_file_new = fopen("installer/deps_check.txt", "w+");
        
        fputs("1", deps_file);
        fclose(deps_file);   
    fclose(sourcesFile);
    }
    
    
    // installing the g77 compiler
    system("sudo apt update");
    const char *g77install = "sudo apt install -y g77";
    printf("installing g77 compiler");
    int result = system(g77install);
    
        
    if (result == -1) {
        printf("Failed to install dependencies. Error code: %d\n", result);
        printf("Try a manual g77 installation:");
        printf("sudo apt install g77");
        fputs("0", deps_file);
        fclose(deps_file);
        exit(0);

    } else {
	//no idea what this bullshit does, thank you to some guy on stackexchange
	int exit_status = WEXITSTATUS(result);
	printf("Raw exit status: %d\n", exit_status);

	
	if (exit_status == 0) {
            printf("g77 installed successfully.\n");
        } else if (exit_status == 1) {
            printf("Command completed successfully, but minor issues occurred (e.g., package already installed) Continuing with ifeffit installation.\n");
        } else {
            printf("Command failed with exit status: %d\n", exit_status);
		return 1;
        }
    }
    

    // move the ifeffit folder into the /usr/lib folder
    const char *copyifeffitcommand = "sudo cp ifeffit-1.2.5_lucas /usr/lib/ -r";
    system(copyifeffitcommand);

    if (result == 0) {
        printf("ifeffit directory copied to /usr/lib/.\n");
    } else {
        printf("Failed to copy ifeffit folder: %d\n", result);
    }
    // all of the installation files are located in /usr/lib/ifeffit-1.2.5_lucas

    system("sudo apt-get install -y libpng-dev");
    system("sudo apt-get install -y libx11-dev");
    printf("installed libpng-dev & libx11-dev \n");
    printf("linking files to correct directories \n");

    system("sudo ln -s /usr/include/png.h /usr/local/pgplot/");
    system("sudo ln -s /usr/include/pngconf.h /usr/local/pgplot/");
    system("sudo ln -s /usr/include/zlib.h /usr/local/pgplot/");
    system("sudo ln -s /usr/include/zconf.h /usr/local/pgplot/");

    system("sudo ln -s /usr/lib/x86_64-linux-gnu/crti.o /usr/lib/");
    system("sudo ln -s /usr/lib/x86_64-linux-gnu/crt1.o /usr/lib/ ");
    system("sudo ln -s /usr/lib/x86_64-linux-gnu/crtn.o /usr/lib/");
    system("sudo ln -s /usr/lib/gcc/x86_64-linux-gnu/5/libgcc_s.so /usr/lib/x86_64-linux-gnu/");

    system("sudo /usr/lib/ifeffit-1.2.5_lucas/PGPLOT_install");
    printf("PGPLOT is now installed!");
    
    system("sudo apt-get install -y libtinfo-dev");

    system("sudo mkdir /lib/termcap/");
    system("sudo ln -s /usr/lib/x86_64-linux-gnu/libtermcap.a /lib/termcap/");
    // the dependancies are now sorted, with that we can install ifeffit
    chdir("/usr/lib/ifeffit-1.2.5_lucas");
    system("sudo ./configure");
    system("sudo make");
    system("sudo make install");

    printf("ifeffit-lucas is now installed!\n Test out by typing 'ifeffit' into a command window and type feffit2 to check that the custom version is working.\n You should see this error: 'feffit2: no chi(k) data array?'");

    return 0;
}
