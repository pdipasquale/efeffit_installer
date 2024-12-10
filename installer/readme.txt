Just a list of things to do for me to follow what needs to be done in the installer executable:

install dependancies:
    edit sources.list

    in C, command sequence will be:
            FILE *sourcesFile = fopen("/etc/apt/sources.list", "a");
            fputs("##################\n", sourcesFile);
            fputs("deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe\n", sourcesFile);
            fputs("deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe\n", sourcesFile);
            fputs("deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe\n", sourcesFile);
            fputs("deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe\n", sourcesFile);
            fputs("##################\n", sourcesFile);


sudo apt update

sudo apt install g77