# ifeffit_lucas

# Need to install a VM with a particular Ubuntu version
# For this we'll install 14.04, once the efeffit install is done you can upgrade your ubuntu freely.

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
# That should update the list of sources, and then the second command will let us install the g77 compiler
# if you want to test that the compiler works, just type 'g77' into a terminal and ou should see a compilation error (rather than a "command not found" error)

