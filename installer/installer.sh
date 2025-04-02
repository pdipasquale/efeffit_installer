#!/bin/bash

FILE="./ifeffit_sources.list"

while IFS= read -r line; do
    if grep -Fxq "$line" /etc/apt/sources.list; then
        echo "$line already in sources.list"
    else
        echo "$line not found in sources.list, adding..."
        echo "$line" | sudo tee -a /etc/apt/sources.list > /dev/null
    fi
done < "$FILE"

sudo apt update
sudo apt install -y -q g77

echo "Copying ifeffit to /usr/lib/"
sudo cp -r ../ifeffit-1.2.5 /usr/lib/

echo "Installing libpng-dev & libx11-dev"
sudo apt-get install -y -q libpng-dev
sudo apt-get install -y -q libx11-dev

sudo mkdir /usr/local/pgplot

echo "Linking files to correct directories"
sudo ln -s /usr/include/png.h /usr/local/pgplot/png.h
sudo ln -s /usr/include/pngconf.h /usr/local/pgplot/pngconf.h
sudo ln -s /usr/include/zlib.h /usr/local/pgplot/zlib.h
sudo ln -s /usr/include/zconf.h /usr/local/pgplot/zconf.h

sudo ln -s /usr/lib/x86_64-linux-gnu/crti.o /usr/lib/
sudo ln -s /usr/lib/x86_64-linux-gnu/crt1.o /usr/lib/
sudo ln -s /usr/lib/x86_64-linux-gnu/crtn.o /usr/lib/
sudo ln -s /usr/lib/gcc/x86_64-linux-gnu/4.8/libgcc_s.so /usr/lib/x86_64-linux-gnu/

echo "Installing PGPLOT"
sudo /usr/lib/ifeffit-1.2.5/PGPLOT_install

sudo apt-get install -y -q libtinfo-dev
sudo mkdir /lib/termcap/
sudo ln -s /usr/lib/x86_64-linux-gnu/libtermcap.a /lib/termcap/

echo "Dependencies have been handled. Now installing ifeffit."
cd /usr/lib/ifeffit-1.2.5
sudo ./configure
sudo make
sudo make install

echo "ifeffit is now installed!\n Test out by typing 'ifeffit' into a command window and type feffit2 to check that the custom version is working. You should see this error: 'feffit2: no chi(k) data array?'"
