cd $1
./configure --prefix=$2
make;
make install;
#cp fitsio.h ../include
