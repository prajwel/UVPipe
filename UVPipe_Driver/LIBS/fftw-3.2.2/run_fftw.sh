cd $1
./configure --prefix=$2
#./configure --prefix=$1
#echo $1
make;
make install;
#cp fitsio.h ../include
