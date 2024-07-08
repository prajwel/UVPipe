cd $1
echo -------------
echo $1
echo $2
echo --------------
./configure --prefix=$2
make;
make install;
#cp fitsio.h ../include
