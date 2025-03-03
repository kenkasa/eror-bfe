#!/bin/bash

make clean
./configure --with-mkl --prefix=`pwd`/intel
make
make install
