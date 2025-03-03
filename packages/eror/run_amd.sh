#!/bin/bash

make clean
./configure --with-mkl --prefix=`pwd`/amd
make
make install
