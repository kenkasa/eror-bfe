#!/bin/bash

prefix=`pwd`
./configure --with-mkl --prefix=$prefix

make 
make install
