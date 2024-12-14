#!/usr/bin/env bash

cwd=`pwd`

echo "Install ANATRA"
#
cd packages/anatra
./run_install.sh
cd $cwd

echo "Install EROR"

cd packages/eror
./run_install.sh
cd $cwd
