#!/usr/bin/env bash

export exe=/data2/kasa/packages/GitHub/vmd_scripts/ermod-0.3.7/amd/bin/ermod
export OMP_NUM_THREADS=1

mpiexec.hydra -bootstrap rsh -np 128 $exe >& ermod.log

