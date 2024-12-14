#!/usr/bin/env bash

export exe=../../../packages/eror/bin/ermod
export OMP_NUM_THREADS=1

mpiexec.hydra -bootstrap rsh -np 128 $exe >& ermod.log
cp engsln.tt engsln.01
cp engint.tt engint.01
