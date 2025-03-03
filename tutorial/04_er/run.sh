#!/usr/bin/env bash

export exe=../../packages/eror/intel/bin/slvfe
export OMP_NUM_THREADS=28

cwd=`pwd`
cd refs/Long_Sample
unzip corref.01.zip
cd $cwd

cat << EOF > parameters_fe
&fevars
inptemp=298.000000
slndnspf  = 'engovl'
aveuvfile = 'avovl.tt'
solndirec = 'soln/Long_Sample'
refsdirec = 'refs/Long_Sample'
numdiv    = 5
/

EOF
$exe >& slvfe.out
