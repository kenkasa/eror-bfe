#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS ${PBS_O_WORKDIR}
cd ${PBS_O_WORKDIR}
exe=../rr_analysis.x

files="zened/run*.zened"

######

ls $files | awk '{print $1}' > flist
mkdir -p dat
mkdir -p btstrp
  dat=dat/test.dat
  rp=btstrp/test
cat << EOF > run.inp
&input_param
  flist = "flist"
/

&output_param
  fhead = "$rp"
/

&option_param
  calcfe        = .true.
  ndim = 2 
  dr            = 0.2d0
  btstrp        = .true.
  ngrid         = 1000
  nsta          = 9875
  temperature   = 298.0d0
  vol0          = 1661.0d0
  urange        = 25.0d0 30d0
  react_range   = -500.0d0 500.0d0 -12.0 -4.0
/

&btstrp_param
  duplicate     = .true.
  seed_input    = .false.
  gen_num       = 1000
  samples       = 1000
/
EOF
  $exe run.inp >& $dat

#rm -f run.inp flist
