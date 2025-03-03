#!/usr/bin/env bash

exe=../../packages/eror/tools/GENESIS/gen_structure_AMBER
prmtop=../01_amber_setup/host-guest_in_water/hg_in_water.prmtop

# 1. Prepare soln and refs directories 
#   soln : working directory for rho in the solution  system
#   refs : working directory for rho and chi in the reference system
python3 $exe  \
  -t $prmtop  \
  -s APR

# 2. Copy files in refs directories to dGcorr
#
rsync -av --exclude "run.sh" --exclude "gen_input.sh" refs/ dGcorr/

