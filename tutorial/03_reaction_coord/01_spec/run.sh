#!/bin/bash

prmtop=../../01_amber_setup/host-guest_in_water_for_insertion/hg_in_water_for_ins.prmtop
inpcrd=../../01_amber_setup/host-guest_in_water_for_insertion/hg_in_water_for_ins.inpcrd

anatra specfile                                              \
  -stype   parm7                                             \
  -sfile   $prmtop                                           \
  -tintype rst7                                              \
  -tin     $inpcrd                                           \
  -fhead   out                                               \
  -sel0    resname APR                                       \
  -sel1    resname BCD                                       \
  -sel2    resname APR and noh                               \
  -sel3    resname BCD and name O15 O16 O17 O18 O19 O20 O21 


# Uattr
cp out_0.spec guest_uij.spec
cp out_1.spec host_uij.spec

# Fitting 
cp out_2.spec guest_fit.spec
cp out_3.spec host_fit.spec

# Ins
cp out_2.spec guest_ins.spec

rm -f out_*.spec
