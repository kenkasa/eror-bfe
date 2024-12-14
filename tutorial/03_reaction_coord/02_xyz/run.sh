#!/bin/bash

anatra trjconv                                 \
  -stype      pdb                              \
  -sfile      refbcd.pdb                       \
  -tintype    pdb                              \
  -tin        refbcd.pdb                       \
  -totype     xyz                              \
  -to         refbcd.xyz                       \
  -sel0       name O15 O16 O17 O18 O19 O20 O21 \
  -outselid   0
