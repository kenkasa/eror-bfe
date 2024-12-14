#!/usr/bin/env bash

# Usage: ./install.sh <compiler type> ("gcc" or "intel")
#
compiler=$1   # "intel" or "gcc"
#
#=====================================================================
#
# ... install xdr library
#
cwd=`pwd`
XDRPATH=$ANATRA_PATH/f90/lib/external/xdr-interface-fortran

if [ "$compiler" != "fugaku" ]; then
  cd $XDRPATH
    if [ ! -e xdrfile-1.1.4 ]; then
      tar xvf xdrfile-1.1.4.tar.gz 
    fi
    cd xdrfile-1.1.4
    ./configure --prefix=$XDRPATH/xdrfile-1.1.4
    make && make install
  cd $cwd
fi
#exit
#
# ... install ANATRA fortran programs
#
if [ "$compiler" == "" ]; then
  compiler=gcc
  echo "no compiler type is specified"
  echo ">> gcc is used for compile"
elif [ "$compiler" == "gcc" ]||[ "$compiler" == "intel" ]; then
  echo "$compiler is used"
fi

list="comcrd_analysis      \
      lipidorder_analysis  \
      zp_analysis     \
      hs_analysis     \
      rd_analysis     \
      dp_analysis     \
      or_analysis     \
      rd_gather       \
      rp_analysis     \
      ed_analysis     \
      pd_analysis     \
      rr_analysis     \
      ed_analysis     \
      pm_analysis     \
      sd_analysis     \
      sd_gather       \
      energy_analysis \
      ss_analysis     \
      rot_analysis    \
      sr_analysis     \
      movave_analysis \
      gderp_analysis  \
      rprate_analysis \
      funcave_analysis \
      slicedx_analysis \
      elecond_analysis"

cwd=`pwd`
mkdir -p bin
for d in $list;do
  echo "o Installing $d ..."
  echo ""
  if [ "$compiler" == "gcc" ]&&[ "$d" == "en_analysis" ]; then
    echo "Compiler: gcc  Program: en_analysis"
    echo ">> Skipped"
    echo ""
    continue
  fi
  cd $d

  make -f Makefile.$compiler

  #if [ -e ${d}.x ]; then
  #  cp ${d}.x ../bin/
  #fi

  cd $cwd 
  echo ">> Finished"
  echo "" 
done

chk=0
for d in $list;do
  if [ ! -e ./bin/${d}.x ];then

    if [ "$compiler" == "gcc" ]&&[ "$d" == "en_analysis" ]; then
      continue
    fi

    echo "Installation of $d is failed."
    echo "Please contact the developers"
    echo "if the problem is due to bugs."
    echo ""
    chk=`expr $chk + 1` 
  fi 
done

if [ $chk -eq 0 ];then
  echo "-------------------------------------------------"
  echo "Installation of ANATRA fortran programs have been"
  echo "succesfully finished!!"
  echo "-------------------------------------------------"
else
  echo "-------------------------------------------------"
  echo "$chk errors occured during the installation."
  echo "Installation terminated abnormally."
  echo "-------------------------------------------------"
fi

