proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                      Angle Analysis"
  puts ""
  puts "============================================================"
}

proc show_tr_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra angle                                                      \\"
    puts "  -stype      <structure file type>                               \\"
    puts "  -sfile      <structure file name>                               \\"
    puts "  -tintype    <input trajectory file type>                        \\"
    puts "  -tin        <input trajectory file name>                        \\"
    puts "  -fhead      <header of output file name>                        \\"
    puts "  -selX       <VMD selection> (X=0,1,2...)                        \\"
    puts "  -npoints    <3: bond-angle, 4: dihedral-angle>                  \\"
    puts "  -prep_only  <where analysis is performed or not (true or false) \\"
    puts "              (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra angle                  \\"
    puts "  -stype     parm7            \\"
    puts "  -sfile     str.prmtop       \\"
    puts "  -tintype   dcd              \\"
    puts "  -tin       inp.dcd          \\"
    puts "  -fhead     out              \\"
    puts "  -sel0      resid 1 and noh  \\"
    puts "  -sel1      resid 2 and noh  \\"
    puts "  -sel2      resid 3 and noh  \\"
    puts "  -npoints   3                \\"
    puts "  -prep_only false"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_cmoptinfo
#! @brief         define CoM option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_cmoptinfo {} {
  
  global cmopt

  set cmopt(fhead)      "out"
  set cmopt(npoints)     3
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_cmoptinfo
#! @brief         read CoM option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_cmoptinfo {arglist} {

  global cmopt

  set cmopt(fhead)      [parse_arguments $arglist \
      "-fhead"      "value" $cmopt(fhead)]
  set cmopt(mode)       [parse_arguments $arglist \
      "-mode"       "value" $cmopt(mode)]
  set cmopt(outmsd)     [parse_arguments $arglist \
      "-outmsd"     "value" $cmopt(outmsd)]
  set cmopt(unwrap)     [parse_arguments $arglist \
      "-unwrap"     "value" $cmopt(unwrap)]
  set cmopt(onlyz)      [parse_arguments $arglist \
      "-onlyz"      "value" $cmopt(onlyz)]
  set cmopt(msddim)     [parse_arguments $arglist \
      "-msddim"     "value" $cmopt(msddim)]
  set cmopt(dt)         [parse_arguments $arglist \
      "-dt"         "value" $cmopt(dt)]
  set cmopt(msdrange)   [parse_arguments $arglist \
      "-msdrange"   "value" $cmopt(msdrange)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_cmoptinfo
#! @brief         show CoM option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_cmoptinfo {} {

  global cmopt

  puts "<< option info >>"
  puts "fhead      = $cmopt(fhead)"
  puts "outmsd     = $cmopt(outmsd)"
  puts "onlyz      = $cmopt(onlyz)"
  puts "msddim     = $cmopt(msddim)"
  puts "dt         = $cmopt(dt)"
  puts "msdrange   = $cmopt(msdrange)"
  puts "mode       = $cmopt(mode)"

  puts ""

}

proc crd_as_array {c varname} {
  upvar $varname crd 
  
  set crd(0) [lindex $c 0]
  set crd(1) [lindex $c 1] 
  set crd(2) [lindex $c 2] 
}

proc cm_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global cmopt
  global common

  global sel 


  set anatra_path $::env(ANATRA_PATH);list
  set cmfort     "${anatra_path}/f90/bin/cm_analysis.x";list

  # read trajectory
  #
  puts ""
  puts "--------------------"
  puts " Read trajectory"
  puts "--------------------"
  puts ""

  set mol [mol load $str(stype) "$str(sfile)"]

  # setup selection
  #
  puts ""
  puts "--------------------"
  puts " Setup selection"
  puts "--------------------"
  puts ""
  for {set isel 0} {$isel < $seltxt(nsel)} {incr isel} {
    puts [format "selection %5d : %s" $isel $seltxt($isel)]
    set sel($isel) [atomselect $mol "$seltxt($isel)"]
  } 

  set rnam [$sel(0) get resname]
  set res  [$sel(0) get resid]
  set mass [$sel(0) get mass]
  set anam [$sel(0) get name]
  set chg  [$sel(0) get charge]
  set ind  [$sel(0) get index]
  set segn [$sel(0) get segname]
  set natm [llength $res]

  # Convert 
  #
  puts ""
  puts "--------------------"
  puts " Start analysis"
  puts "--------------------"
  puts ""
  puts ""

  puts ">> Start CoM calculation"
  puts ""

  set rand [expr int((100000*rand()))]
  set fmolinfo [format "cm%06d.molinfo" $rand]
  set fcminp   [format "cm%06d.inp"     $rand]
  set fcmout   [format "cm%06d.out"     $rand]
  set fdcdtmp  [format "cm%06d.dcd"     $rand]

  if {$common(prep_only)} {
    set fmolinfo [format "cm.molinfo"]
    set fcminp   [format "cm.inp"]
  }

  set f [open $fmolinfo "w"]
  for {set iatm 0} {$iatm < $natm} {incr iatm} {
    puts $f [format "%10d  %6s  %6s  %15.7f  %15.7f  %d  %6s  %3s" \
       [lindex $res  $iatm]           \
	     [lindex $rnam $iatm]           \
	     [lindex $anam $iatm]           \
	     [lindex $mass $iatm]           \
	     [lindex $chg  $iatm]           \
	     [expr [lindex $ind $iatm] + 1] \
       [lindex $segn $iatm]           \
       "END"]
  }
  close $f

  set ntraj [llength $traj(tin)] 

  set f [open $fcminp "w"]
  puts $f " &input_param"
  puts $f  "  ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  puts $f ""

  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead   = \"$cmopt(fhead)\""
  puts $f " /"

  puts $f " &trajopt_param"
  puts $f "   dt      = $cmopt(dt)"
  puts $f "   molinfo = \"$fmolinfo\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   mode     = \"$cmopt(mode)\""
  puts $f "   msdcalc  = .$cmopt(outmsd)."
  puts $f "   onlyz    = .$cmopt(onlyz)."
  puts $f "   unwrap   = .$cmopt(unwrap)."
  puts $f "   msddim   = $cmopt(msddim)"
  puts $f "   msdrange = $cmopt(msdrange)"
  puts $f " /"

  close $f

  #animate write dcd $fdcdtmp beg 0 end -1 waitfor all sel $sel(0) $mol

  if {!$common(prep_only)} {
    puts "CoM is calculated with ANATRA fortran program:"
    puts "$cmfort ..."
    puts "=== INPUT ==="
    set content [exec cat $fcminp]
    puts $content
    puts "============="
    exec $cmfort $fcminp >& $fcmout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fcmout]
    puts $content
    exec rm -f $fmolinfo $fcminp $fcmout $fdcdtmp
  }

  puts "=============="
  puts ">> Finished"

  exit
}
