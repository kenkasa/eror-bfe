proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                  CoM Coordinate Analysis"
  puts ""
  puts "============================================================"
}

proc show_tr_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra comcrd                                                              \\"
    puts "  -stype      <structure file type>                                        \\"
    puts "  -sfile      <structure file name>                                        \\"
    puts "  -tintype    <input trajectory file type>                                 \\"
    puts "  -tin        <input trajectory file name>                                 \\"
    puts "  -flist_traj <trajectory file list (neccesary if tin is not specified)>   \\"
    puts "  -fhead      <header of output file name>                                 \\"
    puts "  -flog       <log file name (optional)>                                   \\"
    puts "  -sel0       <VMD selection> (X=0,1,2...)                                 \\"
    puts "  -mode       <analysis mode (residue or whole or atom)>                   \\"
    puts "              (default: residue)                                           \\"
    puts "  -outcom     <whether com file is generated (true or false)>              \\"
    puts "              (default:false)>                                             \\"
    puts "  -outmsd     <whether msd file is generated (true or false)>              \\"
    puts "              (default:false)>                                             \\"
    puts "  -outngp     <whether ngp file is generated (true or false)>              \\"
    puts "              (default:false)>                                             \\"
    puts "  -outvhf     <whether vhf file is generated (true or false)>              \\"
    puts "              (default:false)>                                             \\"
    puts "  -unwrap     <whether unwrapping trajecotry is performed (true or false)> \\"
    puts "              (default:false)                                              \\"
    puts "  -msddim     <dimension used for MSD analysis (2 or 3)>                   \\"
    puts "              (default:3)                                                  \\"
    puts "  -dt         <time interval>                                              \\"
    puts "  -msdrange   <number of steps analyzed for MSD>                           \\"
    puts "  -t_sparse   <output time interval for VHF>                               \\"
    puts "  -t_range    <time range of VHF output>                                   \\"
    puts "  -dr         <distance interval for VHF>                                  \\"
    puts "  -rmax       <maximum distance for VHF>                                   \\"
    puts "  -prep_only  <where analysis is performed or not (true or false)          \\"
    puts "              (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra comcrd           \\"
    puts "  -stype     parm7      \\"
    puts "  -sfile     str.prmtop \\"
    puts "  -tintype   dcd        \\"
    puts "  -tin       inp.dcd    \\"
    puts "  -fhead     out        \\"
    puts "  -outcom    true       \\"
    puts "  -outmsd    true       \\"
    puts "  -outngp    false      \\"
    puts "  -outvhf    false      \\"
    puts "  -sel0      not water  \\"
    puts "  -mode      residue    \\"
    puts "  -unwrap    false      \\"
    puts "  -msddim    3          \\"
    puts "  -dt        0.1        \\"
    puts "  -msdrange  10         \\"
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
  set cmopt(flog)       ""
  set cmopt(fts)        ""
  set cmopt(outcom)     false
  set cmopt(outmsd)     false
  set cmopt(outngp)     false
  set cmopt(outvhf)     false 
  set cmopt(onlyz)      false
  set cmopt(use_cond)   false
  set cmopt(unwrap)     false
  set cmopt(msddim)     3
  set cmopt(dt)         0.1 
  set cmopt(msdrange)   10
  set cmopt(t_sparse)   -1.0
  set cmopt(t_range)    -1.0
  set cmopt(dr)         0.1
  set cmopt(rmax)       50.0
  set cmopt(rcrange)    "0.0 0.0" 
  set cmopt(mode)       "residue";list 
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
  set cmopt(flog)       [parse_arguments $arglist \
      "-flog"       "value" $cmopt(flog)]
  set cmopt(fts)        [parse_arguments $arglist \
      "-fts"        "value" $cmopt(fts)]
  set cmopt(mode)       [parse_arguments $arglist \
      "-mode"       "value" $cmopt(mode)]
  set cmopt(outcom)     [parse_arguments $arglist \
      "-outcom"     "value" $cmopt(outcom)]
  set cmopt(outmsd)     [parse_arguments $arglist \
      "-outmsd"     "value" $cmopt(outmsd)]
  set cmopt(outngp)     [parse_arguments $arglist \
      "-outngp"     "value" $cmopt(outngp)]
  set cmopt(outvhf)     [parse_arguments $arglist \
      "-outvhf"     "value" $cmopt(outvhf)]
  set cmopt(unwrap)     [parse_arguments $arglist \
      "-unwrap"     "value" $cmopt(unwrap)]
  set cmopt(onlyz)      [parse_arguments $arglist \
      "-onlyz"      "value" $cmopt(onlyz)]
  set cmopt(use_cond)   [parse_arguments $arglist \
      "-use_cond"   "value" $cmopt(use_cond)]
  set cmopt(msddim)     [parse_arguments $arglist \
      "-msddim"     "value" $cmopt(msddim)]
  set cmopt(dt)         [parse_arguments $arglist \
      "-dt"         "value" $cmopt(dt)]
  set cmopt(msdrange)   [parse_arguments $arglist \
      "-msdrange"   "value" $cmopt(msdrange)]
  set cmopt(t_sparse)   [parse_arguments $arglist \
      "-t_sparse"   "value" $cmopt(t_sparse)]
  set cmopt(t_range)    [parse_arguments $arglist \
      "-t_range"    "value" $cmopt(t_range)]
  set cmopt(dr)         [parse_arguments $arglist \
      "-dr"         "value" $cmopt(dr)]
  set cmopt(rmax)       [parse_arguments $arglist \
      "-rmax"       "value" $cmopt(rmax)]
  set cmopt(rcrange)    [parse_arguments $arglist \
      "-rcrange"    "value" $cmopt(rcrange)]
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
  puts "fts        = $cmopt(fts)"
  puts "outcom     = $cmopt(outcom)"
  puts "outmsd     = $cmopt(outmsd)"
  puts "outngp     = $cmopt(outngp)"
  puts "outvhf     = $cmopt(outvhf)"
  puts "onlyz      = $cmopt(onlyz)"
  puts "use_cond   = $cmopt(use_cond)"
  puts "msddim     = $cmopt(msddim)"
  puts "dt         = $cmopt(dt)"
  puts "msdrange   = $cmopt(msdrange)"
  puts "t_sparse   = $cmopt(t_sparse)"
  puts "t_range    = $cmopt(t_range)"
  puts "dr         = $cmopt(dr)"
  puts "rmax       = $cmopt(rmax)"
  puts "mode       = $cmopt(mode)"
  puts "rcrange    = $cmopt(rcrange)"

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
  set cmfort     "${anatra_path}/f90/bin/comcrd_analysis.x";list

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

  if {$cmopt(flog) == ""} { 
    set fcmout   [format "cm%06d.out"     $rand]
  } else {
    set fcmout   $cmopt(flog)
  }

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
  puts $f  "  fts        = \"$cmopt(fts)\""
  puts $f  "  flist_traj = \"$traj(flist_traj)\""
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
  puts $f "   comcalc  = .$cmopt(outcom)."
  puts $f "   msdcalc  = .$cmopt(outmsd)."
  puts $f "   ngpcalc  = .$cmopt(outngp)."
  puts $f "   vhfcalc  = .$cmopt(outvhf)."
  puts $f "   onlyz    = .$cmopt(onlyz)."
  puts $f "   use_cond = .$cmopt(use_cond)"
  puts $f "   unwrap   = .$cmopt(unwrap)."
  puts $f "   msddim   = $cmopt(msddim)"
  puts $f "   msdrange = $cmopt(msdrange)"
  puts $f "   dt       = $cmopt(dt)"
  puts $f "   t_sparse = $cmopt(t_sparse)"
  puts $f "   t_range  = $cmopt(t_range)"
  puts $f "   dr       = $cmopt(dr)"
  puts $f "   rmax     = $cmopt(rmax)"
  puts $f "   rcrange  = $cmopt(rcrange)"
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

    if {$cmopt(flog) == ""} {
      exec rm -f $fmolinfo $fcminp $fcmout $fdcdtmp
    } else {
      exec rm -f $fmolinfo $fcminp $fdcdtmp
    }
  }

  puts "=============="
  puts ">> Finished"

  exit
}
