proc print_title {} {
  puts "============================================================"
  puts ""
  puts "              Electric Conductivity Analysis"
  puts ""
  puts "============================================================"
}

proc show_el_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra elecond                                                             \\"
    puts "  -stype      <structure file type>                                        \\"
    puts "  -sfile      <structure file name>                                        \\"
    puts "  -tintype    <input trajectory file type>                                 \\"
    puts "  -tin        <input trajectory file name>                                 \\"
    puts "  -flist_traj <trajectory file list (neccesary if tin is not specified)>   \\"
    puts "  -fhead      <header of output file name>                                 \\"
    puts "  -dt         <time interval>                                              \\"
    puts "  -calctype   <VEL or COORD>                                               \\"
    puts "  -tcfrange   <number of steps analyzed for TCF>                           \\"
    puts "  -bigtraj    <true or false>                                              \\"
    puts "  -prep_only  <where analysis is performed or not (true or false)          \\"
    puts "              (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra elecond          \\"
    puts "  -stype     parm7      \\"
    puts "  -sfile     str.prmtop \\"
    puts "  -tintype   dcd        \\"
    puts "  -tin       inp.dcd    \\"
    puts "  -fhead     out        \\"
    puts "  -calctype  COORD      \\"
    puts "  -dt        0.1        \\"
    puts "  -tcfrange  10         \\"
    puts "  -prep_only false"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_eloptinfo
#! @brief         define CoM option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_eloptinfo {} {
  
  global elopt

  set elopt(fhead)      "out"
  set elopt(dt)         0.1
  set elopt(calctype)   "COORD"
  set elopt(tcfrange)   10
  set elopt(t_sparse)   -1.0
  set elopt(t_range)    -1.0
  set elopt(bigtraj)    "false";list 
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_eloptinfo
#! @brief         read CoM option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_eloptinfo {arglist} {

  global elopt

  set elopt(fhead)      [parse_arguments $arglist \
      "-fhead"      "value" $elopt(fhead)]
  set elopt(calctype)   [parse_arguments $arglist \
      "-calctype"   "value" $elopt(calctype)]
  set elopt(dt)         [parse_arguments $arglist \
      "-dt"         "value" $elopt(dt)]
  set elopt(tcfrange)   [parse_arguments $arglist \
      "-tcfrange"   "value" $elopt(tcfrange)]
  set elopt(t_sparse)   [parse_arguments $arglist \
      "-t_sparse"   "value" $elopt(t_sparse)]
  set elopt(t_range)    [parse_arguments $arglist \
      "-t_range"    "value" $elopt(t_range)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_eloptinfo
#! @brief         show CoM option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_eloptinfo {} {

  global elopt

  puts "<< option info >>"
  puts "fhead      = $elopt(fhead)"
  puts "dt         = $elopt(dt)"
  puts "calctype   = $elopt(calctype)"
  puts "tcfrange   = $elopt(tcfrange)"
  puts "bigtraj    = $elopt(t_sparse)"
  puts "t_sparse   = $elopt(t_sparse)"
  puts "t_range    = $elopt(t_range)"

  puts ""

}

proc crd_as_array {c varname} {
  upvar $varname crd 
  
  set crd(0) [lindex $c 0]
  set crd(1) [lindex $c 1] 
  set crd(2) [lindex $c 2] 
}

proc el_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global elopt
  global common

  global sel 


  set anatra_path $::env(ANATRA_PATH);list
  set elfort     "${anatra_path}/f90/bin/elecond_analysis.x";list

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
  #puts ""
  #puts "--------------------"
  #puts " Setup selection"
  #puts "--------------------"
  #puts ""
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

  puts ">> Start EleCond calculation"
  puts ""

  set rand [expr int((100000*rand()))]
  set fmolinfo [format "el%06d.molinfo" $rand]
  set felinp   [format "el%06d.inp"     $rand]

  #if {$elopt(flog) == ""} { 
    set felout   [format "el%06d.out"     $rand]
  #} else {
  #  set felout   $elopt(flog)
  #}

  set fdcdtmp  [format "el%06d.dcd"     $rand]

  if {$common(prep_only)} {
    set fmolinfo [format "el.molinfo"]
    set felinp   [format "el.inp"]
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

  set f [open $felinp "w"]
  puts $f " &input_param"
  puts $f  "  flist_traj = \"$traj(flist_traj)\""
  puts $f  "  ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  puts $f ""

  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead   = \"$elopt(fhead)\""
  puts $f " /"

  puts $f " &trajopt_param"
  puts $f "   dt      = $elopt(dt)"
  puts $f "   molinfo = \"$fmolinfo\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   dt       = $elopt(dt)"
  puts $f "   calctype = $elopt(calctype)"
  puts $f "   tcfrange = $elopt(tcfrange)"
  puts $f "   bigtraj  = $elopt(bigtraj)"
  puts $f "   t_sparse = $elopt(t_sparse)"
  puts $f "   t_range  = $elopt(t_range)"
  puts $f " /"

  close $f

  #animate write dcd $fdcdtmp beg 0 end -1 waitfor all sel $sel(0) $mol

  if {!$common(prep_only)} {
    puts "EleCond is calculated with ANATRA fortran program:"
    puts "$elfort ..."
    puts "=== INPUT ==="
    set content [exec cat $felinp]
    puts $content
    puts "============="
    exec $elfort $felinp >& $felout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $felout]
    puts $content

    #if {$elopt(flog) == ""} {
    #  exec rm -f $fmolinfo $fcminp $fcmout $fdcdtmp
    #} else {
      exec rm -f $fmolinfo $felinp $felout $fdcdtmp 
    #}
  }

  puts "=============="
  puts ">> Finished"

  exit
}
