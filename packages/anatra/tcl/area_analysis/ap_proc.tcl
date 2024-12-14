proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                  Area Per Lipid Analysis"
  puts ""
  puts "============================================================"
}

proc show_ap_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra ap -stype      <structure file type>          \\"
    puts "          -sfile      <structure file name>          \\"
    puts "          -tintype    <input trajectory file type>   \\"
    puts "          -tin        <input trajectory file name>   \\"
    puts "          -apfile     <area per lipid file name>     \\"
    puts "          -apdistfile <area per lipid distribution file name> \\"
    puts "          -flog       <log file name (optional)>              \\"
    puts "          -dt         <time interval> \\"
    puts "          -da         <delta area per lipid>"
    puts "          -sel0       <VMD selection> (X=0,1,2...)"
    puts ""
    puts "Usage:"
    puts "anatra ap -stype      parm7             \\"
    puts "          -sfile      str.prmtop        \\"
    puts "          -tintype    dcd               \\"
    puts "          -tin        inp.dcd           \\"
    puts "          -apfile     out.ap            \\"
    puts "          -apdistfile out.apdist        \\"
    puts "          -dt         0.01              \\"
    puts "          -da         0.5               \\"
    puts "          -sel0       segid MEMB"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_apoptinfo
#! @brief         define APL option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_apoptinfo {} {
  
  global apopt

  set apopt(apfile)        "out.ap" 
  set apopt(apdistfile)    "out.ap" 
  set apopt(flog)          "" 
  set apopt(dt)             0.01
  set apopt(da)             0.50 
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_apoptinfo
#! @brief         read APL option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_apoptinfo {arglist} {

  global apopt

  set apopt(apfile)     [parse_arguments $arglist \
      "-apfile"      "value" $apopt(apfile)]
  set apopt(apdistfile) [parse_arguments $arglist \
      "-apdistfile"  "value" $apopt(apdistfile)]
  set apopt(flog)       [parse_arguments $arglist \
      "-flog"        "value" $apopt(flog)]
  set apopt(dt)         [parse_arguments $arglist \
      "-dt"          "value" $apopt(dt)]
  set apopt(da)         [parse_arguments $arglist \
      "-da"          "value" $apopt(da)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_apoptinfo
#! @brief         show APL option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_apoptinfo {} {

  global apopt

  puts "<< option info >>"
  puts "apfile     = $apopt(apfile)"
  puts "apdistfile = $apopt(apdistfile)"
  puts "dt         = $apopt(dt)"
  puts "da         = $apopt(da)"
  puts ""

}

proc ap_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global apopt

  global sel 


  set anatra_path $::env(ANATRA_PATH);list
  set hsfort     "${anatra_path}/f90/bin/hs_analysis.x";list

  # read trajectory
  #
  puts ""
  puts "--------------------"
  puts " Read trajectory"
  puts "--------------------"
  puts ""

  set mol 0;
  read_traj $mol $str(stype) $str(sfile) $traj(tintype) $traj(tin) $traj(stride)
  set nf   [molinfo $mol get numframes]
  set nsel $seltxt(nsel)

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

  # Convert 
  #
  puts ""
  puts "--------------------"
  puts " Start analysis"
  puts "--------------------"
  puts ""

  set rand     [expr int((100000*rand()))]
  set fapinp   [format "ap%06d.inp"     $rand]

  if {$apopt(flog) == ""} {
    set fapout   [format "ap%06d.out"     $rand]
  } else {
    set fapout   $apopt(flog)
  }

  set nres [llength [lsort -unique [$sel(0) get resid]]]
  package require pbctools
  set box [pbc get -all]
  for {set istep 0} {$istep < $nf} {incr istep} {
    set bnow [lindex $box $istep]
    set apl($istep) [expr [lindex $bnow 0] * [lindex $bnow 1] * 2 / $nres ] 
  }

  set f [open $apopt(apfile) "w"]
  for {set istep 0} {$istep < $nf} {incr istep} {
    set t [expr $istep * $apopt(dt)]
    puts $f [format "%20.10f  %20.10f" $t $apl($istep)]
  }
  close $f

  set f [open $fapinp "w"]
  puts $f " &input_param"
  puts $f "   fts        = \"$apopt(apfile)\""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fdist      = \"$apopt(apdistfile)\""
  puts $f " /"
  puts $f " &option_param"
  puts $f "   dx       = $apopt(da)"
  puts $f "   xsta     = 30.0"
  puts $f "   nx       = 1000"
  puts $f " /"

  close $f

  puts "APL is calculated with ANATRA fortran program:"
  puts "$hsfort ..."
  puts "=== INPUT ==="
  set  content [exec cat $fapinp]
  puts $content
  puts "============="
  exec $hsfort $fapinp >& $fapout
  puts ""
  puts "=== OUTPUT ==="
  set  content [exec cat $fapout]
  puts $content
  puts "=============="
  puts ">> Finished"

  if {$apopt(flog) == ""} {
    exec rm -f $fapinp $fapout
  } else {
    exec rm -f $fapinp 
  }
  exit
}
