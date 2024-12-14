proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                  Orientation Analysis"
  puts ""
  puts "============================================================"
}

proc show_or_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra or \\"
    puts "  -stype       <structure file type>                                \\"
    puts "  -sfile       <structure file name>                                \\"
    puts "  -tintype     <input trajectory file type>                         \\"
    puts "  -tin         <input trajectory file name>                         \\"
    puts "  -fhead       <header of output file name>                         \\"
    puts "  -flog        <log file name (optional)>                           \\"
    puts "  -judgeup     <how to judge up region (nmolup or coord)>           \\"
    puts "               (default: coord)                                     \\"
    puts "  -nmolup      <number of molecules in upper leaflet>               \\"
    puts "  -nmollow     <number of molecules in lower leaflet>               \\"
    puts "  -xcoord      <coordinate for distribution (cost or theta)         \\"
    puts "               (default: cost)>                                     \\"
    puts "  -dtheta      <theta grid>                                         \\"
    puts "  -dt          <time interval>                                      \\"
    puts "  -sel0        <VMD selection>                                      \\"
    puts "  -sel1        <VMD selection>                                      \\"
    puts "  -mode        <analysis mode (residue or whole or segment)         \\"
    puts "               (default: residue)>                                  \\"
    puts "  -prep_only   <where analysis is performed or not (true or false)> \\" 
    puts "               (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra or                 \\"
    puts "  -stype       parm7      \\"
    puts "  -sfile       str.prmtop \\"
    puts "  -tintype     dcd        \\"
    puts "  -tin         inp.dcd    \\"
    puts "  -fhead       output     \\"
    puts "  -judgeup     coord      \\"
    puts "  -xcoord      cost       \\"
    puts "  -dtheta      1.0        \\"
    puts "  -dt          0.1        \\"
    puts "  -sel0        name P     \\"
    puts "  -sel1        name N     \\"
    puts "  -mode        residue    \\"
    puts "  -prep_only   false"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_oroptinfo
#! @brief         define Orient option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_oroptinfo {} {
  
  global oropt

  set oropt(fhead)      "output"
  set oropt(flog)       ""
  set oropt(judgeup)    "coord"
  set oropt(nmolup)     100
  set oropt(nmollow)    0
  set oropt(xcoord)     "cos"
  set oropt(dtheta)     1.0 
  set oropt(dcost)      0.1 
  set oropt(dt)         0.1 
  set oropt(mode)       "residue" 
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_oroptinfo
#! @brief         read Dipole option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_oroptinfo {arglist} {

  global oropt

  set oropt(fhead)      [parse_arguments $arglist \
      "-fhead"        "value" $oropt(fhead)]
  set oropt(flog)       [parse_arguments $arglist \
      "-flog"         "value" $oropt(flog)]
  set oropt(judgeup)    [parse_arguments $arglist \
      "-judgeup"      "value" $oropt(judgeup)]
  set oropt(nmolup)     [parse_arguments $arglist \
      "-nmolup"       "value" $oropt(nmolup)]
  set oropt(nmollow)    [parse_arguments $arglist \
      "-nmollow"      "value" $oropt(nmollow)]
  set oropt(mode)        [parse_arguments $arglist \
      "-mode"         "value" $oropt(mode)]
  set oropt(xcoord)      [parse_arguments $arglist \
      "-xcoord"       "value" $oropt(xcoord)]
  set oropt(dtheta)      [parse_arguments $arglist \
      "-dtheta"       "value" $oropt(dtheta)]
  set oropt(dcost)       [parse_arguments $arglist \
      "-dcost"       "value" $oropt(dcost)]
  set oropt(dt)          [parse_arguments $arglist \
      "-dt"           "value" $oropt(dt)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_oroptinfo
#! @brief         show Dipole option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_oroptinfo {} {

  global oropt

  puts "<< option info >>"
  puts "fhead      = $oropt(fhead)"
  puts "judgeup    = $oropt(judgeup)"
  puts "nmolup     = $oropt(nmolup)"
  puts "nmollow    = $oropt(nmollow)"
  puts "xcoord     = $oropt(xcoord)"
  puts "dtheta     = $oropt(dtheta)"
  puts "dt         = $oropt(dt)"
  puts "mode       = $oropt(mode)"
  puts ""

}

proc or_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global oropt
  global common

  global sel 


  set anatra_path $::env(ANATRA_PATH);list
  set orfort     "${anatra_path}/f90/bin/or_analysis.x";list

  # read trajectory
  #
  puts ""
  puts "--------------------"
  puts " Read trajectory"
  puts "--------------------"
  puts ""

  #set mol 0;
  #read_traj $mol $str(stype) $str(sfile) $traj(tintype) $traj(tin) $traj(stride)
  #set nf   [molinfo $mol get numframes]
  
  set mol [mol load $str(stype) "$str(sfile)"]
  
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

  set nf [molinfo $mol get numframes] 
  puts ""

  puts ">> Start CoM calculation"
  puts ""

  set rand [expr int((100000*rand()))]
  set forinp   [format "or%06d.inp"     $rand]

  if {$oropt(flog) == ""} {
    set forout   [format "or%06d.out"     $rand]
  } else {
    set forout   $oropt(flog)
  }

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fdcdtmp($isel)  [format "or%06d_%i.dcd"     $rand $isel]
    set fmolinfo($isel) [format "or%06d_%i.molinfo" $rand $isel]
  }

  if {$common(prep_only)} {
    for {set isel 0} {$isel < $nsel} {incr isel} {
      set fmolinfo($isel) [format "or_%i.molinfo" $isel]
    }
    set forinp   [format "or.inp"]
  }

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set rnam [$sel($isel) get resname]
    set res  [$sel($isel) get resid]
    set mass [$sel($isel) get mass]
    set anam [$sel($isel) get name]
    set chg  [$sel($isel) get charge]
    set ind  [$sel($isel) get index]
    set segn [$sel($isel) get segname]
    set natm [llength $res]

    set f [open $fmolinfo($isel) "w"]
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

    #animate write dcd $fdcdtmp($isel) \
    #  beg 0 end -1 waitfor all sel $sel($isel) $mol 
  }

  set ntraj [llength $traj(tin)] 


  set f [open $forinp "w"]
  puts $f " &input_param"
  puts $f "   ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  puts $f ""

  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead   = \"$oropt(fhead)\""
  puts $f " /"

  puts $f " &trajopt_param"
  puts $f "   dt      = $oropt(dt)"
  puts $f "   molinfo = \"$fmolinfo(0)\" \"$fmolinfo(1)\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   judgeup = \"$oropt(judgeup)\""
  puts $f "   mode    = \"$oropt(mode)\" \"$oropt(mode)\""
  puts $f "   xcoord  = \"$oropt(xcoord)\""
  puts $f "   dcost   = $oropt(dcost)"
  puts $f "   dtheta  = $oropt(dtheta)"
  puts $f "   nmolup  = $oropt(nmolup)"
  puts $f "   nmollow = $oropt(nmollow)"
  puts $f " /"

  close $f


  if {!$common(prep_only)} {
    puts "Orientation is calculated with ANATRA fortran program:"
    puts "$orfort ..."
    puts "=== INPUT ==="
    set content [exec cat $forinp]
    puts $content
    puts "============="
    exec $orfort $forinp >& $forout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $forout]
    puts $content

    if {$oropt(flog) == ""} {
      exec rm -f $fmolinfo(0) $fmolinfo(1) $forinp $forout
    } else {
      exec rm -f $fmolinfo(0) $fmolinfo(1) $forinp
    }
  }

  puts "=============="
  puts ">> Finished"

  exit
}
