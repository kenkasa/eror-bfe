proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                   Z-deviation Analysis"
  puts ""
  puts "============================================================"
}

proc show_zd_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra z_deviation \\"
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
    puts "  -ngrid       <number of grids>                                    \\"
    puts "  -dz          <z grid>                                             \\"
    puts "  -dt          <time interval>                                      \\"
    puts "  -sel0        <VMD selection>                                      \\"
    puts "  -mode        <analysis mode (residue or whole or segment)         \\"
    puts "               (default: residue)>                                  \\"
    puts "  -prep_only   <where analysis is performed or not (true or false)> \\" 
    puts "               (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra z_deviation        \\"
    puts "  -stype       parm7      \\"
    puts "  -sfile       str.prmtop \\"
    puts "  -tintype     dcd        \\"
    puts "  -tin         inp.dcd    \\"
    puts "  -fhead       output     \\"
    puts "  -judgeup     coord      \\"
    puts "  -ngrid       100        \\"
    puts "  -dz          1.0        \\"
    puts "  -dt          0.1        \\"
    puts "  -sel0        name P     \\"
    puts "  -mode        residue    \\"
    puts "  -prep_only   false"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_zdoptinfo
#! @brief         define zdev option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_zdoptinfo {} {
  
  global zdopt

  set zdopt(fhead)      "output"
  set zdopt(flog)       ""
  set zdopt(judgeup)    "coord"
  set zdopt(nmolup)     100
  set zdopt(nmollow)    0
  set zdopt(ngrid)      100.0 
  set zdopt(dz)         1.0 
  set zdopt(dt)         0.1 
  set zdopt(mode)       "residue" 
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_zdoptinfo
#! @brief         read zdev option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_zdoptinfo {arglist} {

  global zdopt

  set zdopt(fhead)      [parse_arguments $arglist \
      "-fhead"        "value" $zdopt(fhead)]
  set zdopt(flog)       [parse_arguments $arglist \
      "-flog"         "value" $zdopt(flog)]
  set zdopt(judgeup)    [parse_arguments $arglist \
      "-judgeup"      "value" $zdopt(judgeup)]
  set zdopt(nmolup)     [parse_arguments $arglist \
      "-nmolup"       "value" $zdopt(nmolup)]
  set zdopt(nmollow)    [parse_arguments $arglist \
      "-nmollow"      "value" $zdopt(nmollow)]
  set zdopt(mode)        [parse_arguments $arglist \
      "-mode"         "value" $zdopt(mode)]
  set zdopt(ngrid)       [parse_arguments $arglist \
      "-ngrid"        "value" $zdopt(ngrid)]
  set zdopt(dz)          [parse_arguments $arglist \
      "-dz"           "value" $zdopt(dz)]
  set zdopt(dt)          [parse_arguments $arglist \
      "-dt"           "value" $zdopt(dt)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_zdoptinfo
#! @brief         show Dipole option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_zdoptinfo {} {

  global zdopt

  puts "<< option info >>"
  puts "fhead      = $zdopt(fhead)"
  puts "judgeup    = $zdopt(judgeup)"
  puts "nmolup     = $zdopt(nmolup)"
  puts "nmollow    = $zdopt(nmollow)"
  puts "ngrid      = $zdopt(ngrid)"
  puts "dz         = $zdopt(dz)"
  puts "dt         = $zdopt(dt)"
  puts "mode       = $zdopt(mode)"
  puts ""

}

proc zd_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global zdopt
  global common

  global sel 


  set anatra_path $::env(ANATRA_PATH);list
  set zdfort     "${anatra_path}/f90/bin/zdev_analysis.x";list

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
  set fzdinp   [format "zd%06d.inp"     $rand]

  if {$zdopt(flog) == ""} {
    set fzdout   [format "zd%06d.out"     $rand]
  } else {
    set fzdout   $zdopt(flog)
  }

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fdcdtmp($isel)  [format "zd%06d_%i.dcd"     $rand $isel]
    set fmolinfo($isel) [format "zd%06d_%i.molinfo" $rand $isel]
  }

  if {$common(prep_only)} {
    for {set isel 0} {$isel < $nsel} {incr isel} {
      set fmolinfo($isel) [format "zd_%i.molinfo" $isel]
    }
    set fzdinp   [format "zd.inp"]
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


  set f [open $fzdinp "w"]
  puts $f " &input_param"
  puts $f "   ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  puts $f ""

  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead   = \"$zdopt(fhead)\""
  puts $f " /"

  puts $f " &trajopt_param"
  puts $f "   dt      = $zdopt(dt)"
  puts $f "   molinfo = \"$fmolinfo(0)\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   judgeup = \"$zdopt(judgeup)\""
  puts $f "   mode    = \"$zdopt(mode)\" \"$zdopt(mode)\""
  puts $f "   ngrid   = $zdopt(ngrid)"
  puts $f "   dz      = $zdopt(dz)"
  puts $f "   nmolup  = $zdopt(nmolup)"
  puts $f "   nmollow = $zdopt(nmollow)"
  puts $f " /"

  close $f


  if {!$common(prep_only)} {
    puts "z-deviation is calculated with ANATRA fortran program:"
    puts "$zdfort ..."
    puts "=== INPUT ==="
    set content [exec cat $fzdinp]
    puts $content
    puts "============="
    exec $zdfort $fzdinp >& $fzdout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fzdout]
    puts $content

    if {$zdopt(flog) == ""} {
      exec rm -f $fmolinfo(0) $fzdinp $fzdout
    } else {
      exec rm -f $fmolinfo(0) $fzdinp
    }
  }

  puts "=============="
  puts ">> Finished"

  exit
}
