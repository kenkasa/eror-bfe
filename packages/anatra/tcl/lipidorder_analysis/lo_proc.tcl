proc print_title {} {
  puts "============================================================"
  puts ""
  puts "               Lipid Order Parameter Analysis"
  puts ""
  puts "============================================================"
  puts "Remark: Current version only supports CHARMM forcefield"
}

proc show_sc_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra sc                                                         \\"
    puts "  -stype     <structure file type>                                \\"
    puts "  -sfile     <structure file name>                                \\"
    puts "  -tintype   <input trajectory file type>                         \\"
    puts "  -tin       <input trajectory file name>                         \\"
    puts "  -fhead     <header of output file name>                         \\"
    puts "  -dt        <time interval>                                      \\"
    puts "  -cindex    <index for analyzed carbons>                         \\"
    puts "  -sel0      <VMD selection> (X=0,1,2...)                         \\"
    puts "  -prep_only <where analysis is performed or not (true or false)> \\"
    puts "             (default: false)>"
    puts ""  
    puts "Usage:"
    puts "anatra sc                                     \\"
    puts "  -stype     parm7                            \\"
    puts "  -sfile     str.prmtop                       \\"
    puts "  -tintype   dcd                              \\"
    puts "  -tin       inp.dcd                          \\"
    puts "  -fhead     scd                              \\"
    puts "  -dt        0.1                              \\"
    puts "  -cindex    1 2 3                            \\"
    puts "  -sel0      name C32  H2X H2Y and segid MEMB \\"
    puts "  -sel1      name C33  H3X H3Y and segid MEMB \\"
    puts "  -sel2      name C34  H4X H4Y and segid MEMB \\"
    puts "  -prep_only false"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_scoptinfo
#! @brief         define Scd option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_scoptinfo {} {
  
  global scopt

  set scopt(fhead)    "out"
  set scopt(dt)        0.01
  set scopt(cindex)    " 1  2  3  4  5  6  7  8  9 10 \
                        11 12 13 14 15 16 17 18 19 20"
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_scoptinfo
#! @brief         read Scd option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_scoptinfo {arglist} {

  global scopt

  set scopt(fhead) [parse_arguments $arglist \
      "-fhead"     "value" $scopt(fhead)]
  set scopt(dt)     [parse_arguments $arglist \
      "-dt"         "value" $scopt(dt)]
  set scopt(cindex) [parse_arguments $arglist \
      "-cindex"     "value" $scopt(cindex)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_scoptinfo
#! @brief         show Scd option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_scoptinfo {} {

  global scopt

  puts "<< option info >>"
  puts "fhead      = $scopt(fhead)"
  puts "dt         = $scopt(dt)"
  puts "cindex     = $scopt(cindex)"
  puts ""

}

proc sc_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global scopt
  global common

  global sel 

  set anatra_path $::env(ANATRA_PATH);list
  set scfort     "${anatra_path}/f90/bin/lipidorder_analysis.x";list

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

  set nsel $seltxt(nsel)
  set rand [expr int((100000*rand()))]
  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fmolinfo($isel) [format "sc_%i_%06d.molinfo" $isel $rand]
  }
  set fscinp   [format "sc%06d.inp"     $rand]
  set fscout   [format "sc%06d.out"     $rand]
  set fdcdtmp  [format "sc%06d.dcd"     $rand]

  if {$common(prep_only)} {
    for {set isel 0} {$isel < $nsel} {incr isel} {
      set fmolinfo($isel) [format "sc_%i.molinfo" $isel]
    }
    set fscinp   [format "sc.inp"]
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
    set nf   [molinfo $mol get numframes]
    set nres [llength [lsort -unique [$sel($isel) get residue]]]

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
  }

  #set i2        [format "%04d" $isel]
  #set fscd      [format "%s_%s.scd"    $scopt(header) $i2]
  #set fscdave   [format "%s_%s.scdave" $scopt(header) $i2]

  #set carbon_id [lindex $scopt(cindex) $isel]

  set ntraj   [llength $traj(tin)]

  set f [open $fscinp "w"]
  puts $f " &input_param"
  puts $f "   ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  puts $f ""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead   = \"$scopt(fhead)\""
  puts $f " /"
  
  puts $f " &trajopt_param"
  puts $f "   dt      = $scopt(dt)"
  puts $f "   molinfo = "
  for {set i 0} {$i < $nsel} {incr i} {
    puts -nonewline $f "    \"$fmolinfo($i)\" "
  }
  puts $f ""
  puts $f " /"
  
  puts $f " &option_param"
  puts $f "   carbon_id = $scopt(cindex)"
  puts $f " /"
  close $f

  #animate write dcd $fdcdtmp beg 0 end -1 waitfor all sel $sel($isel) $mol
  #

  if {!$common(prep_only)} {
    puts "Scd is calculated with ANATRA fortran program:"
    puts "$scfort ..."
    puts "=== INPUT ==="
    set content [exec cat $fscinp]
    puts $content
    puts "============="
    exec $scfort $fscinp >& $fscout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fscout]
    puts $content
    exec rm -f $fscinp $fscout
    for {set i 0} {$i < $nsel} {incr i} {
      exec rm -f  $fmolinfo($i)
    }
  }

  puts "=============="
  puts ">> Finished"

  exit
}


