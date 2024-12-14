proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                    Z-profile Analysis"
  puts ""
  puts "============================================================"
}

proc show_zp_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra zp                                                          \\"
    puts "  -stype      <structure file type>                                \\"
    puts "  -sfile      <structure file name>                                \\"
    puts "  -tintype    <input trajectory file type>                         \\"
    puts "  -tin        <input trajectory file name>                         \\"
    puts "  -zpfile     <output z-profile file name>                         \\"
    puts "  -sel0       <VMD selection> (X=0,1,2...)                         \\"
    puts "  -mode       <analysis mode (residue or whole or atom)>           \\" 
    puts "              (default: residue)                                   \\"
    puts "  -denstype   <density type (number or electron)>                  \\"
    puts "              (default: number)                                    \\"
    puts "  -centertype <box center position (zero or half)>                 \\"
    puts "              (default: zero)                                      \\"
    puts "  -prep_only  <where analysis is performed or not (true or false)> \\"
    puts "              (default: false)>"
    puts ""
    puts "Usage:"
    puts "anatra zp                \\"
    puts "  -stype      parm7      \\"
    puts "  -sfile      str.prmtop \\"
    puts "  -tintype    dcd        \\"
    puts "  -tin        inp.dcd    \\"
    puts "  -zpfile     out.zp     \\"
    puts "  -sel0       water      \\"
    puts "  -mode       residue    \\"
    puts "  -denstype   number     \\"
    puts "  -centertype zero       \\"
    puts "  -prep_only  false"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_zpoptinfo
#! @brief         define CoM option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_zpoptinfo {} {
  
  global zpopt

  set zpopt(zpfile)      "out.zp" 
  set zpopt(mode)        "residue"
  set zpopt(dz)          0.1
  set zpopt(denstype)    "number" 
  set zpopt(centertype)  "zero" 
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_zpoptinfo
#! @brief         read z-profile option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_zpoptinfo {arglist} {

  global zpopt

  set zpopt(zpfile)      [parse_arguments $arglist \
      "-zpfile"     "value" $zpopt(zpfile)]
  set zpopt(mode)        [parse_arguments $arglist \
      "-mode"       "value" $zpopt(mode)]
  set zpopt(dz)          [parse_arguments $arglist \
      "-dz"         "value" $zpopt(dz)]
  set zpopt(denstype)    [parse_arguments $arglist \
      "-denstype"   "value" $zpopt(denstype)]
  set zpopt(centertype)  [parse_arguments $arglist \
      "-centertype" "value" $zpopt(centertype)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_zpoptinfo
#! @brief         show Z-profile option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_zpoptinfo {} {

  global zpopt

  puts "<< option info >>"
  puts "zpfile      = $zpopt(zpfile)"
  puts "mode        = $zpopt(mode)"
  puts "denstype    = $zpopt(denstype)"
  puts "centertype  = $zpopt(centertype)"
  puts ""

}

proc zp_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global zpopt
  global common

  global sel 


  set anatra_path $::env(ANATRA_PATH);list
  set zpfort     "${anatra_path}/f90/bin/zp_analysis.x";list

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

  set nf [molinfo $mol get numframes] 
  puts ""

  puts ">> Start Z-profile calculation"
  puts ""

  set rand [expr int((100000*rand()))]
  set fmolinfo [format "zp%06d.molinfo" $rand]
  set fzpinp   [format "zp%06d.inp"     $rand]
  set fzpout   [format "zp%06d.out"     $rand]
  set fdcdtmp  [format "zp%06d.dcd"     $rand]

  if {$common(prep_only)} {
    set fmolinfo [format "zp.molinfo"]
    set fzpinp   [format "zp.inp"]
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

  set f [open $fzpinp "w"]
  puts $f " &input_param"
  puts $f  "  ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  puts $f ""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fzp    = \"$zpopt(zpfile)\""
  puts $f " /"

  puts $f " &trajopt_param"
  puts $f "   molinfo = \"$fmolinfo\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   mode       = \"$zpopt(mode)\""
  puts $f "   dz         = $zpopt(dz)"
  puts $f "   denstype   = \"$zpopt(denstype)\""
  puts $f "   centertype = \"$zpopt(centertype)\""
  puts $f " /"

  close $f

  #animate write dcd $fdcdtmp beg 0 end -1 waitfor all sel $sel(0) $mol

  if {!$common(prep_only)} {
    puts "Z-profile is calculated with ANATRA fortran program:"
    puts "$zpfort ..."
    puts "=== INPUT ==="
    set content [exec cat $fzpinp]
    puts $content
    puts "============="
    exec $zpfort $fzpinp > $fzpout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fzpout]
    puts $content
    exec rm -f $fmolinfo $fzpinp $fzpout
  }
  puts "=============="
  puts ">> Finished"

  exit
}
