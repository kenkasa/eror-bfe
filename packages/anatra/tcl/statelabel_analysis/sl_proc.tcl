proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                 State Labelling Analysis"
  puts ""
  puts "============================================================"
}

proc show_sl_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra sl -stype      <structure file type>         \\"
    puts "          -sfile      <structure file name>         \\"
    puts "          -tintype    <input trajectory file type>  \\"
    puts "          -tin        <input trajectory file name>  \\"
    puts "          -fhead      <header of output file name>  \\"
    puts "          -mode0      <analysis mode of sel0 (residue or whole or atom) (default: residue)>  \\"
    puts "          -mode1      <analysis mode of sel1 (residue or whole or atom) (default: residue)>  \\"
    puts "          -dr         <delta r (angstrom)>          \\"
    puts "          -identical  <true or false>               \\"
    puts "          -normalize  <true or false>               \\"
    puts "          -separate_self <true or false> (default: false) \\"
    puts "          -use_conditional  <true or false>         \\"
    puts "          -use_symmetric    <true or false>         \\"
    puts "          -cond_type  <X,Y,Z> (default: Z)          \\"
    puts "          -cond_range  <min  max>                   \\"
    puts "          -bond_range  <bond range (angstrom)>      \\"
    puts "          -sel0       <VMD selection> (X=0,1,2...)  \\"
    puts "          -sel1       <VMD selection> (X=0,1,2...)  \\"
    puts ""  
    puts "Usage:"
    puts "anatra sl -stype     parm7                          \\"
    puts "          -sfile     str.prmtop                     \\"
    puts "          -tintype   dcd                            \\"
    puts "          -tin       inp.dcd                        \\"
    puts "          -fhead     out                            \\"
    puts "          -mode0     residue                        \\"
    puts "          -mode1     residue                        \\"
    puts "          -dr        0.2                            \\"
    puts "          -identical false                          \\"
    puts "          -normalize false                          \\"
    puts "          -separate_self false                      \\"
    puts "          -use_conditional  true                    \\"
    puts "          -use_symmetric    true                    \\"
    puts "          -cond_type   Z                            \\"
    puts "          -cond_range  15.0 25.0                    \\"
    puts "          -bond_range  3.5                          \\"
    puts "          -sel0      name CAL                       \\"
    puts "          -sel1      water \\"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_rdoptinfo
#! @brief         define RDF option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_rdoptinfo {} {
  
  global rdopt

  set rdopt(fhead)         "out"
  set rdopt(mode0)          residue
  set rdopt(mode1)          residue
  set rdopt(dr)             0.1
  set rdopt(identical)      false
  set rdopt(normalize)      false
  set rdopt(separate_self)  false
  set rdopt(use_conditional)  false
  set rdopt(cond_symmetric)   false
  set rdopt(cond_type)        false
  set rdopt(cond_range)       "0.0  20.0"
  set rdopt(bond_range)       3.5
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_rdoptinfo
#! @brief         read RDF option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_rdoptinfo {arglist} {

  global rdopt

  set rdopt(fhead)          [parse_arguments $arglist \
      "-fhead"             "value" $rdopt(fhead)]
  set rdopt(mode0)           [parse_arguments $arglist \
      "-mode0"             "value" $rdopt(mode0)]
  set rdopt(mode1)           [parse_arguments $arglist \
      "-mode1"             "value" $rdopt(mode1)]
  set rdopt(dr)              [parse_arguments $arglist \
      "-dr"                "value" $rdopt(dr)]
  set rdopt(identical)       [parse_arguments $arglist \
      "-identical"         "value" $rdopt(identical)]
  set rdopt(normalize)       [parse_arguments $arglist \
      "-normalize"         "value" $rdopt(normalize)]
  set rdopt(separate_self)   [parse_arguments $arglist \
      "-separate_self"     "value" $rdopt(separate_self)]
  set rdopt(use_conditional) [parse_arguments $arglist \
      "-use_conditional"   "value" $rdopt(use_conditional)]
  set rdopt(cond_symmetric)  [parse_arguments $arglist \
      "-cond_symmetric"    "value" $rdopt(cond_symmetric)]
  set rdopt(cond_type)       [parse_arguments $arglist \
      "-cond_type"         "value" $rdopt(cond_type)]
  set rdopt(cond_range)       [parse_arguments $arglist \
      "-cond_range"        "value" $rdopt(cond_range)]
  set rdopt(bond_range)       [parse_arguments $arglist \
      "-bond_range"        "value" $rdopt(bond_range)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_rdoptinfo
#! @brief         show RDF option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_rdoptinfo {} {

  global rdopt

  puts "<< option info >>"
  puts "fhead          = $rdopt(fhead)"
  puts "mode0          = $rdopt(mode0)"
  puts "mode1          = $rdopt(mode1)"
  puts "dr             = $rdopt(dr)"
  puts "identical      = $rdopt(identical)"
  puts "normalize      = $rdopt(normalize)"
  puts "separate_self  = $rdopt(separate_self)"
  puts "use_conditional = $rdopt(use_conditional)"
  if {$rdopt(cond_symmetric)} {
    puts "cond_symmetric  = $rdopt(cond_symmetric)"
    puts "cond_type       = $rdopt(cond_type)"
    puts "cond_range      = $rdopt(cond_range)"
    puts "bond_range      = $rdopt(bond_range)"
    puts ""
  }
 
}

proc rd_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global rdopt

  global sel 

  set anatra_path $::env(ANATRA_PATH);list
  set rdfort     "${anatra_path}/f90/bin/px_analysis.x";list

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

  set nsel $seltxt(nsel)
  set rand [expr int((100000*rand()))]
  set frdinp   [format "rd%06d.inp"     $rand]
  set frdout   [format "rd%06d.out"     $rand]
  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fdcdtmp($isel)  [format "rd%06d_%i.dcd"     $rand $isel]
    set fmolinfo($isel) [format "rd%06d_%i.molinfo" $rand $isel]
  }

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set rnam [$sel($isel) get resname]
    set res  [$sel($isel) get resid]
    set mass [$sel($isel) get mass]
    set anam [$sel($isel) get name]
    set natm [llength $res]
    set nf   [molinfo $mol get numframes]
    set nres [llength [lsort -unique [$sel($isel) get residue]]]

    set f [open $fmolinfo($isel) "w"]
    for {set iatm 0} {$iatm < $natm} {incr iatm} {
      puts $f [format "%10d  %6s  %6s  %15.7f" \
              [lindex $res  $iatm] \
              [lindex $rnam $iatm] \
              [lindex $anam $iatm] \
              [lindex $mass $iatm]]
    } 
    close $f

    animate write dcd $fdcdtmp($isel) \
      beg 0 end -1 waitfor all sel $sel($isel) $mol 
  }

  set f [open $frdinp "w"]
  puts $f " &input_param"
  puts $f "   fdcd = \"$fdcdtmp(0)\" \"$fdcdtmp(1)\""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead  = \"$rdopt(fhead)\""
  puts $f " /"
  
  puts $f " &trajopt_param"
  puts $f "   molinfo = \"$fmolinfo(0)\" \"$fmolinfo(1)\""
  puts $f " /"
  
  puts $f " &option_param"
  puts $f "   mode            = \"$rdopt(mode0)\" \"$rdopt(mode1)\""
  puts $f "   dr              = $rdopt(dr)"
  puts $f "   identical       = .$rdopt(identical)."
  puts $f "   normalize       = .$rdopt(normalize)."
  puts $f "   separate_self   = .$rdopt(separate_self)."
  puts $f "   use_conditional = .$rdopt(use_conditional)."
  puts $f "   cond_symmetric  = .$rdopt(cond_symmetric)."
  puts $f "   cond_type       = \"$rdopt(cond_type)\""
  puts $f "   cond_range      = $rdopt(cond_range)"
  puts $f "   bond_range      = $rdopt(bond_range)"
  puts $f " /"
  close $f

  puts "RDF is calculated with ANATRA fortran program:"
  puts "$rdfort ..."
  puts "=== INPUT ==="
  set content [exec cat $frdinp]
  puts $content
  puts "============="
  exec $rdfort $frdinp >& $frdout
  puts ""
  puts "=== OUTPUT ==="
  set content [exec cat $frdout]
  puts $content
  puts "=============="
  puts ">> Finished"
  exec rm -f $frdinp $frdout $fdcdtmp(0) $fdcdtmp(1) \
             $fmolinfo(0) $fmolinfo(1)

  exit
}

