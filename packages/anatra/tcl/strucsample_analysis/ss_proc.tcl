proc print_title {} {
  puts "============================================================"
  puts ""
  puts "               Structural Sampling Analysis"
  puts ""
  puts "============================================================"
}

proc show_ss_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra ss                                                         \\"
    puts "  -stype      <structure file type>                               \\"
    puts "  -sfile      <structure file name>                               \\"
    puts "  -tintype    <input trajectory file type>                        \\"
    puts "  -tin        <input trajectory file name>                        \\"
    puts "  -totype     <output trajectory file type>                       \\"
    puts "  -to         <output trajectory file name>                       \\"
    puts "  -fhead      <header of output file>                             \\"
    puts "  -nsample    <number of extracted snapshots>                     \\"
    puts "  -out_rst7   <whether rst7 files are generated or not>           \\"
    puts "              (true or false)                                     \\"
    puts "              (default: false)                                    \\"
    puts "  -use_allsnap <whether all the snapshots are used or not>        \\"
    puts "              (true or false)                                     \\"
    puts "              (default: false)                                    \\"
    puts "  -duplicate  <extract the identical snapshots repeatedly or not> \\"
    puts "              (true or false)                                     \\"
    puts "              (default: false)                                    \\"
    puts "  -iseed      <input random seed (integer)>                       \\"
    puts "              (if not specified, determined from timeclock)>      \\"
    puts "  -prep_only  <where analysis is performed or not (true or not)>  \\"
    puts "              (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra ss                          \\"
    puts "  -stype     parm7                 \\"
    puts "  -sfile     inp.prmtop            \\"
    puts "  -tintype   dcd                   \\"
    puts "  -tin       inp.dcd               \\"
    puts "  -totype    dcd                   \\"
    puts "  -to        out.dcd               \\"
    puts "  -fhead     out                   \\"
    puts "  -nsample   50                    \\"
    puts "  -duplicate false                 \\"
    puts "  -iseed     3141592               \\"
    puts "  -prep_only false"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_ssoptinfo
#! @brief         define CoM option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_ssoptinfo {} {
  
  global ssopt

  set ssopt(nsample)     100
  set ssopt(fhead)       "out"
  set ssopt(out_rst7)    "false"
  set ssopt(use_allsnap) "false"
  set ssopt(duplicate)   "false"
  set ssopt(iseed)       -1 
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_ssoptinfo
#! @brief         read CoM option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_ssoptinfo {arglist} {

  global ssopt

  set ssopt(nsample)    [parse_arguments $arglist \
      "-nsample"    "value" $ssopt(nsample)]
  set ssopt(fhead)      [parse_arguments $arglist \
      "-fhead"      "value" $ssopt(fhead)]
  set ssopt(out_rst7)   [parse_arguments $arglist \
      "-out_rst7"   "value" $ssopt(out_rst7)]
  set ssopt(use_allsnap) [parse_arguments $arglist \
      "-use_allsnap" "value" $ssopt(use_allsnap)]
  set ssopt(duplicate)  [parse_arguments $arglist \
      "-duplicate"  "value" $ssopt(duplicate)]
  set ssopt(iseed)      [parse_arguments $arglist \
      "-iseed"      "value" $ssopt(iseed)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_ssoptinfo
#! @brief         show CoM option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_ssoptinfo {} {

  global ssopt

  puts "<< option info >>"
  puts "nsample     = $ssopt(nsample)"
  puts "fhead       = $ssopt(fhead)"
  puts "out_rst7    = $ssopt(out_rst7)"
  puts "use_allsnap = $ssopt(use_allsnap)"
  puts "duplicate   = $ssopt(duplicate)"
  puts "iseed       = $ssopt(iseed)"
  puts ""

}

proc ss_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global ssopt
  global common

  global sel 


  set anatra_path $::env(ANATRA_PATH);list
  set ssfort     "${anatra_path}/f90/bin/ss_analysis.x";list

  # Convert 
  #
  puts ""
  puts "--------------------"
  puts " Start analysis"
  puts "--------------------"
  puts ""

  puts ">> Start calculation"
  puts ""

  set rand [expr int((100000*rand()))]
  set fssinp   [format "ss%06d.inp"     $rand]
  set fssout   [format "ss%06d.out"     $rand]
  set fdcdtmp  [format "ss%06d.dcd"     $rand]

  if {$common(prep_only)} {
    set fssinp   [format "ss.inp"]
  }

  set ntraj [llength $traj(tin)]

  set f [open $fssinp "w"]
  puts $f " &input_param"
  puts $f  "  ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  puts $f ""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   ftraj = \"$traj(to)\""
  puts $f "   fhead = \"$ssopt(fhead)\""
  puts $f " /"
  puts $f " &option_param"
  puts $f "   nsample     = $ssopt(nsample)"
  puts $f "   duplicate   = .$ssopt(duplicate)."
  puts $f "   iseed       = $ssopt(iseed)"
  puts $f "   use_allsnap = .$ssopt(use_allsnap)."
  puts $f "   out_rst7    = .$ssopt(out_rst7)."
  puts $f " /"

  close $f

  if {!$common(prep_only)} {
    puts "StrucSample is calculated with ANATRA fortran program:"
    puts "$ssfort ..."
    puts "=== INPUT ==="
    set content [exec cat $fssinp]
    puts $content
    puts "============="
    exec $ssfort $fssinp >& $fssout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fssout]
    puts $content
    exec rm $fssinp $fssout 
  }

  puts "=============="
  puts ">> Finished"


  #set mol 0;
  #read_traj $mol $str(stype) $str(sfile) dcd $fdcdtmp 1
  #set nf [molinfo $mol get numframes]
  #
  #if {$traj(totype) == "rst7"} {
  #  set ext "inpcrd"
  #} else {
  #  set ext $traj(totype)
  #}

  #for {set istep 0} {$istep < $nf} {incr istep} {
  #  set j [expr $istep + 1]
  #  set fout [format "%s%06d.%s" $ssopt(fhead) $j $ext]
  #  animate write $traj(totype) $fout beg $istep end $istep waitfor all
  #}
  #exec rm $fdcdtmp
  #
  exit
}
