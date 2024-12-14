proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                       RDF Analysis"
  puts ""
  puts "============================================================"
}

proc show_rd_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra rdf                                            \\"
    puts "  -stype         <structure file type>                \\"
    puts "  -sfile         <structure file name>                \\"
    puts "  -tintype       <input trajectory file type>         \\"
    puts "  -tin           <input trajectory file name>         \\"
    puts "  -fhead         <header of output file name>         \\"
    puts "  -flog          <log file name (optional)>           \\"
    puts "  -sel0          <VMD selection> (X=0,1,2...)         \\"
    puts "  -sel1          <VMD selection> (X=0,1,2...)         \\"
    puts "  -mode0         <analysis mode of sel0>              \\"
    puts "                 (residue or whole or atom)           \\"
    puts "                 (default: residue)>                  \\"
    puts "  -mode1         <analysis mode of sel1>              \\"
    puts "                 (residue or whole or atom)           \\"
    puts "                 (default: residue)                   \\"
    puts "  -dr            <delta r (angstrom)>                 \\"
    puts "  -identical     <true or false>                      \\"
    puts "                 (default: false)                     \\"
    puts "  -normalize     <true or false>                      \\"
    puts "                 (default: false)                     \\"
    puts "  -separate_self <true or false>                      \\"
    puts "                 (default: false)                     \\"
    puts "  -prep_only     <where analysis is performed or not> \\"
    puts "                 (true or false)                      \\"
    puts "                 (default: false)>"
    puts ""  
    puts "Usage:"
    puts "anatra rdf                                        \\"
    puts "  -stype         parm7                            \\"
    puts "  -sfile         str.prmtop                       \\"
    puts "  -tintype       dcd                              \\"
    puts "  -tin           inp.dcd                          \\"
    puts "  -fhead         out                              \\"
    puts "  -sel0          name C32  H2X H2Y and segid MEMB \\"
    puts "  -sel1          water                            \\"
    puts "  -mode0         residue                          \\"
    puts "  -mode1         residue                          \\"
    puts "  -dr            0.4                              \\"
    puts "  -identical     false                            \\"
    puts "  -normalize     false                            \\"
    puts "  -separate_self false                            \\"
    puts "  -prep_only     false                            \\"
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

  set rdopt(rdfile)         ""
  set rdopt(rdbinfile)      ""
  set rdopt(fhead)          "out"
  set rdopt(flog)           ""
  set rdopt(mode0)          residue
  set rdopt(mode1)          residue
  set rdopt(dr)             0.1
  set rdopt(identical)      false
  set rdopt(normalize)      true 
  set rdopt(separate_self)  false
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

  set rdopt(rdfile)    [parse_arguments $arglist \
      "-rdfile"     "value" $rdopt(rdfile)]
  set rdopt(rdbinfile) [parse_arguments $arglist \
      "-rdbinfile"  "value" $rdopt(rdbinfile)]
  set rdopt(fhead)     [parse_arguments $arglist \
      "-fhead"      "value" $rdopt(fhead)]
  set rdopt(flog)      [parse_arguments $arglist \
      "-flog"       "value" $rdopt(flog)]
  set rdopt(mode0)     [parse_arguments $arglist \
      "-mode0"      "value" $rdopt(mode0)]
  set rdopt(mode1)     [parse_arguments $arglist \
      "-mode1"      "value" $rdopt(mode1)]
  set rdopt(dr)       [parse_arguments $arglist \
      "-dr"         "value" $rdopt(dr)]
  set rdopt(identical) [parse_arguments $arglist \
      "-identical"  "value" $rdopt(identical)]
  set rdopt(normalize) [parse_arguments $arglist \
      "-normalize"  "value" $rdopt(normalize)]
  set rdopt(separate_self) [parse_arguments $arglist \
      "-separate_self"  "value" $rdopt(separate_self)]
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
  puts ""

}

proc rd_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global rdopt
  global common

  global sel 

  set anatra_path $::env(ANATRA_PATH);list
  set rdfort     "${anatra_path}/f90/bin/rd_analysis.x";list

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
  set frdinp   [format "rd%06d.inp"     $rand]

  if {$rdopt(flog) == ""} {
    set frdout   [format "rd%06d.out"     $rand]
  } else {
    set frdout   $rdopt(flog) 
  }

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fdcdtmp($isel)  [format "rd%06d_%i.dcd"     $rand $isel]
    set fmolinfo($isel) [format "rd%06d_%i.molinfo" $rand $isel]
  }

  if {$common(prep_only)} {
    for {set isel 0} {$isel < $nsel} {incr isel} {
      set fmolinfo($isel) [format "rd_%i.molinfo" $isel]
    }
    set frdinp   [format "rd.inp"]
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

    #animate write dcd $fdcdtmp($isel) \
    #  beg 0 end -1 waitfor all sel $sel($isel) $mol 
  }

  set ntraj [llength $traj(tin)] 

  set f [open $frdinp "w"]
  puts $f " &input_param"
  #puts $f "   fdcd = \"$fdcdtmp(0)\" \"$fdcdtmp(1)\""
  puts $f "   ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  puts $f ""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   frd    = \"$rdopt(rdfile)\""
  puts $f "   frdbin = \"$rdopt(rdbinfile)\""
  puts $f "   fhead  = \"$rdopt(fhead)\""
  puts $f " /"
  
  puts $f " &trajopt_param"
  puts $f "   molinfo = \"$fmolinfo(0)\" \"$fmolinfo(1)\""
  puts $f " /"
  
  puts $f " &option_param"
  puts $f "   mode      = \"$rdopt(mode0)\" \"$rdopt(mode1)\""
  puts $f "   dr        = $rdopt(dr)"
  puts $f "   identical = .$rdopt(identical)."
  puts $f "   normalize = .$rdopt(normalize)."
  puts $f "   separate_self = .$rdopt(separate_self)."
  puts $f " /"
  close $f

  if {!$common(prep_only)} {
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

    if {$rdopt(flog) == ""} {
      exec rm -f $frdinp $frdout $fmolinfo(0) $fmolinfo(1)
    } else {
      exec rm -f $frdinp $fmolinfo(0) $fmolinfo(1)
    }
  }

  puts "=============="
  puts ">> Finished"

  exit
}

