proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                  PBC Distance Analysis"
  puts ""
  puts "============================================================"
}

proc show_pd_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra distance                                                              \\"
    puts "  -stype           <structure file type>                                     \\"
    puts "  -sfile           <structure file name>                                     \\"
    puts "  -tintype         <input trajectory file type>                              \\"
    puts "  -tin             <input trajectory file name>                              \\"
    puts "  -fhead           <header of output file name>                              \\"
    puts "  -flog            <log file name (optional)>                                \\"
    puts "  -pbc             <treat pbc or not (true or false)>                        \\"
    puts "                   (default: false)                                          \\"
    puts "  -mode0           <analysis mode of sel0 (residue or whole or atom>         \\"
    puts "  -mode1           <analysis mode of sel1 (residue or whole or atom>         \\"
    puts "  -distance_type   <standard or minimum or intra>                            \\"
    puts "                   (default: standard)>                                      \\"
    puts "  -mindist_type0   <mindist type for species 0 (site or com)>                \\" 
    puts "                   (default: site)>                                          \\"
    puts "  -mindist_type1   <mindist type for species 1 (site or com)                 \\"
    puts "                   (default: site)>                                          \\"
    puts "  -dt              <time step>                                               \\"
    puts "  -t_sta           <time at which analysis start (default: 0)>               \\"
    puts "  -t_end           <time at which analysis stop (default: 0 (till the end))> \\"
    puts "  -sel0            <VMD selection> (X=0,1,2...)                              \\"
    puts "  -sel1            <VMD selection> (X=0,1,2...)                              \\"
    puts "  -prep_only       <where analysis is performed or not (true or false)       \\"
    puts "                   (default: false)"
    puts ""
    puts "Remark : "
    puts "  o if distance_type = standard => standard distance is calculated"
    puts "  o if distance_type = minimum  => minimum distance between pairs is calculated"
    puts ""
    puts "  o mindist_typeX specifies the minimum distance type"
    puts "    if mindist_typeX = site => atomic sites are analyzed"
    puts "    if mindist_typeX = com  => CoMs are analyzed"
    puts ""  
    puts "Usage:"
    puts "anatra distance                                    \\"
    puts "  -stype          parm7                            \\"
    puts "  -sfile          str.prmtop                       \\"
    puts "  -tintype        dcd                              \\"
    puts "  -tin            inp.dcd                          \\"
    puts "  -fhead          out                              \\"
    puts "  -pbc            true                             \\"
    puts "  -distance_type  standard                         \\"
    puts "  -mode0          residue                          \\"
    puts "  -mode1          residue                          \\"
    puts "  -dt             0.1                              \\"
    puts "  -sel0           name C32  H2X H2Y and segid MEMB \\"
    puts "  -sel1           water                            \\"
    puts "  -prep_only      false"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_pdoptinfo
#! @brief         define PD option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_pdoptinfo {} {
  
  global pdopt

  set pdopt(fhead)          "out"
  set pdopt(flog)           ""
  # Remark: dsfile is deprecated option. 
  #         Instead of dsfile, please use fhead
  #set pdopt(dsfile)         ""
  set pdopt(dt)             0.1
  set pdopt(t_sta)          0.0
  set pdopt(t_end)          0.0
  set pdopt(pbc)            false
  set pdopt(distance_type)  "standard" 
  set pdopt(mindist_type0)  "site" 
  set pdopt(mindist_type1)  "site" 
  set pdopt(mode0)          "residue"
  set pdopt(mode1)          "residue"
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_pdoptinfo
#! @brief         read PD option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_pdoptinfo {arglist} {

  global pdopt

  set pdopt(fhead)         [parse_arguments $arglist \
      "-fhead"         "value" $pdopt(fhead)]
  set pdopt(flog)          [parse_arguments $arglist \
      "-flog"          "value" $pdopt(flog)]
  #set pdopt(dsfile)        [parse_arguments $arglist \
  #    "-dsfile"        "value" $pdopt(dsfile)]
  set pdopt(dt)            [parse_arguments $arglist \
      "-dt"            "value" $pdopt(dt)]
  set pdopt(t_sta)         [parse_arguments $arglist \
      "-t_sta"         "value" $pdopt(t_sta)]
  set pdopt(t_end)         [parse_arguments $arglist \
      "-t_end"         "value" $pdopt(t_end)]
  set pdopt(dt)            [parse_arguments $arglist \
      "-dt"            "value" $pdopt(dt)]
  set pdopt(pbc)           [parse_arguments $arglist \
      "-pbc"           "value" $pdopt(pbc)]
  set pdopt(distance_type) [parse_arguments $arglist \
      "-distance_type" "value" $pdopt(distance_type)]
  set pdopt(mindist_type0) [parse_arguments $arglist \
      "-mindist_type0" "value" $pdopt(mindist_type0)]
  set pdopt(mindist_type1) [parse_arguments $arglist \
      "-mindist_type1" "value" $pdopt(mindist_type1)]
  set pdopt(mode0)         [parse_arguments $arglist \
      "-mode0"         "value" $pdopt(mode0)]
  set pdopt(mode1)         [parse_arguments $arglist \
      "-mode1"         "value" $pdopt(mode1)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_pdoptinfo
#! @brief         show PD option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_pdoptinfo {} {

  global pdopt

  puts "<< option info >>"
  puts "fhead          = $pdopt(fhead)"
  #if {$pdopt(dsfile) != ""} {
  #  puts "dsfile         = $pdopt(dsfile)"
  #}
  puts "dt             = $pdopt(dt)"
  puts "t_sta          = $pdopt(t_sta)"
  puts "t_end          = $pdopt(t_end)"
  puts "pbc            = $pdopt(pbc)"
  puts "distance_type  = $pdopt(distance_type)"
  puts "mindist_type0  = $pdopt(mindist_type0)"
  puts "mindist_type1  = $pdopt(mindist_type1)"
  puts "mode0          = $pdopt(mode0)"
  puts "mode1          = $pdopt(mode1)"
  puts ""

}

proc pd_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global pdopt
  global common

  global sel 

  set anatra_path $::env(ANATRA_PATH);list
  set pdfort     "${anatra_path}/f90/bin/pd_analysis.x";list

  # read trajectory
  #
  puts ""
  puts "--------------------"
  puts " Read trajectory"
  puts "--------------------"
  puts ""

  set mol 0;
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
  set fpdinp   [format "pd%06d.inp"     $rand]

  if {$pdopt(flog) == ""} {
    set fpdout   [format "pd%06d.out"     $rand]
  } else {
    set fpdout   $pdopt(flog)
  }


  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fdcdtmp($isel)  [format "pd%06d_%i.dcd"     $rand $isel]
    set fmolinfo($isel) [format "pd%06d_%i.molinfo" $rand $isel]
  }

  if {$common(prep_only)} {
    for {set isel 0} {$isel < $nsel} {incr isel} {
      set fmolinfo($isel) [format "pd_%i.molinfo" $isel]
    }
    set fpdinp   [format "pd.inp"]
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

  #  animate write dcd $fdcdtmp($isel) \
  #    beg 0 end -1 waitfor all sel $sel($isel) $mol 
  }
  
  set ntraj [llength $traj(tin)] 

  set f [open $fpdinp "w"]
  puts $f " &input_param"
  #puts $f "   fdcd = \"$fdcdtmp(0)\" \"$fdcdtmp(1)\""
  puts $f  "  ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead = \"$pdopt(fhead)\""
  #if {$pdopt(dsfile) != ""} {
  #  puts $f "   fds   = \"$pdopt(dsfile)\""
  #}
  puts $f " /"
  
  puts $f " &trajopt_param"
  puts $f "   dt        = $pdopt(dt)"
  puts $f "   molinfo   = \"$fmolinfo(0)\" \"$fmolinfo(1)\""
  puts $f " /"
  
  puts $f " &option_param"
  puts $f "   pbc            = .$pdopt(pbc)."
  puts $f "   mode           = \"$pdopt(mode0)\" \"$pdopt(mode1)\""
  puts $f "   distance_type  = \"$pdopt(distance_type)\""
  puts $f "   mindist_type   = \"$pdopt(mindist_type0)\" \"$pdopt(mindist_type1)\""
  puts $f "   t_sta          = $pdopt(t_sta)"
  puts $f "   t_end          = $pdopt(t_end)"
  puts $f " /"
  close $f

  if {!$common(prep_only)} {
    puts "PBC Distance is calculated with ANATRA fortran program:"
    puts "$pdfort ..."
    puts "=== INPUT ==="
    set content [exec cat $fpdinp]
    puts $content
    puts "============="
    exec $pdfort $fpdinp >& $fpdout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fpdout]
    puts $content
    #exec rm -f $fpdinp $fpdout $fdcdtmp(0) $fdcdtmp(1) $fmolinfo(0) $fmolinfo(1)
    if {$pdopt(flog) == ""} {
      exec rm -f $fpdinp $fpdout $fmolinfo(0) $fmolinfo(1)
    } else {
      exec rm -f $fpdinp $fmolinfo(0) $fmolinfo(1)
    }
  }

  puts "=============="
  puts ">> Finished"

  exit
}

