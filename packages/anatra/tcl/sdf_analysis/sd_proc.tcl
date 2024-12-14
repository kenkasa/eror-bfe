proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                       SDF Analysis"
  puts ""
  puts "============================================================"
}

proc show_tr_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra sd                                                                                      \\"
    puts "  -stype             <structure file type>                                                     \\"
    puts "  -sfile             <structure file name>                                                     \\"
    puts "  -tintype           <input trajectory file type>                                              \\"
    puts "  -tin               <input trajectory file name>                                              \\"
    puts "  -flist             <input list file that contains list of cv files>                          \\"
    puts "                     (necessary if use_restriction is true)                                    \\"
    puts "  -flist_weight      <input list file that contains list of weight files>                      \\"
    puts "  -fhead             <header of output file name>                                              \\"
    puts "  -sel0              <VMD selection> (X=0,1,2...)                                              \\"
    puts "  -mode              <analysis mode (residue or whole or atom)                                 \\"
    puts "                     (default: residue)>                                                       \\"
    puts "  -ng3               <number of grids for x, y, z axes>                                        \\"
    puts "  -del               <grid spacing for x, y, z axes>                                           \\"
    puts "  -origin            <origin of 3d-grids>                                                      \\"
    puts "  -use_spline        <whether spline is performed or not (true or false)>                      \\"
    puts "                     (default: false)>                                                         \\"
    puts "  -spline_resolution <spline resolution (integer)>                                             \\"
    puts "                     (default: 4)                                                              \\"
    puts "  -use_restriction   <whether restricted sampling is used or not>                              \\"
    puts "                     (true or false) (default: false)                                          \\"
    puts "  -normalize_reacsnap <whether normalize is done with reactive config>                         \\"
    puts "                     (true or false) (default: false)                                          \\"
    puts "  -ndim              <dimensions of reaction coords (neccesary if use_restriction is true)>    \\"
    puts "  -react_range       <range of sampled reaction coords (neccesary if use_restriction is true)> \\"
    puts "  -prep_only         <where analysis is performed or not (true or false)>                      \\"
    puts "                     (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra sd                                    \\"
    puts "  -stype             parm7                   \\"
    puts "  -sfile             str.prmtop              \\"
    puts "  -tintype           dcd                     \\"
    puts "  -tin               inp.dcd                 \\"
    puts "  -fhead             out                     \\"
    puts "  -sel0              not water               \\"
    puts "  -mode              residue                 \\"
    puts "  -ng3               50 50 50                \\"
    puts "  -del               0.4 0.4 0.4             \\"
    puts "  -origin            0.0 0.0 0.0             \\"
    puts "  -use_spline        false                   \\"
    puts "  -spline_resolution 4                       \\"
    puts "  -prep_only         false"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_sdoptinfo
#! @brief         define SDF option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_sdoptinfo {} {
  
  global sdopt

  set sdopt(fhead)             "run" 
  set sdopt(mode)              "residue" 
  set sdopt(ng3)               "50 50 50" 
  set sdopt(del)               "0.5 0.5 0.5" 
  set sdopt(origin)            "0.0 0.0 0.0"
  set sdopt(use_pbcwrap)       false
  set sdopt(centertype)        "ZERO"
  set sdopt(use_restriction)   false
  set sdopt(ndim)              1
  set sdopt(react_range)       "0.0 0.0"
  set sdopt(fcv)               ""
  set sdopt(flist)             ""
  set sdopt(flist_weight)      ""
  set sdopt(use_spline)        false 
  set sdopt(spline_resolution) 4 
  set sdopt(is_whole_snap)     true 
  set sdopt(normalize_reacsnap) false
  set sdopt(nstep_whole)       0
  set sdopt(count_threshold)   1.0e-10


}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_sdoptinfo
#! @brief         read SDF option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_sdoptinfo {arglist} {

  global sdopt

  set sdopt(fhead)                [parse_arguments $arglist \
      "-fhead"                    "value" $sdopt(fhead)]
  set sdopt(mode)                 [parse_arguments $arglist \
      "-mode"                     "value" $sdopt(mode)]
  set sdopt(ng3)                  [parse_arguments $arglist \
      "-ng3"                      "value" $sdopt(ng3)]
  set sdopt(del)                  [parse_arguments $arglist \
      "-del"                      "value" $sdopt(del)]
  set sdopt(origin)               [parse_arguments $arglist \
      "-origin"                   "value" $sdopt(origin)]

  # Hidden options
  set sdopt(use_pbcwrap)          [parse_arguments $arglist \
      "-use_pbcwrap"              "value" $sdopt(use_pbcwrap)]
  set sdopt(centertype)           [parse_arguments $arglist \
      "-centertype"               "value" $sdopt(centertype)]

  set sdopt(use_restriction)      [parse_arguments $arglist \
      "-use_restriction"          "value" $sdopt(use_restriction)]
  set sdopt(ndim)                 [parse_arguments $arglist \
      "-ndim"                     "value" $sdopt(ndim)]
  set sdopt(react_range)          [parse_arguments $arglist \
      "-react_range"              "value" $sdopt(react_range)]
  set sdopt(fcv)                 [parse_arguments $arglist \
      "-fcv"                      "value" $sdopt(fcv)]
  set sdopt(flist)               [parse_arguments $arglist \
      "-flist"                    "value" $sdopt(flist)]
  set sdopt(flist_weight)        [parse_arguments $arglist \
      "-flist_weight"             "value" $sdopt(flist_weight)]
  set sdopt(use_spline)           [parse_arguments $arglist \
      "-use_spline"               "value" $sdopt(use_spline)]
  set sdopt(spline_resolution)    [parse_arguments $arglist \
      "-spline_resolution"        "value" $sdopt(spline_resolution)]
  set sdopt(is_whole_snap)        [parse_arguments $arglist \
      "-is_whole_snap"            "value" $sdopt(is_whole_snap)]
  set sdopt(normalize_reacsnap)   [parse_arguments $arglist \
      "-normalize_reacsnap"       "value" $sdopt(normalize_reacsnap)]
  set sdopt(nstep_whole)          [parse_arguments $arglist \
      "-nstep_whole"              "value" $sdopt(nstep_whole)]
  set sdopt(count_threshold)      [parse_arguments $arglist \
      "-count_threshold"          "value" $sdopt(count_threshold)]

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_sdoptinfo
#! @brief         show CoM option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_sdoptinfo {} {

  global sdopt

  puts "<< option info >>"
  puts "fhead             = $sdopt(fhead)"
  puts "mode              = $sdopt(mode)"
  puts "ng3               = $sdopt(ng3)"
  puts "del               = $sdopt(del)"
  puts "origin            = $sdopt(origin)"

  if {$sdopt(use_pbcwrap)} {
    puts "use_pbcwrap       = $sdopt(use_pbcwrap)"
    puts "centertype        = $sdopt(centertype)"
  }

  puts "use_restriction   = $sdopt(use_restriction)"
  puts "ndim              = $sdopt(ndim)"
  puts "react_range       = $sdopt(react_range)"
  puts "fcv               = $sdopt(fcv)"
  puts "flist             = $sdopt(flist)"
  puts "flist_weight      = $sdopt(flist_weight)"
  puts "use_spline        = $sdopt(use_spline)"
  puts "spline_resolution = $sdopt(spline_resolution)"
  puts "is_whole_snap     = $sdopt(is_whole_snap)"
  puts "normalize_reacsnap= $sdopt(normalize_reacsnap)"
  puts "nstep_whole       = $sdopt(nstep_whole)"
  puts "count_threshold   = $sdopt(count_threshold)"
  puts ""

}

proc sd_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global sdopt
  global common

  global sel 


  set anatra_path $::env(ANATRA_PATH);list
  set sdfort     "${anatra_path}/f90/bin/sd_analysis.x";list

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

  #set nf [molinfo $mol get numframes] 
  puts ""

  puts ">> Start CoM calculation"
  puts ""

  set rand [expr int((100000*rand()))]
  set fmolinfo [format "sd%06d.molinfo" $rand]
  set fsdinp   [format "sd%06d.inp"     $rand]
  set fsdout   [format "sd%06d.out"     $rand]
  set fdcdtmp  [format "sd%06d.dcd"     $rand]


  if {$common(prep_only)} {
    set fmolinfo [format "sd.molinfo"]
    set fsdinp   [format "sd.inp"]
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

  set f [open $fsdinp "w"]
  puts $f " &input_param"
  #puts $f "   fdcd0 = \"$fdcdtmp\""
  puts $f  "  ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  
  puts $f "   fcv          = \"$sdopt(fcv)\""
  puts $f "   flist        = \"$sdopt(flist)\""
  puts $f "   flist_weight = \"$sdopt(flist_weight)\""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead     = \"$sdopt(fhead)\""
  puts $f " /"

  puts $f " &trajopt_param"
  puts $f "   molinfo = \"$fmolinfo\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   mode              = \"$sdopt(mode)\""
  puts $f "   ng3               = $sdopt(ng3)"
  puts $f "   del               = $sdopt(del)"
  puts $f "   origin            = $sdopt(origin)"
  puts $f "   use_pbcwrap       = .$sdopt(use_pbcwrap)."
  puts $f "   centertype        = \"$sdopt(centertype)\""
  puts $f "   use_restriction   = .$sdopt(use_restriction)."
  puts $f "   ndim              = $sdopt(ndim)"
  puts $f "   react_range       = $sdopt(react_range)"
  puts $f "   use_spline        = .$sdopt(use_spline)."
  puts $f "   spline_resolution = $sdopt(spline_resolution)"
  puts $f "   is_whole_snap     = .$sdopt(is_whole_snap)."
  puts $f "   normalize_reacsnap= .$sdopt(normalize_reacsnap)."
  puts $f "   nstep_whole       = $sdopt(nstep_whole)"
  puts $f "   count_threshold   = $sdopt(count_threshold)"
  puts $f " /"

  close $f

  
  #animate write dcd $fdcdtmp beg 0 end -1 waitfor all sel $sel(0) $mol

  if {!$common(prep_only)} {
    puts "SDF is calculated with ANATRA fortran program:"
    puts "$sdfort ..."
    puts "=== INPUT ==="
    set content [exec cat $fsdinp]
    puts $content
    puts "============="
    exec $sdfort $fsdinp >& $fsdout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fsdout]
    puts $content
    
    exec rm $fmolinfo $fsdinp $fsdout 
  } 
  puts "=============="
  puts ">> Finished"

  exit
}
