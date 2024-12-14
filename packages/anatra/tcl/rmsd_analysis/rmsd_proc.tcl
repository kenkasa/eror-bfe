proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                     RMSD Analysis"
  puts ""
  puts "============================================================"
}

proc show_rmsd_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra rmsd                                               \\"
    puts "  -stype        <structure file type>                     \\"
    puts "  -sfile        <structure file name>                     \\"
    puts "  -tintype      <input trajectory file type>              \\"
    puts "  -tin          <input trajectory file name>              \\"
    puts "  -fout         <output file name>                        \\"
    puts "  -selX         <X-th VMD selection> (X=0,1,2...)         \\"
    puts "  -fit          <fit is performed or not (true or false)> \\"
    puts "                (default: false)                          \\"
    puts "  -refpdb       <reference pdb file name>                 \\"
    puts "  -fitselid     <selection id for fitting>                \\"
    puts "  -refselid     <selection id for reference>              \\"
    puts "  -rmsdselid    <selection id for rmsd>                   \\"
    puts "  -rmsdrefselid <selection id for rmsd reference>"
    puts ""
    puts "Usage:"
    puts "anatra rmsd                                   \\"
    puts "  -stype        parm7                         \\"
    puts "  -sfile        str.prmtop                    \\"
    puts "  -tintype      dcd                           \\"
    puts "  -tin          inp.dcd                       \\"
    puts "  -fout         out.rmsd                      \\"
    puts "  -sel0         resid 1 to 275 and name CA    \\"
    puts "  -sel1         resid 1 to 275 and name CA    \\"
    puts "  -sel2         resid 1 to 275 and name CA    \\"
    puts "  -sel3         resid 1 to 275 and name CA    \\"
    puts "  -fit          true                          \\"
    puts "  -refpdb       ref.pdb                       \\"
    puts "  -fitselid     0                             \\"
    puts "  -refselid     1                             \\"
    puts "  -rmsdselid    2                             \\"
    puts "  -rmsdrefselid 3                             \\"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_rmsdoptinfo
#! @brief         define RMSD option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_rmsdoptinfo {} {
  
  global rmsdopt

  set rmsdopt(fout)          "out.rmsd"
  set rmsdopt(fit)           false
  set rmsdopt(wrap)          false
  set rmsdopt(centering)     false

  set rmsdopt(wrapcenter)    origin
  set rmsdopt(wrapcomp)      fragment
  set rmsdopt(fitselid)      0 
  set rmsdopt(refselid)      0 
  set rmsdopt(rmsdselid)     0 
  set rmsdopt(rmsdrefselid)  0 
  set rmsdopt(refpdb)        ""
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_rmsdoptinfo
#! @brief         read trajectory option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_rmsdoptinfo {arglist} {

  global rmsdopt

  set rmsdopt(fout)          [parse_arguments $arglist \
      "-fout"         "value" $rmsdopt(fout)]
  set rmsdopt(fit)           [parse_arguments $arglist \
      "-fit"          "value" $rmsdopt(fit)]
  set rmsdopt(fitselid)      [parse_arguments $arglist \
      "-fitselid"     "value" $rmsdopt(fitselid)]
  set rmsdopt(refselid)      [parse_arguments $arglist \
      "-refselid"     "value" $rmsdopt(refselid)]
  set rmsdopt(rmsdselid)     [parse_arguments $arglist \
      "-rmsdselid"    "value" $rmsdopt(rmsdselid)]
  set rmsdopt(rmsdrefselid)  [parse_arguments $arglist \
      "-rmsdrefselid" "value" $rmsdopt(rmsdrefselid)]
  set rmsdopt(refpdb)        [parse_arguments $arglist \
      "-refpdb"       "value" $rmsdopt(refpdb)]

  # Combination error check
  #

  if {$rmsdopt(refpdb) == ""} {
    puts "ERROR: refpdb should be specified."
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      fitting 
#! @brief         perform fitting the structure to reference structure 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc fitting {molid seltxt_fit seltxt_out refid seltxt_ref} {
  set refsel [atomselect $refid $seltxt_ref]
  set fitsel [atomselect $molid $seltxt_fit]
  set outsel [atomselect $molid $seltxt_out]

  set nf [molinfo $molid get numframes]
   
  for {set i 0} {$i < $nf} {incr i} {
    if {[expr $i % 100] == 0 || [expr $i + 1] == $nf} {
      puts [format "%10d / %10d" $i [expr $nf - 1]] 
    } 
    $fitsel frame $i
    $outsel frame $i
    $outsel move [measure fit $fitsel $refsel weight mass]
  } 

}	

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_rmsdoptinfo
#! @brief         show trajectory option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_rmsdoptinfo {} {

  global rmsdopt

  puts "<< option info >>"
  puts "fout         = $rmsdopt(fout)"
  puts "fit          = $rmsdopt(fit)"
  puts "fitselid     = $rmsdopt(fitselid)"
  puts "refselid     = $rmsdopt(refselid)"
  puts "rmsdselid    = $rmsdopt(rmsdselid)"
  puts "rmsdrefselid = $rmsdopt(rmsdrefselid)"
  puts "refpdb       = $rmsdopt(refpdb)"
  puts ""

}

proc rmsd_analyze {} {
  # in
  global str
  global traj
  global seltxt
  global rmsdopt

  global sel 

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

  # setup reference 
  #
  puts ""
  puts "--------------------"
  puts " Setup reference"
  puts "--------------------"
  puts ""
  set ref        [mol load pdb "$rmsdopt(refpdb)"]
  set refsel     [atomselect $ref "$seltxt($rmsdopt(refselid))"]
  set rmsdrefsel [atomselect $ref "$seltxt($rmsdopt(rmsdrefselid))"]
  set comref     [measure center $refsel weight mass]

  set fitsel     [atomselect $mol "$seltxt($rmsdopt(fitselid))"]

  # Convert 
  #
  puts ""
  puts "--------------------"
  puts " Start analysis"
  puts "--------------------"
  puts ""

  set nf [molinfo $mol get numframes] 

  if {$rmsdopt(fit)} {
    puts "o Start fitting"
    fitting $mol \
            "$seltxt($rmsdopt(fitselid))" \
            "all"                         \
	          $ref                          \
	          "$seltxt($rmsdopt(refselid))"
    puts "> Finish fitting" 
     
  }

  for {set i 0} {$i < $nf} {incr i} {
    $sel($rmsdopt(rmsdselid)) frame $i
    set rmsval($i) [measure rmsd                 \
                    $sel($rmsdopt(rmsdselid))    \
                    $rmsdrefsel                  \
                    weight mass]
  }

  set fid [open "$rmsdopt(fout)" w]
  for {set i 0} {$i < $nf} {incr i} {
    puts $fid [format "%10d  %15.7f" $i $rmsval($i)]
  }
  close $fid

  puts ""
  puts ">> Finish all Analysis"
  puts ""

}
