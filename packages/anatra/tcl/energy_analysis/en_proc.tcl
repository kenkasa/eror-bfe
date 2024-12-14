proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                     Energy Analysis"
  puts ""
  puts "============================================================"
}

proc show_en_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra en                                                          \\"
    puts "  -stype      <structure file type>                                \\"
    puts "  -sfile      <structure file name>                                \\"
    puts "  -tintype    <input trajectory file type>                         \\"
    puts "  -tin        <input trajectory file name>                         \\"
    puts "  -fhead      <header of output file name>                         \\"
    puts "  -parmformat <parameter file format (prmtop or anaparm)>          \\"
    puts "              (default: prmtop)                                    \\"
    puts "  -fanaparm   <anaparam file>                                      \\"
    puts "              (necessary if parmformat = anaparm)                  \\"
    puts "  -pbc        <treat pbc or not (true or false)>                   \\"
    puts "              (default: false)                                     \\"
    puts "  -calc_vdw   <whether vdw is calculated (true or false)>          \\"
    puts "              (default: true)                                      \\"
    puts "  -calc_elec  <whether elec is calculated (true or false)>         \\"
    puts "              (default: false)                                     \\"
    puts "  -calc_torque <whether torque is calculated (true or false)>      \\"
    puts "              (default: false)                                     \\"
    puts "  -calc_siteesp <whether electrostatic potential is calculated or not> \\"
    puts "              (default: false)                                     \\"
    puts "  -vdw        <vdw interaction type (standard or attractive)>      \\" 
    puts "  -elec       <elec interaction type (bare or pme)>                \\"
    puts "              (default: bare)>                                     \\" 
    puts "  -dt         <time step>                                          \\"
    puts "  -rljcut     <LJ cutoff distance (A)>                             \\"
    puts "  -relcut     <ELEC cutoff distance (A)>                           \\"
    puts "  -pme_alpha  <PME screening parameter (A^-1)>                     \\"
    puts "  -pme_grids  <PME grids>                                          \\"
    puts "  -pme_rigid  <PME solute rigidity (true or false)>                \\"
    puts "              (default: false)                                     \\"
    puts "  -pme_dual   <two different solute potentials are used>           \\"
    puts "              (true or false)                                      \\"
    puts "              (default: false)                                     \\"
    puts "  -pme_fprm   <additional potential parameter file for dual calc.> \\"
    puts "  -sel0       <VMD selection> (X=0,1,2...)                         \\"
    puts "  -sel1       <VMD selection> (X=0,1,2...)                         \\"
    puts "  -mode0      <analysis mode of sel0 (residue or whole or atom>    \\"
    puts "  -mode1      <analysis mode of sel1 (residue or whole or atom>    \\"
    puts "  -prep_only  <where analysis is performed or not (true or false)> \\"
    puts "              (default: false)"
    puts ""  
    puts "Usage:"
    puts "anatra en                                        \\"
    puts "  -stype        parm7                            \\"
    puts "  -sfile        str.prmtop                       \\"
    puts "  -tintype      dcd                              \\"
    puts "  -tin          inp.dcd                          \\"
    puts "  -fhead        out                              \\"
    puts "  -parmformat   prmtop                           \\"
    puts "  -fanaparm     complex.anaparm                  \\"
    puts "  -pbc          true                             \\"
    puts "  -calc_vdw     true                             \\"
    puts "  -calc_elec    false                            \\"
    puts "  -calc_torque  false                            \\"
    puts "  -calc_espsite false                            \\"
    puts "  -vdw          standard                         \\"
    puts "  -elec         bare                             \\"
    puts "  -mode0        residue                          \\"
    puts "  -mode1        residue                          \\"
    puts "  -dt           0.1                              \\"
    puts "  -rljcut       12.0                             \\"
    puts "  -relcut       1.0e10                           \\"
    puts "  -pme_alpha    0.35e0                           \\"
    puts "  -pme_grids    64 64 64                         \\"
    puts "  -pme_rigid    false                            \\"
    puts "  -pme_dual     false                            \\"
    puts "  -pme_fprm     excited.prmtop                   \\"
    puts "  -sel0         name C32  H2X H2Y and segid MEMB \\"
    puts "  -sel1         water                            \\"
    puts "  -prep_only    false"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_enoptinfo
#! @brief         define EN option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_enoptinfo {} {
  
  global enopt

  set enopt(fhead)         "out"
  set enopt(parmformat)    "prmtop"
  set enopt(fanaparm)      "complex.anapram"
  set enopt(dt)             0.1
  set enopt(pbc)            false
  set enopt(calc_vdw)       true 
  set enopt(calc_elec)      false 
  set enopt(calc_torque)    false 
  set enopt(calc_siteesp)   false 
  set enopt(vdw)            "standard"
  set enopt(elec)           "bare"
  set enopt(mode0)          "residue"
  set enopt(mode1)          "residue"
  set enopt(rljcut)         12.0 
  set enopt(relcut)         1.0e10 
  set enopt(pme_alpha)      0.35e0
  set enopt(pme_grids)      "64 64 64"
  set enopt(pme_rigid)      false
  set enopt(pme_dual)       false
  set enopt(pme_fprm)       "excited.prmtop" 
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_enoptinfo
#! @brief         read PD option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_enoptinfo {arglist} {

  global enopt

  set enopt(fhead)      [parse_arguments $arglist \
      "-fhead"       "value" $enopt(fhead)]
  set enopt(parmformat) [parse_arguments $arglist \
      "-parmformat"  "value" $enopt(parmformat)]
  set enopt(fanaparm)   [parse_arguments $arglist \
      "-fanaparm"    "value" $enopt(fanaparm)]
  set enopt(dt)         [parse_arguments $arglist \
      "-dt"          "value" $enopt(dt)]
  set enopt(pbc)        [parse_arguments $arglist \
      "-pbc"         "value" $enopt(pbc)]
  set enopt(calc_vdw)   [parse_arguments $arglist \
      "-calc_vdw"    "value" $enopt(calc_vdw)]
  set enopt(calc_elec)  [parse_arguments $arglist \
      "-calc_elec"   "value" $enopt(calc_elec)]
  set enopt(calc_torque) [parse_arguments $arglist \
      "-calc_torque" "value" $enopt(calc_torque)]
  set enopt(calc_siteesp) [parse_arguments $arglist \
      "-calc_siteesp" "value" $enopt(calc_siteesp)]
  set enopt(vdw)        [parse_arguments $arglist \
      "-vdw"         "value" $enopt(vdw)]
  set enopt(elec)       [parse_arguments $arglist \
      "-elec"        "value" $enopt(elec)]
  set enopt(mode0)      [parse_arguments $arglist \
      "-mode0"       "value" $enopt(mode0)]
  set enopt(mode1)      [parse_arguments $arglist \
      "-mode1"       "value" $enopt(mode1)]
  set enopt(rljcut)     [parse_arguments $arglist \
      "-rljcut"      "value" $enopt(rljcut)]
  set enopt(relcut)     [parse_arguments $arglist \
      "-relcut"      "value" $enopt(relcut)]
  set enopt(pme_alpha)  [parse_arguments $arglist \
      "-pme_alpha"   "value" $enopt(pme_alpha)]
  set enopt(pme_grids)  [parse_arguments $arglist \
      "-pme_grids"   "value" $enopt(pme_grids)]
  set enopt(pme_rigid)  [parse_arguments $arglist \
      "-pme_rigid"   "value" $enopt(pme_rigid)]
  set enopt(pme_dual)   [parse_arguments $arglist \
      "-pme_dual"    "value" $enopt(pme_dual)]
  set enopt(pme_fprm)   [parse_arguments $arglist \
      "-pme_fprm"    "value" $enopt(pme_fprm)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_enoptinfo
#! @brief         show PD option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_enoptinfo {} {

  global enopt

  puts "<< option info >>"
  puts "fhead          = $enopt(fhead)"
  puts "parmformat     = $enopt(parmformat)"
  puts "fanaparm       = $enopt(fanaparm)"
  puts "dt             = $enopt(dt)"
  puts "pbc            = $enopt(pbc)"
  puts "calc_vdw       = $enopt(calc_vdw)"
  puts "calc_elec      = $enopt(calc_elec)"
  puts "calc_torque    = $enopt(calc_torque)"
  puts "calc_siteesp   = $enopt(calc_siteesp)"
  puts "vdw            = $enopt(vdw)"
  puts "elec           = $enopt(elec)"
  puts "mode0          = $enopt(mode0)"
  puts "mode1          = $enopt(mode1)"
  puts "rljcut         = $enopt(rljcut)"
  puts "relcut         = $enopt(relcut)"
  puts "pme_alpha      = $enopt(pme_alpha)"
  puts "pme_grids      = $enopt(pme_grids)"
  puts "pme_rigid      = $enopt(pme_rigid)"
  puts "pme_dual       = $enopt(pme_dual)"
  puts "pme_fprm       = $enopt(pme_fprm)"
  puts ""

}

proc en_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global enopt
  global common

  global sel 

  set anatra_path $::env(ANATRA_PATH);list
  set enfort     "${anatra_path}/f90/bin/energy_analysis.x";list


  # check control parameter check
  #
  #if {$str(stype) != "parm7"} {
  #  puts "Error: only parm7 structure file is supported in this analysis."
  #  exit
  #}

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
  set feninp   [format "en%06d.inp"     $rand]
  set fenout   [format "en%06d.out"     $rand]
  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fdcdtmp($isel)  [format "en%06d_%i.dcd"     $rand $isel]
    set fmolinfo($isel) [format "en%06d_%i.molinfo" $rand $isel]
  }

  if {$common(prep_only)} {
    for {set isel 0} {$isel < $nsel} {incr isel} {
      set fmolinfo($isel) [format "en_%i.molinfo" $isel]
    }
    set feninp   [format "en.inp"]
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

  set f [open $feninp "w"]
  puts $f " &input_param"
  puts $f "   fprmtop  = \"$str(sfile)\""
  puts $f "   ftraj ="
  for {set i 0} {$i < $ntraj} {incr i} {
    set t [lindex $traj(tin) $i]
    puts -nonewline $f "    \"$t\" "
  }
  puts $f ""

  if {$enopt(pme_dual)} {
    puts $f "   fprmtop2   = \"$enopt(pme_fprm)\""
    puts $f "   fanaparm2  = \"$enopt(pme_fprm)\""
  }
  puts $f "   fanaparm = \"$enopt(fanaparm)\""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead = \"$enopt(fhead)\""
  puts $f " /"
  
  puts $f " &trajopt_param"
  puts $f "   dt        = $enopt(dt)"
  puts $f "   molinfo   = \"$fmolinfo(0)\" \"$fmolinfo(1)\""
  puts $f " /"
  
  puts $f " &option_param"
  puts $f "   parmformat   = \"$enopt(parmformat)\""
  puts $f "   pbc          = .$enopt(pbc)."
  puts $f "   calc_vdw     = .$enopt(calc_vdw)."
  puts $f "   calc_elec    = .$enopt(calc_elec)."
  puts $f "   calc_torque  = .$enopt(calc_torque)."
  puts $f "   calc_siteesp = .$enopt(calc_siteesp)."
  puts $f "   vdw          = \"$enopt(vdw)\""
  puts $f "   elec         = \"$enopt(elec)\""
  puts $f "   mode         = \"$enopt(mode0)\" \"$enopt(mode1)\""
  puts $f "   rljcut       = $enopt(rljcut)"
  puts $f "   relcut       = $enopt(relcut)"
  puts $f "   pme_alpha    = $enopt(pme_alpha)"
  puts $f "   pme_grids    = $enopt(pme_grids)"
  puts $f "   pme_rigid    = .$enopt(pme_rigid)."
  puts $f "   pme_dual     = .$enopt(pme_dual)."
  puts $f " /"
  close $f

  if {!$common(prep_only)} {
    puts "Interaction Energy is calculated with ANATRA fortran program:"
    puts "$enfort ..."
    puts "=== INPUT ==="
    set content [exec cat $feninp]
    puts $content
    puts "============="
    exec $enfort $feninp >& $fenout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fenout]
    puts $content

    exec rm -f $feninp $fenout $fdcdtmp(0) $fdcdtmp(1) $fmolinfo(0) $fmolinfo(1)
  }

  puts "=============="
  puts ">> Finished"

  exit
}

