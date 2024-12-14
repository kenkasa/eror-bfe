package provide anatra_selection 1.0

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_selinfo
#! @brief         define selection paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_selinfo {} {

  global seltxt 
  global nselmax 
  
  set nselmax 10000

  set seltxt(nsel) 0
  for {set i 0} {$i < $nselmax} {incr i} {
    set seltxt($i) ""
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure  read_selinfo
#! @brief     read selection parameters 
#! @authors   KK
#! @param[in] arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_selinfo {arglist} {

  global seltxt 
  global nselmax

  set seltxt(nsel) 0  
  for {set i 0} {$i < $nselmax} {incr i} {
    set selopt     [format "-sel%d" $i]
    set seltxt($i) [parse_arguments $arglist "$selopt" "value" $seltxt($i)];list
    if {$seltxt($i) != ""} {
      incr seltxt(nsel)
    }
  } 

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure  show_selinfo 
#! @brief     show selection parameters 
#! @authors   KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_selinfo {} {

  global seltxt 
 
  puts "<< selection info >>"
  set i 0
  while {$seltxt($i) != ""} {
    set sels [format "sel%5d = " $i]
    puts "$sels $seltxt($i)"
    incr i
  }

}

