package provide anatra_structure 1.0

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_structure
#! @brief         define structure paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_strinfo {} {

  global str 

  set str(stype) "";list
  set str(sfile) "";list

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure  read_strinfo
#! @brief     read structure parameters 
#! @authors   KK
#! @param[in] arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_strinfo {arglist} {

  global str 
 
  set    str(stype) [parse_arguments $arglist "-stype" "value" $str(stype)];list
  set    str(sfile) [parse_arguments $arglist "-sfile" "value" $str(sfile)];list

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure  show_structure 
#! @brief     show structure parameters 
#! @authors   KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_strinfo {} {

  global str 
 
  puts "<< structure info >>"
  puts "stype = $str(stype)"
  puts "sfile = $str(sfile)"
  puts ""
}
