package provide anatra_common 1.0

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_common
#! @brief         define common paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_commoninfo {} {

  global common 

  set common(prep_only)  false;list

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure  read_commoninfo
#! @brief     read common parameters 
#! @authors   KK
#! @param[in] arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_commoninfo {arglist} {

  global common 

  set common(prep_only)   [parse_arguments $arglist "-prep_only"   "value" $common(prep_only)];list

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure  show_commoninfo 
#! @brief     show common parameters 
#! @authors   KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_commoninfo {} {

  global common 
 
  puts "<< common info >>"
  puts "prep_only  = $common(prep_only)"

}
