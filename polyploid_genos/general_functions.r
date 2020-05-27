# General purpose functions used by many scripts

add_slash <- function(dir_string){
  # Add '/' to end of a string, particularly directory names that will be
  #  used with paste() to generate file names
  # INPUTS
  # dir_string = character string; ex: '/FULL/DIR/STRING'
  # OUPUT
  # character string
  #   if dir_string already ends with '/', then returns dir_string
  #   if dir_string does not end with '/', returns dir_string ending with '/'
  #########
  final_string <- dir_string
  last_char <- rev(unlist(strsplit(dir_string, split = '')))[1]
  if(last_char != '/'){
    final_string <- paste(final_string, '/', sep = '')
  }
  return(final_string)
}
