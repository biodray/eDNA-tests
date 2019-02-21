# Function to remove the first argument of ref fasta file

# when taxonomy, keep element 7

keep_sp <- function(COL){
  if(str_detect(COL[1], pattern = ";")){
    res <- sapply(str_split(COL, pattern=";"), `[`, 7)
  } else {
    res <- paste(sapply(str_split(COL, pattern=" "), `[`, 2),
                 sapply(str_split(COL, pattern=" "), `[`, 3),
                 sep = " ")
  }
  
  if(is.na(res[1])){
    stop("Il y a deja seulement le nom de l'espece, action non effectuee", call. = F)
    
    
  }
  return(res)
  
}


# Then the exact contrary

keep_ID <- function(COL){
  if(str_detect(COL[1], pattern = ";")){
      stop("Il n'y a plus l'ID, action non effectuee", call. = F)
  } else {
    res <- sapply(str_split(COL, pattern=" "), `[`, 1)
  }

  return(res)
  
}
