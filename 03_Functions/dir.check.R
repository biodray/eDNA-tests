# Simple function to create a directory if it doesn't existe only

dir.check <- function(directory){
  
  if(!dir.exists(directory)){
     
    dir.create(directory)
  
    }
  
}



