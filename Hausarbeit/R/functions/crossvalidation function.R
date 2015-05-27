## Crossvalidation function


createFolds = function(data, k){
  
  #randomly shuffle samples
  if(is.data.frame(data)) {
    index = sample(1:nrow(data))
  } else {
    index = sample(1:length(data))
  }
  
  folds = split(index, 1:k)
  
  
  learning.list = vector("list", length(folds))
  test.list = vector("list", length(folds))
  
  for(i in 1:length(folds)){
    
    learning.data = data[-folds[[i]],]
    test.data = data[folds[[i]],]
    
    learning.list[[i]] = learning.data
    test.list[[i]] = test.data
    
  }
  
  
  
  
  return(list(folds = folds, learning = learning.list, test = test.list))
}





pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
