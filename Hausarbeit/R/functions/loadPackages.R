# Load Packages Function
loadPackages = function(packageList){
  if (length(setdiff(packageList, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packageList, rownames(installed.packages())))  
  }
  
  sapply(packageList, require, character.only = T)
  
}