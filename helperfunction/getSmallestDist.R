getSmallestDist = function(x) {
  if(length(x)>1){
    y = c(x[2:length(x)], max(x)+1000000)
    z = min(y-x)
  }else{
    z=1000000
  }
  return(z)
}