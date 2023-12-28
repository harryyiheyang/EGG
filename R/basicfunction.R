colMedian=function(A){
  a=apply(A,2,median)
  return(a)
}

rowMedian=function(A){
  a=apply(A,1,median)
  return(a)
}

vec=function(a){
  as.vector(a)
}
