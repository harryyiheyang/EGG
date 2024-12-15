EGG.subsampling.Pvalue=function(A,B){
  p=dim(A)[2];P=SE=diag(p)
  for(i in 1:p){
    for(j in i:p){
      a=A[i,j];b=sd(B[,i,j])
      SE[i,j]=SE[j,i]=b
      P[i,j]=P[j,i]=pchisq((a/b)^2,1,lower.tail=F)
    }
  }
  return(list(Pvalue=P,SE=SE))
}
