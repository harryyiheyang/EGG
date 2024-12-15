reliability.adj=function(X,R,thres=0.8){
p=ncol(R)
total.var=c(1:p)
for(i in 1:p){
total.var[i]=var(X[,i])
}
r=diag(R)
R=cov2cor(R)
reliability=(total.var-r)/total.var
ind=which(reliability<thres)
if(length(ind)>0){
r[ind]=total.var[ind]*(1-thres)
}
r=sqrt(r)
R=t(R*r)*r
return(R)
}
