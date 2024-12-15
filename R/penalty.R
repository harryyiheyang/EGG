soft=function(a,b){
  c=abs(a)-b
  c[c<0]=0
  c=c*sign(a)
  return(c)
}

scad=function(a,lam,ga=3.7){
  b=abs(a)
  z=soft(a,lam)
  z[which(b>(2*lam)&b<=(ga*lam))]=soft(a[which(b>(2*lam)&b<=(ga*lam))],ga*lam/(ga-1))/(1-1/(ga-1))
  z[which(b>(ga*lam))]=a[which(b>(ga*lam))]
  return(z)
}

mcp=function(a,lam,ga=3){
  b=abs(a)
  z=soft(a,lam)/(1-1/ga)
  z[which(b>(ga*lam))]=a[which(b>(ga*lam))]
  return(z)
}
