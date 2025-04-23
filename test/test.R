library(igraph)
library(RColorBrewer)
devtools::load_all()
par(mfrow=c(2,3))
data("EAS_Error_COV")
data("EAS_Zscore")
BETA=EASZZ[,-1]
R=EASZR
fit1=entropy.mcp.spearman.sampling(BETA,Rnoise=R,lamvec=c(25:52)/1125,max.eps=0.005,max.iter=25,rho=0.05,mineig=0.01,subfrac=0.5,subthres=0.5,subtime=300,alpha=0)
plot(rowMeans(fit1$cv.error),main="CV Error of EGG")

S=cor(BETA)
n=dim(BETA)[1]
E1=matrix(0,300,30)
B=array(0,c(300,30,20,20))
for(k in 1:300){
  ind=sample(n,0.5*n,replace=F)
  BETA1=BETA[ind,]
  BETA2=BETA[-ind,]
  S1=cor(BETA1)
  S2=cor(BETA2)
  for(i in 1:30){
    fit=glasso::glasso(s=S1,rho=i/1000)
    w=fit$wi
    df=(20+sum(w!=0))/2
    E1[k,i]=entropyloss(A=S2,B=w)
    B[k,i,,]=w
  }
}
istar=which.min(colMeans(E1))
B=B[,istar,,]
K=S*0
for(i in 1:300){
  K=K+(B[i,,]!=0)/300
}
fit=glasso::glasso(s=S,rho=istar/1000)
wcor=fit$wi*(K>0.95)
plot(colMeans(E1),main="CV Error of graphical lasso (Pearson)")

S=spearmancov(BETA)
n=dim(BETA)[1]
E2=matrix(0,300,30)
B=array(0,c(300,30,20,20))
for(k in 1:300){
  ind=sample(n,0.5*n,replace=F)
  BETA1=BETA[ind,]
  BETA2=BETA[-ind,]
  S1=spearmancov(BETA1)
  S2=spearmancov(BETA2)
  for(i in 1:30){
    fit=glasso::glasso(s=S1,rho=i/500)
    w=fit$wi
    df=(20+sum(w!=0))/2
    E2[k,i]=entropyloss(A=S2,B=w)
    B[k,i,,]=w
  }
}
istar=which.min(colMeans(E2))
B=B[,istar,,]
K=S*0
for(i in 1:300){
  K=K+(B[i,,]!=0)/300
}
fit=glasso::glasso(s=S,rho=istar/500)
wspearman=fit$wi*(K>0.95)
plot(colMeans(E2),main="CV Error of graphical lasso (Spearman)")

ThetaEUR=-fit1$Theta*(fit1$K>0.99)
KEUR=fit1$K
nam=rownames(ThetaEUR)
colnames(ThetaEUR)=rownames(ThetaEUR)=nam
nam=sort(nam);
ThetaEUR=ThetaEUR[nam,nam]
gEUR <- graph_from_adjacency_matrix(ThetaEUR, weighted = TRUE, diag = FALSE, mode = "undirected")
V(gEUR)$size <- 12
E(gEUR)$color <- ifelse(E(gEUR)$weight > 0, "black", "grey80")
plot(gEUR, layout = layout.circle,
     vertex.size = V(gEUR)$size,
     vertex.label = V(gEUR)$name,
     vertex.label.cex = 2,
     vertex.label.color="black",
     edge.color = E(gEUR)$color,
     vertex.frame.color=V(gEUR)$color,
     edge.width=2,main="Network estimate of EGG")

wspearman=-as.matrix(wspearman)
rownames(wspearman)=colnames(wspearman)=colnames(ThetaEUR)
gEUR <- graph_from_adjacency_matrix(wspearman, weighted = TRUE, diag = FALSE, mode = "undirected")
V(gEUR)$size <- 12
V(gEUR)$color="#2EC4B6"
E(gEUR)$color <- ifelse(E(gEUR)$weight > 0, "black", "grey80")
plot(gEUR, layout = layout.circle,
     vertex.size = V(gEUR)$size,
     vertex.label = V(gEUR)$name,
     vertex.label.cex = 2,
     vertex.label.color="black",
     edge.color = E(gEUR)$color,
     vertex.frame.color=V(gEUR)$color,
     edge.width=2,main="Network estimate of graphical lasso (Pearson)")

wcor=-as.matrix(wcor)
rownames(wcor)=colnames(wcor)=colnames(ThetaEUR)
gEUR <- graph_from_adjacency_matrix(wcor, weighted = TRUE, diag = FALSE, mode = "undirected")
V(gEUR)$size <- 12
V(gEUR)$color="#ee6e9f"
E(gEUR)$color <- ifelse(E(gEUR)$weight > 0, "black", "grey80")
plot(gEUR, layout = layout.circle,
     vertex.size = V(gEUR)$size,
     vertex.label = V(gEUR)$name,
     vertex.label.cex = 2,
     vertex.label.color="black",
     edge.color = E(gEUR)$color,
     vertex.frame.color=V(gEUR)$color,
     edge.width=2,main="Network estimate of graphical lasso (Spearman)")
