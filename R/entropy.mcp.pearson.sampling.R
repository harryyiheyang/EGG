#' Estimate genetic network using MCP-penalized entropy loss function and Pearson's r correlaiton matrix
#'
#' The \code{entropy.mcp.spearman.sampling} function estimates the genetic network by applying a penalized entropy loss function using MCP (Minimax Concave Penalty) and Pearson's r correlation matrix. It's designed to work with genetic data, specifically effect size estimates or Z-scores for various exposures across multiple variants.
#'
#' @param BETA A matrix of dimensions mxp, where m is the number of variants and p is the number of exposures. It represents effect size estimates or Z-scores.
#' @param Rnoise The covariance matrix of estimation errors in BETA. If BETA contains Z-scores, Rnoise should be a correlation matrix.
#' @param lamvec A vector of tuning parameters for MCP. Default is c(1:20)/100.
#' @param max.eps The threshold for stopping the algorithm. Default is 0.01.
#' @param max.iter The maximum number of iterations for the algorithm. Since ADMM converges quickly but not with high precision, this should not be too large. Default is 25.
#' @param rho The penalty parameter for the ADMM algorithm. Default is 0.05.
#' @param mineig The minimum matrixEigenvalue for the precision matrix estimate. Default is 0.01.
#' @param subtime The number of times to resample. Default is 100.
#' @param subfrac The fraction of the data to resample each time. Default is 0.5.
#' @param subthres The threshold for stability selection. Default is 0.95.
#' @param alpha The additional parameter for MCP + alpha * Ridge. Default is 0.
#' @param bic.factor An additional penalty on the cross-validation error. Default is 0.5.
#' @param reliability.thres The threshold for rescaling Rnoise to ensure that (var(BETAj)-Rnoisejj)/var(BETAj) is greater than this value. Default is 0.8.
#' @param PenaMatrix A penalty weight matrix to rescale the tuning parameter in MCP. Default is a matrix of ones.
#'
#' @return A list containing:
#'   - Theta: The estimated precision matrix.
#'   - Pvalue: P-values for the entries in Theta.
#'   - K: Stability selection matrix indicating the proportion of subsamples where each entry in Theta was non-zero.
#'   - cv.error: The cross-validation error for each lambda in lamvec.
#'   - R: The Spearman rank correlation matrix of BETA.
#'   - ThetaSE: The standard error of Theta based on subsampling.
#'
#' @examples
#' # Assuming BETA and Rnoise are already defined:
#' result <- entropy.mcp.pearson.sampling(BETA, Rnoise)
#'
#' @details The function performs subsampling and uses the MCP method to estimate a sparse precision matrix (inverse covariance matrix) representing the genetic network. It uses Spearman's rank correlation to handle non-linear relationships and applies stability selection to enhance the reliability of the selected network.
#' @importFrom stats cor cov2cor median pchisq sd var
#' @import glasso
#' @import CppMatrix
#' @export
#'
entropy.mcp.pearson.sampling=function(BETA,Rnoise,lamvec=c(1:20)/100,max.eps=0.01,max.iter=25,rho=0.05,mineig=0.01,subtime=100,subfrac=0.5,subthres=0.95,alpha=0,bic.factor=0.5,reliability.thres=0.8,PenaMatrix="none"){
  m=dim(BETA)[1];p=dim(BETA)[2]
  alpha=alpha*rho
  Rnoise=reliability.adj(X=BETA,R=Rnoise,thres=reliability.thres)
  if(PenaMatrix[1]=="none") PenaMatrix=matrix(1,p,p)
  #########################tuning parameter selection##############################
  subvecerror=matrix(0,length(lamvec),subtime)
  Thetalist=array(NA,c(length(lamvec),subtime,p,p))

  for(j in 1:subtime){
    indsub=sample(m,round(subfrac*m),replace=F)
    BETAS=BETA[indsub,]
    S=t(BETAS)%*%BETAS
    M=Rnoise*dim(BETAS)[1]
    S1=cov2cor1(S-M)
    Theta0=glasso::glasso(S1,0.1)$wi

    BETAS2=BETA[-indsub,]
    S2=t(BETAS2)%*%BETAS2
    M2=Rnoise*dim(BETAS2)[1]
    S2=cov2cor1(S2-M2)
    S2=MCPthreshold(S2,2*sqrt(log(p)/m),ga=3)
    for(i in 1:length(lamvec)){
      Theta=Theta0
      Theta1=Theta0*0
      Delta1=Theta
      Gamma1=Delta1*0
      Delta2=Theta
      Gamma2=Delta2*0
      error=norm(Theta-Theta1,"f")/sqrt(p)
      iter=0
      while(error>max.eps & iter<max.iter){
        Theta1=Theta
        Q=S1+Gamma1-rho*Delta1
        Q=(Q+t(Q))/2
        Theta=(-Q+matrixsqrt(matrixMultiply(Q,Q)+4*rho*diag(p))$w)/(2*rho+alpha)
        Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[i]/rho,ga=3)
        Delta1=matrix(Delta1,p,p)*PenaMatrix
        Gamma1=Gamma1+rho*(Theta-Delta1)
        if(iter %% 5==0 & min(matrixEigen(Theta)$values)<mineig){
          Q=S1+Gamma1+Gamma2-rho*(Delta1+Delta2)
          Theta=(-Q+matrixsqrt(matrixMultiply(Q,Q)+8*rho*diag(p))$w)/(4*rho+alpha)
          Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[i]/rho,ga=3)
          Delta1=matrix(Delta1,p,p)*PenaMatrix
          Gamma1=Gamma1+rho*(Theta-Delta1)
          Delta2=positiveadj(Theta+Gamma2/rho,min.eps=mineig)
          Gamma2=Gamma2+rho*(Theta-Delta2)
        }
        iter=iter+1
        error=norm(Theta-Theta1,"f")/sqrt(p)
      }
      df=(sum(Delta1!=0)-p)/2
      subvecerror[i,j]=entropyloss(S2,(Delta1!=0)*Theta)+(log(p*(p-1)/2)+log(m))/m*df*bic.factor
      Thetalist[i,j,,]=(Delta1!=0)*Theta
    }
  }

  istar=which.min(rowMedian(subvecerror))
  Thetalist=Thetalist[istar,,,]
  K=Thetalist[1,,]*0
  Thetasub=K
  for(i in 1:subtime){
    K=K+(Thetalist[i,,]!=0)/subtime
  }

  S=t(BETA)%*%BETA
  M=Rnoise*dim(BETA)[1]
  S1=cov2cor1(S-M)
  Theta0=glasso::glasso(S1,0.1)$wi
  Theta=Theta0
  Theta1=Theta0*0
  Delta1=Theta
  Gamma1=Delta1*0
  Delta2=Theta
  Gamma2=Delta2*0
  error=norm(Theta-Theta1,"f")/sqrt(p)
  iter=0
  while(error>max.eps & iter<(2*max.iter)){
    Theta1=Theta
    Q=S1+Gamma1-rho*Delta1
    Q=(Q+t(Q))/2
    Theta=(-Q+matrixsqrt(matrixMultiply(Q,Q)+4*rho*diag(p))$w)/(2*rho+alpha)
    Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[istar]/rho,ga=3)
    Delta1=matrix(Delta1,p,p)*PenaMatrix*(K>subthres)
    Gamma1=Gamma1+rho*(Theta-Delta1)
    if(iter %% 5==0 & min(matrixEigen(Theta)$values)<mineig){
      Q=S1+Gamma1+Gamma2-rho*(Delta1+Delta2)
      Theta=(-Q+matrixsqrt(matrixMultiply(Q,Q)+8*rho*diag(p))$w)/(4*rho+alpha)
      Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[istar]/rho,ga=3)
      Delta1=matrix(Delta1,p,p)*PenaMatrix
      Gamma1=Gamma1+rho*(Theta-Delta1)
      Delta2=positiveadj(Theta+Gamma2/rho,min.eps=mineig)
      Gamma2=Gamma2+rho*(Theta-Delta2)
    }
    iter=iter+1
    error=norm(Theta-Theta1,"f")/sqrt(p)
  }
  Theta0=(Delta1!=0)*Theta
  AA=EGG.subsampling.Pvalue(Theta0,Thetalist)
  Pvalue=AA$Pvalue
  ThetaSE=AA$SE
  colnames(Theta0)=rownames(Theta0)=colnames(K)=rownames(K)=colnames(ThetaSE)=rownames(ThetaSE)=colnames(Pvalue)=rownames(Pvalue)=colnames(BETA)

  A=list()
  A$Theta=Theta0
  A$Pvalue=Pvalue
  A$K=K
  A$cv.error=subvecerror
  A$R=S1
  A$ThetaSE=ThetaSE
  return(A)
}
