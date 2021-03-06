#' Multiple allocation clustering of actor-event networks
#'
#' This function infers K multiple allocation cluster for actor-event network data.
#' @param y A n x d actor-event adjacency matrix, whereby y_ij is 1 if actor i attended event j -- 0 otherwise.
#' @param K Number of multiple clusters. Default is set to 2.
#' @param maxT Number of MCMC iterations. Default is set to 5000.
#' @param seed Random seed. Default is 1.
#' @param link Method to combine the parameters of the parent clusters into the parameter for the heir cluster. Default is "min". The alternative is "max".
#' @param verbose Set to TRUE if you want to see the steps of the MCMC iterations. Defaults is FALSE.
#'@return A manet object consisting of a list with five outputs:
#'\itemize{
#' \item{p.allocation.chain}{ A maxT x n x 2^K array with the posterior probabilities of allocation to the heir clusters.}
#' \item{p.event.chain}{ A maxT x K x d array with the cluster - posterior probabilities of attendance to events.}
#' \item{p.community.chain}{ A maxT x 2^K matrix with the heir cluster proportions.}
#' \item{parent.heir.cluster}{ A 2^K x K matrix, which indicates the relationship between parent and heir clusters.}
#' \item{adj} {The original adjacency matrix.}
#' \item{proc.time}{ The computational time.}
#' }
#' @keywords actor-event, mixture models, networks
#' @export
#' @import MCMCpack combinat
#' @importFrom stats dbinom rbeta rmultinom
#' @examples
#' data(deepsouth)
#' ds<-manet(deepsouth,K=2,maxT=100)
#' plot(ds)
#' summary(ds)
manet=function(y,K=2,maxT=5000,seed=1,link="min",verbose=FALSE)
{
  ptm <- proc.time()
  set.seed(seed)
  d=ncol(y) ## number of events
  n=nrow(y) ## number of units
  K_star=2^K ## number of clusters
  u=z_ext(0:1,K) ## connection matrix
  prim=which(rowSums(u)==1) ### quali sono primary
  ## corrispondenza fra i primary e l'ordinamento k=1,2,3,...,K
  which.prim=rep(0,K)
  for(k in 1:K) which.prim[k]=which(u[prim[k],]==1)
  ###
  z=matrix(0,n,K) ## original allocation vectors
  z_star=matrix(0,n,K_star) ## multiple allocation vectors


  # initialize alpha_star
  alpha_star.chain=matrix(0,maxT,K_star) ## canvas
  alpha_star=rep(0,K_star) ## star vector
  alpha_star.prior=rep(1,K_star)
  alpha_star=rdirichlet(1,alpha_star.prior)

  # initialize p

  pi.greco=matrix(0,K,d) ## matrix of attendance probs
  pi.greco.chain=array(0,c(maxT,K,d)) ## canvas
  pi.greco.prior=array(rep(1,K*d*2),c(K,d,2)) ## prior on attendance probs
  for (j in 1:d) for(k in 1:K) pi.greco[k,j]=rbeta(1,1,1) ## initialize with priors
  index=rep(0,n)
  for (i in 1:n) {index[i]=sample(1:K_star,1,prob=alpha_star)
  z_star[i,index[i]]=1}
  z=u[index,]
  s=array(0,c(d,n,K)) ## auxiliary variable
  for(j in 1:d){
    for (i in 1:n){
      if(sum(z[i,])<=1) s[j,i,]=z[i,] else s[j,i,which.max(z[i,]*pi.greco[,j])]=1
    }
  }

  prob.z_star.chain=array(0,c(maxT,n,K_star)) ### canvas for post.prob of allocation
  z_star.chain=array(0,c(maxT,n,K_star))

  for (tt in 1:maxT) {

    ## Step 1. Sample z
    prob.z_star=matrix(0,n,K_star)

    for(h in 1:K_star){
      f.z.y=rep(0,n)
      for(j in 1:d) {
        prbl=ifelse(sum(u[h,])==0,0.00000001,psi.fun(pi.greco[,j],u[h,],link))
        f.z.y=f.z.y+apply(as.matrix(y[,j]),1,dbinom,1,prob=prbl,log=TRUE)
      }
      prob.z_star[,h]=log(alpha_star[h])+f.z.y
    }
    prob.z_star=exp(prob.z_star)
    prob.z_star=prob.z_star/rowSums(prob.z_star)
    if(any(is.na(prob.z_star))==1) print("FLAG: NA on prob.z")
    prob.z_star=ifelse(is.na(prob.z_star),1/K_star,prob.z_star)
    z_star=t(apply(prob.z_star,1,function(x) {rmultinom(1,1,prob=x)}))
    index=apply(z_star,1,which.max)


    ## Step 2. Allocate units according to sampled values
    z=u[index,]

    ### Step 6. compute s
    for(j in 1:d){
      for (i in 1:n){
        if(sum(z[i,])<=1){
          s[j,i,]=z[i,]
        }else{
          hlp=rep(0,K)
          if (link=="min") hlp[which.min(pi.greco[,j]^z[i,])]=1 else hlp[which.max(z[i,]*pi.greco[,j])]=1
          s[j,i,]=hlp
        }
      }
    }

    ## Step 3. Retrieve n_star
    n_star=colSums(z_star)

    # Step 4. Sample alpha_star
    alpha_star=rdirichlet(1,n_star+alpha_star.prior)
    # alpha_star=true.param$alpha_star

    # Step 5. Sample pi.greco
    for(k in 1:K){
      for (j in 1:d){
        pi.greco[k,j]<-rbeta(1,y[,j]%*%s[j,,k]+pi.greco.prior[k,j,1],sum(s[j,,k])-y[,j]%*%s[j,,k]+pi.greco.prior[k,j,2])
      }
    }

    ## saving chains into canvases
    prob.z_star.chain[tt,,]=prob.z_star
    pi.greco.chain[tt,,]=pi.greco
    alpha_star.chain[tt,]=alpha_star

    ## print iteration number
    if(verbose)
      print(tt)
  }

  out=(list(p.allocation.chain=prob.z_star.chain,
            p.event.chain=pi.greco.chain,
            p.community.chain=alpha_star.chain,
            parent.heir.clusters=u,
            adj=y,
            proc.time=proc.time()-ptm))
  class(out)<-"manet"
  return(out)

}

########### Functions
## combining scheme
psi.fun<-function(x,u,link){
  if(link=="min") return(min(x^u)) else return(max(x*u))
}


### computes u connection matrix
z_ext <-function(x,nfac){
  nq <- length(x)
  zx <- hcube(rep(nq,nfac))
  zx <- zx[,dim(zx)[2]:1]
  z2 <- matrix(x[zx],dim(zx)[1],dim(zx)[2])
  return(z2)
}
