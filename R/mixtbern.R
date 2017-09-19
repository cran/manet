#' Single allocation clustering in networks
#'
#' This function infers K single allocation cluster for actor-event network data.
#' @param y A n x d actor-event adjacency matrix, whereby y_ij is 1 if actor i attended event j -- 0 otherwise.
#' @param K Number of single clusters. Default is set to 4.
#' @param maxT Number of MCMC iterations. Default is set to 5000.
#' @param seed Random seed. Default is 1.
#' @param verbose Set to TRUE if you want to see the steps of the MCMC iterations. Defaults is FALSE.
#' @keywords actor-event, mixture models, networks
#' @export
#' @import MCMCpack combinat
#' @return A manet object consisting of a list with five outputs:
#'\itemize{
#' \item{p.allocation.chain}{ A maxT x n x K array with the posterior probabilities of allocation to the heir clusters.}
#' \item{p.event.chain}{ A maxT x K x d array with the cluster - posterior probabilities of attendance to events.}
#' \item{p.community.chain}{ A maxT x K matrix with the heir cluster proportions.}
#' \item{adj} {The original adjacency matrix.}
#' \item{proc.time}{ The computational time.}
#' }
#' @examples
#' data(deepsouth)
#' ds<-mixtbern(deepsouth,K=2,maxT=100)
#' plot(ds)
#' summary(ds)
mixtbern=function(y,K=4,maxT=5000,seed=1,verbose=FALSE)
{
  ptm <- proc.time()
  set.seed(seed)
  d=ncol(y) ## number of events
  n=nrow(y) ## number of units
  z=matrix(0,n,K) ## original allocation vectors


  # initialize alpha
  alpha.chain=matrix(0,maxT,K) ## canvas
  alpha.prior=rep(1,K) ## prior
  alpha=rep(1/K,K) ## vector

  # initialize p

  pi.greco=matrix(0,K,d) ## matrix of attendance probs
  pi.greco.chain=array(0,c(maxT,K,d)) ## canvas
  pi.greco.prior=array(rep(1,K*d*2),c(K,d,2)) ## prior on attendance probs
  for (j in 1:d) for(k in 1:K) pi.greco[k,j]=rbeta(1,1,1) ## initialize with priors
  index=rep(0,n)
  for (i in 1:n) {index[i]=sample(1:K,1,prob=alpha)
  z[i,index[i]]=1}

  prob.z.chain=array(0,c(maxT,n,K)) ### canvas for post.prob of allocation
  z.chain=array(0,c(maxT,n,K))

  for (tt in 1:maxT) {

    ## Step 1. Sample z
    prob.z=matrix(0,n,K)

    for(k in 1:K){
      f.z.y=rep(0,n)
      for(j in 1:d) {
        f.z.y=f.z.y+apply(as.matrix(y[,j]),1,dbinom,1,prob=pi.greco[k,j],log=TRUE)
      }
      prob.z[,k]=log(alpha[k])+f.z.y
    }
    prob.z=exp(prob.z)
    prob.z=prob.z/rowSums(prob.z)
    prob.z=ifelse(is.na(prob.z),1/K,prob.z)
    z=t(apply(prob.z,1,function(x) {rmultinom(1,1,prob=x)}))

    ## Step 3. Retrieve z and n_k
    n_k=colSums(z)

    # Step 5. Sample alpha_star
    alpha=rdirichlet(1,n_k+alpha.prior)


    # Step 6. Sample pi.greco
    for(k in 1:K){
      for (j in 1:d){
        pi.greco[k,j]<-rbeta(1,y[,j]%*%z[,k]+pi.greco.prior[k,j,1],sum(z[,k])-y[,j]%*%z[,k]+pi.greco.prior[k,j,2])
      }
    }


    ## saving chains into canvases
    prob.z.chain[tt,,]=prob.z
    pi.greco.chain[tt,,]=pi.greco
    alpha.chain[tt,]=alpha
    z.chain[tt,,]=z

    ## print iteration number
    if(verbose)
    print(tt)
  }


  out=(list(p.allocation.chain=prob.z.chain,
            p.event.chain=pi.greco.chain,
            p.community.chain=alpha.chain,
            parent.heir.clusters = diag(rep(1,K)),
            adj=y,
            proc.time=proc.time()-ptm))
  class(out)<-"manet"
  return(out)
}
