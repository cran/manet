#'Printing the output from the multiple allocation clustering
#'
#'This function prints the output of the manet function
#'
#'@param x A manet object.
#'@param digits Number of digits. Default is 3.
#'@param ... Additional arguments to the print function.
#'@export
#' @examples
#' data(deepsouth)
#' ds<-manet(deepsouth,K=2,maxT=100)
#' print(ds)
print.manet<-function(x,digits=3, ...){
  manet.obj<-x
  n.digits<-digits
  #allocation
  actor<-(1:nrow(manet.obj$adj))
  alloc<-apply(apply(manet.obj$p.allocation.chain,c(2,3),mean),1,which.max)
  prob.alloc<-round(apply(apply(manet.obj$p.allocation.chain,c(2,3),mean),1,max),n.digits)
  allocation<-cbind(actor,alloc,prob.alloc)

  #subgroup
  cluster.frac<-round(apply(manet.obj$p.community.chain,2,mean),n.digits)
  cluster<-cbind(cluster.frac,1:length(cluster.frac),manet.obj$parent.heir.clusters)

  #event attendance parameters
  event<-(1:ncol(manet.obj$adj))
  prob.event<-round(apply(manet.obj$p.event.chain,c(2,3),mean),n.digits)

  #give names
  colnames(allocation)<-c("actor", "heir cluster","probability")
  nms<-lapply(1:nrow(prob.event),function(x){paste("parent cluster",x)})
  rownames(prob.event)<-nms
  event<-rbind(event,prob.event)
  colnames(cluster)<-c("fraction","heir cluster",nms)

  out<-list(actor.allocations=allocation, event.probabilities=event,allocation.fractions=cluster)
  return(out)
}
