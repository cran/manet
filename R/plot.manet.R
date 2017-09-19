#' Plotting the output from the multiple allocation clustering.
#'
#' This function plots the output of the manet function.
#' @param x A manet object.
#' @param seed Random seed. Default is 1.
#' @param layout Layout of the network from the igraph package. Default is layout_nicely.
#' @param ... Additional inputs to the igraph function.
#' @return An actor-event network with events as round circles and actors as squared circles with the different colours corresponding to the identified communities.
#' @export
#' @import  igraph
#' @importFrom grDevices rainbow
#' @importFrom graphics plot
#' @examples
#' data(deepsouth)
#' ds<-manet(deepsouth,K=2,maxT=100)
#' plot(ds)

plot.manet<-function(x,seed=1,layout = layout_nicely, ...){
  manet.obj<-x
  set.seed(seed)
  alloc<-apply(apply(manet.obj$p.allocation.chain,c(2,3),sum),1,which.max)
  y<-manet.obj$adj
  edges<-NULL
  n.actor<-nrow(y)
  n.event<-ncol(y)
  types<-rep(c(0,1),c(n.actor,n.event))
  for (i in 1:n.actor){
    for (j in 1:n.event){
      if (y[i,j]==1){
        edges<-c(edges,c(i,n.actor+j))
      }
    }
  }
  g<-make_bipartite_graph(types,edges,directed = FALSE)
  label<-rep(c("square","circle"),c(n.actor,n.event))
  #col<-rep(c("yellow","white"),c(n.actor,n.event))
  col<-rainbow(n=max(alloc),s=.5)[alloc]
  plot(g,vertex.shape=label,vertex.label = c(1:n.actor,1:n.event),vertex.color=col,layout=layout, ...)
  return(g)
}
