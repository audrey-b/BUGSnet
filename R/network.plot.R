#' Network Plot
#' @description Produces network plot where nodes represent treatments and edges represent direct
#' evidence (e.g an RCT) comparing treatments.
#' @param data.nma A data object produced by \code{data.prep()}
#' @param node.scale Size of the nodes (default=5)
#' @param edge.scale Thickness of the edges (default=2)
#' @param flag Used to highlight direct comparisons to a particular treatment (optional.
#' Set this value to treatment of interest and it will highlight, in red, all of the edges
#' going into this treatment.
#' @param label.offset1 Node label location (x-axis) relative to node. Default=0
#' @param label.offset2 Node label location (y-axis) relative to node. Default=1
#' @examples
#' # use default settings
#' network.plot(my.slr)
#' 
#' # Highlight all direct comparisons with Placebo. Adjust node and edge size, centre node labels
#' network.plot(my.slr, node.scale=4, edge.scale=1.5, flag="Placebo", label.offset1=0, label.offset2=0)

net.plot <- function(data.nma,
                         node.scale=5, 
                         edge.scale=2, flag=NULL, 
                         label.offset1=0, label.offset2=1, 
                         graph.scale=T,
                         node.lab.cex = 1,
                         edge.lab.cex = 1,
                         node.colour = "#f69c54",
                         edge.colour = "#4a5b71",
                         edge.lab.colour = "blue",
                         flag.edge.colour = "lightpink",
                         ...) {
  

 edgesANDnodes <- network.structure(data.nma)
 edges <- edgesANDnodes[[1]]
 nodes <- edgesANDnodes[[2]]
 
 names(nodes)[names(nodes) == data.nma$varname.t] <- "trt"
 
 net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
 
 
## Source for function to offset node labels...
## https://stackoverflow.com/questions/23209802/placing-vertex-label-outside-a-circular-layout-in-igraph
 
 radian.rescale <- function(x, start=0, direction=1) {
   c.rotate <- function(x) (x + start) %% (2 * pi) * direction
   c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
 }
 
  lab.locs <- radian.rescale(x=1:nrow(nodes), direction=-1, start=0)
  
  lab.offset <- data.frame(x = abs(lab.locs), y=label.offset1, z = label.offset2) %>%
    mutate(offset = ifelse(x < pi/6 | x > 5*pi/6 & x < 7*pi/6 | x > 11*pi/6, y, y/z)) %>%
    select(offset) %>% pull 
      

## -----------------------------------------------------------------
  # following can be used in the future to specify a 
  # star layout with a single var in the middle
  # layout=layout_as_star(net, center = V(net)$trt==center.trt),

  if(!is.null(flag)) {
    
   inc.edges <- incident(net, V(net)[trt==flag], mode="all")
   
   ecol <- rep("grey", ecount(net))
   ecol[inc.edges] <- flag.edge.colour
   vcol <- rep("grey", vcount(net))
   vcol[V(net)$trt==flag] <- node.colour
   
   plot(net, 
        vertex.size=node.scale*V(net)$node.weight,
        edge.width=edge.scale*E(net)$edge.weight,
        
        vertex.color=vcol,
        vertex.frame.color=vcol,
        vertex.label.cex = node.lab.cex,
        edge.color=ecol,
        
        vertex.label=V(net)$trt,
        vertex.label.color="black",
        vertex.label.family="sans",
        
        
        layout= layout_in_circle(net),
        edge.label= ifelse(ecol==flag.edge.colour, E(net)$study, NA),
        edge.label.family="sans",
        edge.label.cex=edge.lab.cex,
        edge.label.color=edge.lab.colour,
        edge.label.dist=0,
        edge.label.deg=0,
        
        vertex.label.dist=lab.offset,
        vertex.label.degree=lab.locs) 
   }
   else if (graph.scale==F) {
     plot(net, 
          vertex.size=node.scale,
          edge.width=edge.scale,
          vertex.label=V(net)$trt, 
          vertex.label.color="black",
          vertex.color=node.colour,
          vertex.frame.color=node.colour,
          vertex.label.cex = node.lab.cex,

          edge.color=edge.colour,
          
          vertex.label.family="sans",
          
          layout= layout_in_circle(net),
          vertex.label.dist=lab.offset,
          vertex.label.degree=lab.locs,
          rescale = FALSE)
     
   }
 else
     
 plot(net, 
      vertex.label=V(net)$trt, 
      edge.width=edge.scale*E(net)$edge.weight,
      vertex.size=node.scale*V(net)$node.weight,
      
      vertex.label.color="black",
      vertex.color=node.colour,
      vertex.frame.color=node.colour,
      edge.color=edge.colour,
      
      vertex.label.family="sans",

      layout= layout_in_circle(net),
      vertex.label.dist=lab.offset,
      vertex.label.degree=lab.locs,
      
      vertex.label.cex = node.lab.cex)
  
}


