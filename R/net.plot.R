#' Network Plot
#' @description Produces network plot where nodes represent treatments and edges represent direct
#' evidence (e.g an RCT) comparing treatments.
#' @param data A \code{BUGSnetData} object produced by \code{data.prep()}.
# @param outcome A string of the outcome variable to use for the plot (missing values are discarded).
#' @param node.scale Size of the nodes (default=5)
#' @param edge.scale Thickness of the edges (default=2).
#' @param graph.scale Whether to make edges and nodes proportionnaly larger with the number of studies/arms. Default is TRUE.
#' @param flag Used to highlight direct comparisons to particular treatments (optional).
#' Set this value to treatment(s) of interest and it will highlight, in red, all of the edges
#' going into the specified treatment(s).
#' @param study.counts If TRUE, prints the number of studies on each edge.
#' @param label.offset1 Node label location (x-axis) relative to node. Default=0
#' @param label.offset2 Node label location (y-axis) relative to node. Default=1
#' @param node.lab.cex Size of node labels
#' @param edge.lab.cex Size of edge labels
#' @param node.colour Node colour (string)
#' @param edge.colour Edge colour (string)
#' @param edge.lab.colour Edge label colour (string)
#' @param flag.edge.colour Color of flagged edges (string)
#' @param layout Specifies how the nodes should be layed out in the graph. Default is "layout_in_circle". See \code{igraph::layout_} 
#' for layout options
#' @param layout.params Additional parameters to be passed to the chosen 'layout' function. See \code{igraph::layout_} for more information
# @param ... extra parameters for plotting
#' @examples
#' 
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(arm.data = diabetes.sim, 
#' varname.t = "Treatment", 
#' varname.s = "Study")
#' 
#' # use default settings
#' net.plot(diabetes.slr)
#' 
#' # Highlight all direct comparisons with Placebo. Adjust node and edge size, centre node labels
#' net.plot(diabetes.slr, 
#' node.scale=4, 
#' edge.scale=1.5, 
#' flag="Placebo", 
#' label.offset1=0, 
#' label.offset2=0)
#' @export
#' @seealso \code{\link{data.prep}}, \code{igraph::layout_}




net.plot <- function(data,
                     #outcome,
                     node.scale=5, 
                     edge.scale=2, 
                     flag=NULL, 
                     study.counts = FALSE,
                     label.offset1=0, 
                     label.offset2=1, 
                     graph.scale=TRUE,
                     node.lab.cex = 1,
                     edge.lab.cex = 1,
                     node.colour = "#f69c54",
                     edge.colour = "grey",
                     edge.lab.colour = "blue",
                     flag.edge.colour = "lightpink",
                     layout = "layout_in_circle",
                     layout.params = NULL) {
  
  # Bind variables to function
  x <- NULL
  y <- NULL
  z <- NULL
  
  if(class(data) != "BUGSnetData")
    stop("\'data\' must be a valid BUGSnetData object created using the data.prep function.")
  
  if (!is.character(layout) || length(layout) != 1)
    stop("\'layout\' must be a character vector of length 1 specifying the name of an igraph layout function")
  
  #To specify the outcome
  # data %<>% data.prep(arm.data=filter(data$arm.data, !is.na(outcome)), 
  #                     varname.t=data$varname.t, 
  #                     varname.s=data$varname.s)
  
  edgesANDnodes <- suppressMessages(network.structure(data))
  edges <- edgesANDnodes[[1]]
  nodes <- edgesANDnodes[[2]]
  
  net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
  
  ## Source for function to offset node labels...
  ## https://stackoverflow.com/questions/23209802/placing-vertex-label-outside-a-circular-layout-in-igraph
  
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(rescale(x, c(0, 2 * pi), range(x)))
  }
  
  lab.locs <- radian.rescale(x=1:nrow(nodes), direction=-1, start=0)
  
  lab.offset <- data.frame(x = abs(lab.locs), y=label.offset1, z = label.offset2) %>%
    mutate(offset = ifelse(x < pi/6 | x > 5*pi/6 & x < 7*pi/6 | x > 11*pi/6, y, y/z)) %>%
    select(offset) %>% pull 
  
  if (graph.scale == TRUE) {
    vsize <- node.scale * V(net)$node.weight
    ewidth <- edge.scale * E(net)$edge.weight
    rscl <- TRUE
  } else {
    vsize <- node.scale
    ewidth <- edge.scale
    rscl <- FALSE
  }
  
  if (!is.null(flag) || study.counts == TRUE) {
    
    if (!is.null(flag)) {
      inc.edges <- incident(net, V(net)[flag[1]], mode="all")
      if (length(flag) > 1) {
        for (i in 2:length(flag))
          inc.edges <- c(inc.edges, inc.edges <- incident(net, V(net)[flag[i]], mode="all"))
      }
      
      ecol <- rep("grey", ecount(net))
      ecol[inc.edges] <- flag.edge.colour
      vcol <- rep("grey", vcount(net))
      vcol[V(net)[flag]] <- node.colour
    } else {
      ecol <- edge.colour
      vcol <- node.colour
    }
    
    if (study.counts == TRUE) {
      elab <- as.character(E(net)$edge.weight)
    } else if (!is.null(flag)) {
      elab <- ifelse((1:length(E(net))) %in% inc.edges, as.character(E(net)$study), NA)
    } else {
      elab <- NA
    }
    
    plot(net, 
         vertex.size=vsize,
         edge.width=ewidth,
         
         vertex.color=vcol,
         vertex.frame.color=vcol,
         vertex.label.cex = node.lab.cex,
         edge.color=ecol,
         
         vertex.label=names(V(net)),
         vertex.label.color="black",
         vertex.label.family="sans",
         
         
         layout=do.call(get(layout, asNamespace("igraph")), c(list(net), layout.params)),
         edge.label= elab,
         edge.label.family="sans",
         edge.label.cex=edge.lab.cex,
         edge.label.color=edge.lab.colour,
         edge.label.dist=0,
         edge.label.deg=0,
         
         vertex.label.dist=lab.offset,
         vertex.label.degree=lab.locs,
         
         rescale = rscl) 
  } else {
    plot(net, 
         vertex.size=vsize,
         edge.width=ewidth,
         vertex.label=names(V(net)), 
         vertex.label.color="black",
         vertex.color=node.colour,
         vertex.frame.color=node.colour,
         vertex.label.cex = node.lab.cex,
         
         edge.color=edge.colour,
         
         vertex.label.family="sans",
         
         layout= do.call(get(layout, asNamespace("igraph")), c(list(net), layout.params)),
         vertex.label.dist=lab.offset,
         vertex.label.degree=lab.locs,
         rescale = rscl)
  }
}


