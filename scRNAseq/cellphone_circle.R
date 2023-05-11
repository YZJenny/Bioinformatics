circle_plot <- function (net, color.use = NULL, from = NULL, to = NULL, bidirection = FALSE, 
          remove.isolate = FALSE, top = 1, 
          weight.scale = FALSE, label.edge = FALSE, edge.curved = 0.2, 
          shape = "circle", layout = in_circle(), vertex.size = 20, 
          margin = 0.2, vertex.label.cex = 1.2, vertex.label.color = "black", 
          arrow.width = 1, arrow.size = 1, edge.label.color = "black", 
          edge.label.cex = 0.5, edge.max.width = 8, alpha.edge = 0.6) 
{
  options(warn = -1)

  thresh <- stats::quantile(net, probs = 1 - top)
  net[net < thresh] <- 0
  if (!is.null(from)) {
    if (!is.numeric(from)) {
      group.names.all <- rownames(net)
      if (is.null(group.names.all)) {
        stop("No rownames exist in the input network matrix. Please define 'from' using numeric vector. ")
      }
      from <- match(from, group.names.all)
      if (sum(is.na(from)) > 0) {
        message("Some input cell group names in 'from' do not exist!")
        from <- from[!is.na(from)]
      }
    }
    if (is.null(to)) {
      to <- 1:ncol(net)
    }
    else {
      if (!is.numeric(to)) {
        to <- match(to, group.names.all)
        if (sum(is.na(to)) > 0) {
          message("Some input cell group names in 'to' do not exist!")
          to <- to[!is.na(to)]
        }
      }
    }
    if (bidirection) {
      from2 <- c(from, to)
      to2 <- c(to, from)
      from.c <- setdiff(1:nrow(net), from2)
      to.c <- setdiff(1:nrow(net), to2)
      net[from.c, ] <- 0
      net[, to.c] <- 0
      net[from, from] <- 0
      net[to, to] <- 0
    }
    else {
      from.c <- setdiff(1:nrow(net), from)
      to.c <- setdiff(1:nrow(net), to)
      net[from.c, ] <- 0
      net[, to.c] <- 0
    }
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(V(g)))
  }
  vertex.size <- vertex.size/max(vertex.size) * 15 + 5
  loop.angle <- ifelse(coords_scale[V(g), 1] > 0, -atan(coords_scale[V(g), 
                                                                     2]/coords_scale[V(g), 1]), pi - atan(coords_scale[V(g), 
                                                                                                                       2]/coords_scale[V(g), 1]))
  V(g)$size <- vertex.size
  V(g)$color <- color.use[V(g)]
  V(g)$frame.color <- color.use[V(g)]
  V(g)$label.color <- vertex.label.color
  V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    E(g)$label <- E(g)$weight
  }
  if (weight.scale == TRUE) {
    E(g)$width <- 0.3 + edge.max.width/(max(E(g)$weight) - 
                                          min(E(g)$weight)) * (E(g)$weight - min(E(g)$weight))
  }else {
    E(g)$width <- 0.3 + edge.max.width * E(g)$weight
  }
  E(g)$arrow.width <- arrow.width
  E(g)$arrow.size <- arrow.size
  E(g)$label.color <- edge.label.color
  E(g)$label.cex <- edge.label.cex
  E(g)$color <- grDevices::adjustcolor(V(g)$color[edge.start[, 
                                                             1]], alpha.edge)
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
                                                        1])] <- loop.angle[edge.start[which(edge.start[, 
                                                                                                       2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(V(g)), direction = -1, 
                               start = 0)
  label.dist <- vertex.size/max(vertex.size) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Arial")
}
