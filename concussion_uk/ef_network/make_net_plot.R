library(ggnetwork)
times_ten <- scales::trans_new(
  name = "times_ten",
  transform = function(x) abs(x * 10),
  inverse = function(x) abs(x / 10)
)

make_net_plot <- 
  function(net, 
           node_labels = NULL, 
           node_size = 45,
           node_label_col = "black",
           layouts = "circle",
           text_size_node = 12,
           text_size_edge = 13,
           title_size = 24,
           font = "Times",
           title = "")
  {
      
      if (is.null(node_labels)) node_labels <- colnames(net$Y)
      adj <- net$pcor_mat
      edge_weights_unfrmtd <- round(adj[upper.tri(adj)], 2)
      edge_weights <- format(edge_weights_unfrmtd, nsmall = 2)
      
      net_summ <- summary(net)
      lower_bounds_unfrmtd <- round(net_summ$dat_results$Cred.lb, 2)
      upper_bounds_unfrmtd <- round(net_summ$dat_results$Cred.ub, 2)
      
      lower_bounds <- trimws(format(lower_bounds_unfrmtd, nsmall = 2))
      upper_bounds <- trimws(format(upper_bounds_unfrmtd, nsmall = 2))
      
      
      
      # --- Network Attributes ---
      adj_net <- network::as.network(adj, directed = FALSE)
      network::set.vertex.attribute(adj_net, attrname = "node_labs", value = node_labels)
      network::set.edge.attribute(adj_net, attrname = "weights", value = edge_weights_unfrmtd)
      network::get.edge.attribute(adj_net, 'weights')
      
      edge_labels <- paste0(edge_weights, "\n(", lower_bounds, ", ", upper_bounds, ")", sep = "")
      network::set.edge.attribute(adj_net, attrname = "edge_labs", value = edge_labels)
      
      adj_gg <- ggnetwork(adj_net, layout = layouts) 
      
      # --- Network Plot --- 
      net_plot <-
        ggplot(adj_gg,
             aes(x = x, 
                 y = y, 
                 xend = xend, 
                 yend = yend)) +
        geom_edges(aes(size = weights),
                   col = "gray70") +
        geom_nodes(size = node_size,
                   shape = 21,
                   col = "black",
                   fill = "white") +
        geom_nodetext(aes(label = node_labs), 
                      col = node_label_col, 
                      family = font,
                      size = text_size_node) +
        geom_edgetext(aes(label = edge_labs), 
                      family = font,
                      size = text_size_edge,
                      label.padding = unit(0.5, "lines")) +
        scale_size_identity(trans = times_ten) +
        ggnetwork::theme_blank(base_size = title_size, 
                               base_family = font) +
        theme(legend.position = "",
              plot.margin = unit(c(0, 0, 0, 0), "cm"),
              plot.title = element_text(face = "bold",
                                        hjust = 0.1)) +
        ggthemes::scale_fill_few() +
        coord_cartesian(clip = "off")
      
      return(net_plot)
  }

