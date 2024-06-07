`%notin%` <- Negate(`%in%`)

RegNetVis <- function(regulators, edge.list, nodes.show = NULL, size.by = "num_target", topN = 10){
  regulators <- arrange(regulators, cluster, desc(num_target))
  regulators$rank <- factor(regulators$TF, levels = regulators$TF)
  regulators <- regulators %>%
    group_by(cluster) %>%
    mutate(rank_within_cluster = rank(-num_target, ties.method = "first")) %>%
    ungroup()

  to_label <- as.vector(regulators[regulators$rank_within_cluster <= topN,]$TF)
  regulators.plot <- subset(regulators, !is.na(regulators$cluster))
  edge.list.plot <- edge.list
  if(!is.null(nodes.show)) {
    regulators.plot <- subset(regulators.plot, TF %in% nodes.show)
    edge.list.plot <- subset(edge.list.plot, TF %in% nodes.show & target %in% nodes.show)
  }
  ggplot() +
    geom_segment(data = edge.list.plot,
                 mapping = aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 alpha = 0.05) +
    geom_point(data = regulators.plot,
               mapping = aes(x = fr1, y = fr2, color = cluster, size = get(size.by))) +
    guides(size = guide_legend(title = size.by)) +
    ggrepel::geom_text_repel(data = regulators.plot %>% subset(TF %in% to_label),
                             mapping = aes(x = fr1, y = fr2, color = cluster, label = TF),
                             max.overlaps=Inf, show.legend = F) +
    ggsci::scale_fill_d3() +
    ggsci::scale_color_d3() +
    scale_size_area(max_size = 4)  +
    coord_fixed() +
    theme_void() +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
}

RegModuleBar <- function(regulators, topN = 10) {
  regulators <- arrange(regulators, cluster, desc(num_target))
  regulators$rank <- factor(regulators$TF, levels = regulators$TF)

  regulators.plot <- regulators %>%
    group_by(cluster) %>%
    slice_head(n = topN) %>%
    ungroup()

  regulators.plot <- subset(regulators.plot, !is.na(regulators.plot$cluster))

  ggplot(regulators.plot) +
    geom_bar(aes(x = rev(rank), y = num_target, fill = cluster), stat = "identity") +
    geom_text(aes(x = rev(rank), y = 1, label = TF), hjust = 0) +
    coord_flip() +
    expand_limits(y = -100) +
    labs(x = "") +
    ggsci::scale_fill_d3() +
    theme_classic(base_size = 15) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())
}

RegulonGraphVis <- function(data, tf.show, targets.show, edge.alpha = .6, colors = c("grey", "orange"),
                            edge.color = "lightblue", layout = "stress", prop = 0.05, n = NULL) {
  edges <- data[, c("TF", "target", "importance")]
  names(edges)[1:2] <- c("from", "to")

  edges <- edges %>%
    filter(from %in% tf.show) %>%
    group_by(from) %>%
    arrange(desc(importance))

  if (!is.null(prop)) {
    edges <- edges %>% slice_head(prop = prop)
  } else {
    edges <- edges %>% slice_head(n = n)
  }

  targets.show.df <- edges %>%
    filter(to %in% targets.show)

  nodes <- data.frame(name = unique(union(edges$from, edges$to)))
  nodes$annot <- ifelse(nodes$name %in% unique(data$TF), "TF", "target")
  g <- tbl_graph(nodes = nodes, edges = edges)

  lw.breaks <- quantile(edges$importance, c(.01, .5, .99))
  lw.ranges <- c(.3, 2)

  ggraph(g, layout = layout) +
    geom_edge_link(aes(width = importance), alpha = edge.alpha, color = edge.color) +
    geom_node_point(aes(filter = name %in% targets.show.df$to), size = 5, color = "black") +
    geom_node_point(aes(color = annot, size = annot)) +
    geom_node_text(aes(filter = annot == "TF" & name %notin% targets.show, label = name), size = 4, color = "black", repel = T) +
    geom_node_text(aes(filter = annot == "TF" & name %in% targets.show, label = name), size = 4, color = "red", repel = T) +
    geom_node_text(aes(filter = name %in% targets.show.df$to & name %notin% tf.show, label = name), size = 3, color = "red", repel = T) +
    scale_color_manual(values = colors) +
    scale_size_manual(values = c(4, 10)) +
    scale_edge_width_continuous(breaks = lw.breaks, labels = round(lw.breaks, 1), range = lw.ranges) +
    guides(color = guide_legend(title = "", override.aes = list(size = 4, alpha = 1)),
           size = guide_legend(title = "")) + # not work yet
    theme_graph()
}
