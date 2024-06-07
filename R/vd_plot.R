VDPlot <- function(vd.res, y, group.by, group.show) {
  vd.res <- arrange(vd.res, desc(get(y)))
  vd.res <- subset(vd.res, !is.na(cluster))
  vd.plots <- vd.res %>%
    mutate(pt.size = ifelse(get(group.by) == group.show, 2, .5),
           type = ifelse(get(group.by) == group.show, group.show, "others"),
           rank = 1:nrow(.))
  ggplot(vd.plots, aes(rank, get(y))) +
    geom_point(aes(color = type), size = vd.plots$pt.size) +
    ggrepel::geom_text_repel(data = subset(vd.plots, get(group.by) == group.show),
                             aes(rank, get(y), label = gene), color = "red") +
    scale_color_manual(values = c("red", "grey")) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(x = "Rank", y = sprintf("Fraction of variance across %s", y)) +
    theme_classic(base_size = 16) +
    theme(
      legend.title = element_blank(),
      legend.position = c(1,1),
      legend.justification = c(1,1),
      axis.text = element_text(color = "black")
    )
}

