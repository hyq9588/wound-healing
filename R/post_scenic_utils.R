PlotRegulonRank <- function(rssMat, cell.type, topn=5) {
  data <- data.frame(
    Regulons = 1:nrow(rssMat),
    RSS = sort(rssMat[, cell.type], decreasing = T),
    label = sub("(+)", "", names(sort(rssMat[, cell.type], decreasing = T)), fixed = T)
  )

  data$pt.col <- ifelse(data$Regulons <= topn, "#007D9B", "#BECEE3")
  data <- head(data, n=200)
  data.label <- head(data, n=topn)

  ggplot(data, aes(Regulons, RSS)) +
    geom_point(size=3, color=data$pt.col) +
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = data.label, aes(Regulons, RSS, label=label), size=4) +
    ggtitle(cell.type) + ylab("Specificity score") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5)
    )
}

DimPlot2 <- function(cell.info, dim.1="tSNE_1", dim.2="tSNE_2", group.by = "celltype", group.by.name=NULL, regulon=NULL) {
  if (is.null(regulon)){
    data <- cell.info[, c(dim.1, dim.2, group.by)]
    data$pt.col <- ifelse(data[[group.by]] == group.by.name, "red", "#DFDFDF")
    data$pt.size <- ifelse(data[[group.by]] == group.by.name, 0.2, 0.1)
    title <- paste0(group.by.name)
    col.title <- "red"
  } else {
    data <- cell.info[, c(dim.1, dim.2, regulon)]
    data$pt.col <- ifelse(data[, regulon], "#006464", "#DFDFDF")
    data$pt.size <- ifelse(data[, regulon], 0.2, 0.1)
    title <- paste0("Regulon: ", regulon)
    col.title = "#006464"
  }
  ggplot(data, aes(get(dim.1), get(dim.2))) +
    geom_point(size=data$pt.size, color=data$pt.col) +
    theme_bw(base_size = 12) +
    ggtitle("") +
    xlab(dim.1) + ylab(dim.2) +
    annotate("text",x=Inf,y=Inf,hjust=1.1,vjust=1.5,label=title,color=col.title,size=6) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black")
    )
}
