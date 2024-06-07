calIndMat <- function(celltype.vector) {
  cell.types <- as.character(unique(celltype.vector))
  ctMat <- lapply(cell.types, function(ii) {
    as.numeric(celltype.vector == ii)
  }) %>% do.call(cbind, .)
  colnames(ctMat) <- cell.types
  rownames(ctMat) <- names(celltype.vector)
  return(ctMat)
}


calRSSMat <- function(rasMat, ctMat){
  rssMat <- pbapply::pblapply(colnames(rasMat), function(i) {
    sapply(colnames(ctMat), function(j) {
      suppressMessages(
        1 - philentropy::JSD(rbind(rasMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
      )
    })
  }) %>% do.call(rbind, .)
  rownames(rssMat) <- colnames(rasMat)
  colnames(rssMat) <- colnames(ctMat)
  return(rssMat)
}

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

DimPlot2 <- function(seu, reduction = "umap", group.by = "celltype", group.highlight=NULL, regulon=NULL, threshold=NULL, text.size = 4) {
  dim.1 <- paste(reduction, 1, sep = "_")
  dim.2 <- paste(reduction, 2, sep = "_")
  if (is.null(regulon)){
    data <- FetchData(seu, vars = c(dim.1, dim.2, group.by))
    data$pt.col <- ifelse(data[[group.by]] == group.highlight, "red", "#DFDFDF")
    data$pt.size <- ifelse(data[[group.by]] == group.highlight, 0.2, 0.1)
    title <- paste0(group.highlight)
    col.title <- "red"
  } else {
    data <- FetchData(seu, vars = c(dim.1, dim.2, regulon))
    mean.val <- mean(data[, regulon])
    sd.val <- sd(data[, regulon])
    if (is.null(threshold)) {
      threshold <- mean.val + 2*sd.val
      threshold <- round(threshold, 2)
    }
    data$pt.col <- ifelse(data[, regulon] > threshold, "#006464", "#DFDFDF")
    data$pt.size <- ifelse(data[, regulon] > threshold, 0.2, 0.1)
    title <- paste0("Regulon: ", regulon)
    title <- glue::glue("{title}\nThreshold: {threshold}")
    col.title = "#006464"
  }
  ggplot(data, aes(get(dim.1), get(dim.2))) +
    geom_point(size=data$pt.size, color=data$pt.col) +
    theme_bw(base_size = 12) +
    ggtitle("") +
    xlab(dim.1) + ylab(dim.2) +
    annotate("text",x=Inf,y=Inf,hjust=1.1,vjust=1.5,label=title,color=col.title,size=text.size) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black")
    )
}
