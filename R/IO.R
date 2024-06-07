LoadpySCENICOutput <- function(regulon.gmt, adj.mat.file) {
  message(glue::glue("Loading {regulon.gmt} ..."))
  tf2target <- clusterProfiler::read.gmt(regulon.gmt)
  tf2target$TF <- sub("\\([0-9]+g\\)", "", tf2target$term)
  colnames(tf2target)[1:2] <- c("regulon", "target")
  tf2target$id <- paste0(tf2target$TF, "-", tf2target$target)
  head(tf2target)
  message(glue::glue("Loading {adj.mat.file} ..."))
  adj <- data.table::fread(adj.mat.file, sep = "\t", header = T)
  adj$id <- paste0(adj$TF, "-", adj$target)
  adj <- subset(adj, id %in% tf2target$id)
  data <- left_join(adj, tf2target, by = c("id", "TF", "target"))
  return(data)
}
