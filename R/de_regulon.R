DERegulon <- function(seu, celltype, group, test.use = "wilcox") {
  seu$new.group <- paste(seu[[celltype, drop=T]], seu[[group, drop=T]], sep = "_")
  Idents(seu) <- factor(seu$new.group)

  celltypes <- levels(seu[[celltype, drop=T]])
  groups <- levels(seu[[group, drop=T]])
  seu2 <- subset(seu, group == groups[1])
  Idents(seu2) <- seu2[[celltype, drop=T]]
  baseline.levels <- AverageExpression(seu2, assays = "AUCell")$AUCell

  de.list <- lapply(celltypes, function(ct) {
    message(glue::glue("processing {ct} ..."))
    ct1 <- paste(ct, groups[1], sep = "_")
    ct2 <- paste(ct, groups[2], sep = "_")
    de <- FindMarkers(seu, ident.1 = ct1, ident.2 = ct2, test.use = "wilcox", fc.name = "avg_diff", logfc.threshold = 0)
    de$change <- ifelse(de$avg_diff > 0,
                        paste("up in", groups[1]),
                        paste("up in", groups[2]))
    de$avg_diff <- -de$avg_diff
    de$diff_rate <- de$avg_diff / baseline.levels[rownames(de), ct]
    de$group <- ct
    return(de)
  })
  names(de.list) <- celltypes
  de.list
}

#' regulon activity ~ (0,1)
format_change.df <- . %>%
  mutate(change = case_when(p_val_adj < 1e-6 & avg_diff > 0.3  ~ '+++',
                            p_val_adj < 1e-6 & avg_diff > 0.1 & avg_diff <= 0.3  ~ '++',
                            p_val_adj < 1e-6 & avg_diff > 0.05 & avg_diff <= 0.1  ~ '+',
                            p_val_adj < 1e-6 & avg_diff < -0.3  ~ '---',
                            p_val_adj < 1e-6 & avg_diff < -0.1 & avg_diff >= -0.3 ~ '--',
                            p_val_adj < 1e-6 & avg_diff < -0.05 & avg_diff >= -0.1  ~ '-',
                            TRUE ~ "")) %>%
  mutate(gene = rownames(.)) %>%
  dplyr::select(gene, change, group)

#'
format_change <- function(mylist, celltype.levels) {
  f.change <- lapply(mylist, format_change.df) %>%
    Reduce(rbind, .) %>%
    pivot_wider(names_from = "group", values_from = "change")
  f.change <- f.change %>% dplyr::select(c("gene", all_of(celltype.levels)))
  f.change[is.na(f.change)] <- ""
  f.change
}

