
# Statistics information on Astroviridae.

statistics <- function(seq.type = "All") {
  # stats <- readRDS("data/Statistics.rds")
  if (!(seq.type %in% c("All", "RefSeq", "NonRef")))
    stop("Please input the right parameter seq.type, and it must be one of 'All', 'RefSeq' or 'NonRef'!")
  if (seq.type == "All")
    content <- stats
  if (seq.type == "RefSeq")
    content <- stats[grep("NCBI RefSeq", stats$`Annotation Name`), ]
  if (seq.type == "NonRef")
    content <- stats[-grep("NCBI RefSeq", stats$`Annotation Name`), ]
  rownames(content) <- 1:nrow(content)
  print(content)
}
