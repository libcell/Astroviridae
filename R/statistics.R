
# Statistics information on Astroviridae.

statistics <- function(seq.type = "all") {
  # stats <- readRDS("data/Statistics.rds")
  if (seq.type == "all")
    content <- stats
  if (seq.type == "RefSeq")
    content <- stats[grep("NCBI RefSeq", stats$`Annotation Name`), ]
  if (seq.type == "others")
    content <- stats[-grep("NCBI RefSeq", stats$`Annotation Name`), ]
  print(content)
}
