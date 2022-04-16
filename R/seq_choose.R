
# Some useful keyboard shortcuts for package authoring:
# genome_sequences <- readRDS("inst/extdata/AsVDB.rds", package = "Astroviridae")
# all.seq <- readRDS(genome_sequences)

# AsV.lnm <- names(all.seq)

seq.pick <- function(FullGnm = TRUE, RefSeq = TRUE) {

  if (FullGnm) fg <- "Yes" else fg <- "No"
  if (RefSeq) rs <- "Yes" else rs <- "No"

  select.seq <- all.seq[stats$is.FullGnm == fg & stats$is.RefSeq == rs]
  return(select.seq)
}





