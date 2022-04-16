################################################################################
#    &&&....&&&    % Project: Data analysis and mining of student score        #
#  &&&&&&..&&&&&&  % Author: Bo Li, Xiaoxi Zhang                               #
#  &&&&&&&&&&&&&&  % Date: Apr. 15th, 2022                                     #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 4.1.1;                             #
#       &&&&       % x86_64-w64-mingw32/x64 (64-bit)bit)                       #
#        &         %                                                           #
################################################################################

### ****************************************************************************
### code chunk number 01: Run this R code to update the database information.
### ****************************************************************************

### ------------------------------------------------------------------------ ###
### Step-01. Reading the score stored in *.xlsx file in R environment

## (1) Download the newest version statistics file from NCBI website. 
# https://.../datasets/genomes/?taxon=39733&utm_source=genome&utm_medium=referral
# save as "Statistics.tsv"
# And meanwhile, downloading the newest version of ncbi_dataset.zip

## (2) Set up the working directory. 
user <- "admin"
wkdir <- getwd()

## (3) Firstly, update the object - Statistics.rds
library(readr)
stats_path <- paste("C:/Users", user, "Downloads/Statistics.tsv", sep = "/")
stats <- read_tsv(stats_path)
stats <- as.data.frame(stats)

# library(DT)
# datatable(stats)

seq.refseq <- stats[grep("NCBI RefSeq", stats$`Annotation Name`), ]
seq.speseq <- stats[-grep("NCBI RefSeq", stats$`Annotation Name`), ]

# saveRDS(stats, file = "Statistics.rds")

## (4) Exploring the information for Astroviridae. 

# - genome size distribution. 

asv.len <- stats$`Assembly Stats Total Sequence Length`

hist(asv.len, breaks = 10, 
     main = "Distribution of Genome Size for Astroviridae")

# - The timeline on discovery of viruses from Astroviridae. 

time.line <- stats$`Assembly Submission Date`

time.line <- sapply(as.character(time.line), function(x) strsplit(x, "-")[[1]][1])

time.line <- as.numeric(time.line)

tl.x <- seq(min(time.line), max(time.line), by = 1)

tl.y <- rep(NA, times = length(tl.x))

for (i in 1:length(tl.y)) {
  
  count <- length(which(as.character(as.character(time.line)) == as.character(tl.x[i])))
  
  tl.y[i] <- count
  
}

plot(tl.x, tl.y, type = "h", lwd = 4)

pie(sort(table(time.line), decreasing = TRUE), init.angle = 90, clockwise = TRUE)

## (5) Firstly, update the genome sequence files. 

library(utils)

dir.create("genome_data")

# For windows, work directory.
setwd("genome_data")

file.copy("C:/Users/libo//Downloads/ncbi_dataset.zip", ".")

unzip("ncbi_dataset.zip")

# help(package = "utils")


# BiocManager::install("Biostrings")
# BiocManager::install("msa")

# Read the genome sequences of all viruses. 

library(Biostrings)
# library(msa)

AsV.id <- stats$`Assembly Accession`

dna.list <- list()

p <- 0

for (i in AsV.id) {
  p <- p + 1
  file.path <- paste("./ncbi_dataset/data", i, sep = "/")
  setwd(file.path)
  fas <- readDNAStringSet(dir()[1])
  setwd(".."); setwd(".."); setwd("..")
  dna.list[[p]] <- fas
}

all.seq <- do.call(c, dna.list)

all.seq
saveRDS(all.seq, file = "AsvDB.rds")


# Reupdate the statistics.rds after obtaning the names for all viruses. 

AsV.anno <- names(all.seq)

full.gnm <- rep("No", length(AsV.anno))
full.gnm[grep("complete genome", AsV.anno)] <- "Yes"
full.gnm[grep("complete sequence", AsV.anno)] <- "Yes"

ref.seq <- rep("No", length(AsV.anno))
ref.seq[grep("NCBI RefSeq", stats$`Annotation Name`)] <- "Yes"

stats$is.FullGnm <- full.gnm
stats$is.RefSeq <- ref.seq

splt <- function(x) {
  strsplit(x, " ")[[1]][1]
}

res <- unlist(sapply(AsV.anno, splt))

names(res) <- NULL

stats$Identifier <- res

saveRDS(stats, file = "Statistics.rds")


# (3) Choose the viruses for constructing the dendrogram. 
# Update the names of all viruses, i. e., the future labels of evolutionary tree. 

names(all.seq) <- stats$`Organism Name`

select.seq <- all.seq[stats$is.FullGnm == "Yes" & stats$is.RefSeq == "Yes"]

library(msa)

align.seq <- msa(select.seq, method = "ClustalW") # Muscle

align.seq

#. print(align.seq, show = "complete")
#.
msaPrettyPrint(
  align.seq,
  output = "pdf",
  showNames = "none",
  showLogo = "none",
  askForOverwrite = FALSE,
  verbose = FALSE
)

 msaPrettyPrint(
   align.seq,
   y = c(164, 213),
   output = "asis",
   showNames = "none",
   showLogo = "none",
   askForOverwrite = FALSE
 )


library(bios2mds)

DNA <- msaConvert(align.seq, type = "bios2mds::align")

export.fasta(DNA, outfile = "test_alignment.fa", ncol = 60, open = "w")

library(adegenet)
dna <- fasta2DNAbin(file="test_alignment.fa")

dna

library(ape)
D <- dist.dna(dna, model = "TN93")
length(D) #number of pairwise distances, computed as n(n-1)/2
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg = , clabel.row=.5, clabel.col=.5) #
tre <- njs(D)
class(tre) #all trees created using {ape} package will be of class phylo
tre <- ladderize(tre)
tre # te
plot(tre, cex = 0.6)
title("A Simple NJ Tree")

myBoots <- boot.phylo(tre, dna, function(e) root(nj(dist.dna(e, model = "TN93")),1))

plot(tre, show.tip=FALSE, edge.width=2)
title("NJ tree + bootstrap values")
tiplabels(frame="none", pch=20, col=transp(num2col(annot$year, col.pal=myPal),.7), cex=3, fg="transparent")

axisPhylo()
temp <- pretty(1993:2008, 5)
legend("topright", fill=transp(num2col(temp, col.pal=myPal),.7), leg=temp, ncol=2)
nodelabels(myBoots, cex=.6)


library(ggmsa)
protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
ggmsa(protein_sequences, start = 221, end = 280, char_width = 0.5, seq_name = T) + geom_seqlogo() + geom_msaBar()

all.fas <- msaConvert(align.seq, type = "ape::DNAbin")

ggmsa(all.fas, start = 221, end = 280, char_width = 0.5, seq_name = T) + 
  geom_seqlogo() + 
  geom_msaBar()










