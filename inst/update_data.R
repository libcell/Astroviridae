
################################################################################
#    &&&....&&&    % Project: Data analysis and mining of student score        #
#  &&&&&&..&&&&&&  % Author: Bo Li, Xiaoxi Zhang, Xiner Nie                    #
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
user <- "libo"
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

stats$Year <- time.line # added the year into the stats object. 

tl.x <- seq(min(time.line), max(time.line), by = 1)
tl.y <- rep(NA, times = length(tl.x))

for (i in 1:length(tl.y)) {
  count <- length(which(as.character(as.character(time.line)) == as.character(tl.x[i])))
  tl.y[i] <- count
}

plot(tl.x, tl.y, type = "h", lwd = 4)
pie(sort(table(time.line), decreasing = TRUE), 
    init.angle = 0, 
    clockwise = TRUE, 
    srt = 0)

## (5) Firstly, update the genome sequence files. 

library(utils)

if (!dir.exists("genome_data")) 
  dir.create("genome_data")

# For windows, work directory.
setwd("genome_data")

if (!dir.exists("ncbi_dataset")) {
  file.copy("C:/Users/libo//Downloads/ncbi_dataset.zip", ".")
  unzip("ncbi_dataset.zip")
}

# BiocManager::install("Biostrings")
# BiocManager::install("msa")

## (6) Re-update the statistics.rds after obtaning the names for all viruses. 
# 
# - Read the genome sequences of all viruses, and generate a list object in R. 

library(Biostrings)

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

# - Add two new features on all genome sequences in statistics.rds. 

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

DT::datatable(stats)

# saveRDS(stats, file = "Statistics.rds")

save(stats, file = "Statistics.RData")

file.copy("Statistics.RData", "D:/00-GitHub/Astroviridae/data")

## (7) Save the all genome sequences as one file in fasta format. 

writeXStringSet(all.seq, 
                "Astroviridae_genome.fas", 
                append = FALSE,
                compress = FALSE, 
                compression_level = NA, 
                format = "fasta")

file.copy("Astroviridae_genome.fas", "D:/00-GitHub/Astroviridae/inst/extdata/")

# (8) Choose the viruses for constructing the dendrogram. 
# Update the names of all viruses, i. e., the future labels of evolutionary tree. 

# names(all.seq) <- stats$`Organism Name`
names(all.seq) <- stats$Identifier

### Pick sequences. 

select.seq <- all.seq[stats$is.FullGnm == "Yes" & stats$is.RefSeq == "Yes"]
# select.anno <- stats[stats$is.FullGnm == "Yes" & stats$is.RefSeq == "Yes", ]
# rownames(select.anno) <- 1:nrow(select.anno)

library(msa)

align.seq <- msa(select.seq, 
                 method = "ClustalW")

align.seq

#. print(align.seq, show = "complete")

## (9) Extract the specific reqion, from sequences aligned or not aligned. 

subseq(all.seq, start = 5000, end  = 5100)
subseq(all.seq[[1]], start = 1, end = 50)

msaConsensusSequence(align.seq)
# saveRDS(all.seq, file = "AsvDB.rds")

## (10) Convert the aligned genome to fasta format and save it as a file. 

library(bios2mds)

DNA <- msaConvert(align.seq, 
                  type = "bios2mds::align")
export.fasta(DNA, 
             outfile = "aligned_genomes.fas", 
             ncol = 60, open = "w")

file.copy("aligned_genomes.fas", "D:/00-GitHub/Astroviridae/inst/extdata/")

## (11) Read the aligned genome sequences and construct the evolution tree. 

library(adegenet)
dna <- fasta2DNAbin(file = "aligned_genomes.fas")
anno <- stats[match(rownames(dna), stats$Identifier), ]

library(ape)
D <- dist.dna(dna, model = "TN93")
length(D) #number of pairwise distances, computed as n(n-1)/2
temp <- as.data.frame(as.matrix(D))
DT::datatable(temp)
temp[is.na(temp)] <- max(temp, na.rm = TRUE) + 1
table.paint(temp, cleg = 0, clabel.row = .5, clabel.col = .5)
abline(v = 4.5, lwd = 3, col = "red")

tre <- njs(D)
class(tre)
# Reorganizes internal structure of tree to get ladderized effect when plotted.
tre <- ladderize(tre) 
tre
plot(tre, cex = 0.6)
title("A Simple NJ Tree")

# Construct the dendrogram based on bootstrapping technique. 

myBoots <- boot.phylo(tre, 
                      dna, 
                      function(e) root(njs(dist.dna(e, model = "TN93")), 1))

plot(tre, show.tip = TRUE, edge.width = 2)
title("NJ tree + bootstrap values")

myPal <- colorRampPalette(c("red","yellow","green","blue"))

tiplabels(frame = "none", pch = 20, cex = 3, 
          col = transp(num2col(anno$Year, col.pal = myPal), 1.0), 
          fg = "transparent")

axisPhylo()
temp <- pretty(min(anno$Year):max(anno$Year), 10)
legend("topright", 
       fill = transp(num2col(temp, col.pal = myPal), .7), 
       leg = temp, 
       ncol = 2)

nodelabels(myBoots, 
           bg = NULL, 
           frame = "none", 
           col = "red", 
           cex = .8)

## 
# Taking protein sequences from ggmsa as example, to illustrate visualization.

library(ggmsa)

all.fas <- system.file("inst/extdata", 
                         "aligned_genomes.fas", 
                         package = "Astroviridae")

all.fas <- "C:/Users/libo/Documents/genome_data/aligned_genomes.fas"

ggmsa(all.fas, 
      start = 2221, 
      end = 2280, 
      char_width = 0.5, 
      seq_name = T) + 
  geom_seqlogo() + 
  geom_msaBar()


