
################################################################################
#    &&&....&&&    % Project: Update data sets used in Astroviridae package    #
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
### Step-01. Install all R packages in this project.

#. pkg_cran <- c("BiocManager", "readr", "ape")
#. pkg_bioc <- c("BioStrings", "msa", "ggmsa", "adegenet", "bios2mds", )

if (!require("BiocManager"))
  install.packages("BiocManager")
if (!require("readr"))
  install.packages("readr")
if (!require("ape"))
  install.packages("ape")

if (!require("Biostrings"))
  BiocManager::install("Biostrings")
if (!require("msa"))
  BiocManager::install("msa")

if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
if (!require("ggmsa"))
  devtools::install_github("YuLab-SMU/ggmsa")

if (!require("adegenet"))
  BiocManager::install("adegenet")
if (!require("bios2mds"))
  BiocManager::install("bios2mds")

### ------------------------------------------------------------------------ ###
### Step-02. Download the newest version statistics file from NCBI website.

# Enter https://www.ncbi.nlm.nih.gov/genome/?term=Astroviridae
# Click "Genomes", and get into the following web link.
# https://.../datasets/genomes/?taxon=39733&utm_source=genome&utm_medium=referral
# save as "Statistics.tsv"
# And meanwhile, downloading the newest version of ncbi_dataset.zip


### ------------------------------------------------------------------------ ###
### Step-03. Set up the working directory.

wkdir <- getwd()
user <- strsplit(x = wkdir,
                 split = "/")[[1]][3]


### ------------------------------------------------------------------------ ###
### Step-04. Firstly, update the object - Statistics.rds

library(readr)
stats_path <- file.path("C:",
                        "Users",
                        user,
                        "Downloads",
                        "Statistics.tsv")

stats <- read_tsv(stats_path)
stats <- as.data.frame(stats)

seq.refseq <- stats[grep("NCBI RefSeq", stats$`Annotation Name`), ]
seq.speseq <- stats[-grep("NCBI RefSeq", stats$`Annotation Name`), ]

### ------------------------------------------------------------------------ ###
### Step-05. Exploring the information for Astroviridae.

# - genome size distribution.
asv.len <- stats$`Assembly Stats Total Sequence Length`

hist(asv.len,
     breaks = 10,
     col = terrain.colors(15),
     main = "Distribution of Genome Size for Astroviridae")

# - The timeline on discovery of viruses from Astroviridae.
time.line <- stats$`Assembly Submission Date`
time.line <- sapply(as.character(time.line),
                    function(x) strsplit(x, "-")[[1]][1])
time.line <- as.numeric(time.line)

stats$Year <- time.line # added the year into the stats object.

tl.x <- seq(min(time.line),
            max(time.line),
            by = 1)
tl.y <- rep(NA,
            times = length(tl.x))

for (i in 1:length(tl.y)) {
  count <- length(which(as.character(as.character(time.line)) == as.character(tl.x[i])))
  tl.y[i] <- count
}

plot(tl.x,
     tl.y,
     xlab = "Year",
     ylab = "Genome Count, in Astroviridae",
     type = "h",
     col = terrain.colors(length(tl.x)),
     lwd = 4)

pie(sort(table(time.line),
         decreasing = TRUE),
    init.angle = 0,
    clockwise = TRUE,
    srt = 0)


### ------------------------------------------------------------------------ ###
### Step-06. Firstly, update the genome sequence files.

library(utils)

if (!dir.exists("genome_data"))
  dir.create("genome_data")

# For windows, work directory.
setwd("genome_data")

if (!dir.exists("ncbi_dataset")) {

  dtst_path <- file.path("C:",
                         "Users",
                         user,
                         "Downloads",
                         "ncbi_dataset.zip")


  file.copy(dtst_path, ".")
  unzip("ncbi_dataset.zip")
}


### ------------------------------------------------------------------------ ###
### Step-07. Re-update the statistics.rds after obtaning the names for all viruses.
#
# - 1) Read the genome sequences of all viruses, and generate a list object in R.

library(Biostrings)

AsV.id <- stats$`Assembly Accession`

dna.list <- list()

p <- 0

for (i in AsV.id) {
  p <- p + 1
  file.path <- paste("./ncbi_dataset/data",
                     i,
                     sep = "/")
  setwd(file.path)
  fas <- readDNAStringSet(dir()[1])
  setwd(".."); setwd(".."); setwd("..")
  dna.list[[p]] <- fas
}

all.seq <- do.call(c, dna.list)

# - 2) Add two new features on all genome sequences in statistics.rds.

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

# DT::datatable(stats)

save(stats,
     file = "Statistics.RData")

file.copy("Statistics.RData",
          "D:/00-GitHub/Astroviridae/data",
          overwrite = TRUE)


### ------------------------------------------------------------------------ ###
### Step-08. Save the all genome sequences as one file in fasta format.

writeXStringSet(all.seq,
                "Astroviridae_genomes.fas",
                append = FALSE,
                compress = FALSE,
                compression_level = NA,
                format = "fasta")

file.copy("Astroviridae_genomes.fas",
          "D:/00-GitHub/Astroviridae/inst/extdata/",
          overwrite = TRUE)


### ------------------------------------------------------------------------ ###
### Step-09. Choose the viruses for constructing the dendrogram.
# Update the names of all viruses, i. e., the future labels of evolutionary tree.

# - 1) Re-name all sequences with their identifiers (or Organism Name).

# names(all.seq) <- stats$`Organism Name`
names(all.seq) <- stats$Identifier

# - 2) Pick sequences.

select.seq <- all.seq[stats$is.FullGnm == "Yes" & stats$is.RefSeq == "Yes"]

# rownames(select.anno) <- 1:nrow(select.anno)

# - 3) Implement Multiple sequence alignment with msa package.

library(msa)
align.seq <- msa(select.seq,
                 method = "ClustalW")
print(align.seq)
print(align.seq, show = "complete")


### ------------------------------------------------------------------------ ###
### Step-10. Extract the specific reqion, from sequences aligned or not aligned.

subseq(all.seq, start = 5000, end  = 5100)
subseq(all.seq[[1]], start = 1, end = 50)

consv.motif <- msaConsensusSequence(align.seq)

# saveRDS(all.seq, file = "AsvDB.rds")


### ------------------------------------------------------------------------ ###
### Step-11. Convert the aligned genome to fasta format and save it as a file.

library(bios2mds)

DNA <- msaConvert(align.seq,
                  type = "bios2mds::align")
export.fasta(DNA,
             outfile = "aligned_genomes.fas",
             ncol = 60, open = "w")

file.copy("aligned_genomes.fas",
          "D:/00-GitHub/Astroviridae/inst/extdata/",
          overwrite = TRUE)


### ------------------------------------------------------------------------ ###
### Step-12. Read the aligned genome sequences and construct the evolution tree.

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


out.id <- 1:2
tre2 <- root(tre, out = out.id)
tre2 <- ladderize(tre2)
plot(tre2,
     show.tip = TRUE,
     edge.width = 2)
tre.title <- paste("Rooted NJ tree",
                   ", outgroup =",
                   tre$tip.label[out.id],
                   sep = " ")


### ------------------------------------------------------------------------ ###
### Step-13. Read genome sequences from Astroviridae, to illustrate visualization.

library(ggmsa)

all.fas <- system.file("extdata",
                         "aligned_genomes.fas",
                         package = "Astroviridae")

ggmsa(all.fas,
      start = 2221,
      end = 2280,
      char_width = 0.5,
      seq_name = T) +
  geom_seqlogo() +
  geom_msaBar()

### ------------------------------------------------------------------------ ###
### Step-14. Reture the primary working directory.

setwd(wkdir)


### ************************************************************************ ###
### End of Here!
### ************************************************************************ ###
# method = average is used for UPGMA, members can be equal to NULL or a vector
# with a length of size D

D[is.na(D)] <- 3

h_cluster <- hclust(D)

plot(h_cluster, cex = 0.6)



data(bird.orders)
plot(root(bird.orders, 1))
plot(root(bird.orders, 2))
plot(root(bird.orders, 3))
plot(root(bird.orders, 4))
plot(root(bird.orders, 5))
plot(root(bird.orders, 1:5))

tr <- root(bird.orders, 1)
is.rooted(bird.orders) # yes
is.rooted(tr)          # no
### This is because the tree has been unrooted first before rerooting.
### You can delete the outgroup...
is.rooted(drop.tip(tr, "Struthioniformes"))
### ... or resolve the basal trichotomy in two ways:
is.rooted(multi2di(tr))
is.rooted(root(bird.orders, 1, r = TRUE))
### To keep the basal trichotomy but forcing the tree as rooted:
tr$root.edge <- 0
is.rooted(tr)

x <- setNames(rmtree(10, 10), LETTERS[1:10])
is.rooted(x)




