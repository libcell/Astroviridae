
## (1) install R packages.

#. if (!require("BiocManager")) install.packages("BiocManager")
#. BiocManager::install("webchem")
#. BiocManager::install("ChemmineR")

## (2) get the CID of all compounds, taking convolamine as example.

library(webchem)

# help(package = "webchem")

#. drug <- c("convolamine",
#.           "esculetin",
#.           "tinidazole")

drug <- read.csv("CID for approved Scz drugs.txt",
                 header = TRUE,
                 sep = "\t")

drug <- drug$Drug.name



op <- par(mfrow = c(2, 2))

## (3) obtain the sdf files for all compounds.

library(ChemmineR)

SCZ_drug <- list()

length(SCZ_drug) <- length(drug)

names(SCZ_drug) <- drug

SCZ_did <- paste("candidate",
                 1:length(drug),
                 sep = "-")

cmpd <- new("SDFset",
            SDF = SCZ_drug,
            ID = SCZ_did)


for (d in 1:length(drug)) {
  cid <- get_cid(query = drug[d],
                 from = "name",
                 match = "first")
  cid <- cid$cid[1]

  Sys.sleep(3)

  print(cid)

  sdf_url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID",
                   cid,
                   "record/SDF/?record_type=2d&response_type=display",
                   sep = "/")

  sdfset <- read.SDFset(sdf_url)

  plot(sdfset)

  cmpd@SDF[d] <- sdfset@SDF

  cmpd@ID[d] <- drug[d]

  # cid(cmpd) <- makeUnique(cid(cmpd))

}

par(op)

plot(cmpd)


