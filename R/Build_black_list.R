

# get the block gene lists from SignatuR

require(SignatuR)

# vector of list of genes to include
block.list.members <- c("cellCycle.G1S","cellCycle.G2M",
                        "Mito", "Ribo", "TCR", "Immunoglobulins",
                        "Pseudogenes", "HSP", "Non-coding"
                )

# load all signatures from SignatuR
sigs <- GetSignature(SignatuR$Hs$Programs)
sigs <- c(sigs, GetSignature(SignatuR$Hs$Blocklists))
sigs <- c(sigs, GetSignature(SignatuR$Hs$Compartments))

# filter the indicated
sigs <- sigs[names(sigs) %in% block.list.members]


# add Y-chrosome genes

# genes located in the X chromosome that have been reported to escape
# X-inactivation
# http://bioinf.wehi.edu.au/software/GenderGenes/index.html
sigs[["Xgenes"]] <- c("ARHGAP4","STS","ARSD", "ARSL", "AVPR2", "BRS3", "S100G", "CHM",
                     "CLCN4", "DDX3X","EIF1AX","EIF2S3", "GPM6B", "GRPR", "HCFC1",
                     "L1CAM", "MAOA", "MYCLP1", "NAP1L3", "GPR143", "CDK16", "PLXNB3",
                     "PRKX", "RBBP7", "RENBP", "RPS4X", "TRAPPC2", "SH3BGRL", "TBL1X",
                     "UBA1", "KDM6A", "XG", "XIST", "ZFX", "PUDP", "PNPLA4", "USP9X",
                     "KDM5C", "SMC1A", "NAA10", "OFD1", "IKBKG", "PIR", "INE2", "INE1",
                     "AP1S2", "GYG2", "MED14", "RAB9A", "ITM2A", "MORF4L2", "CA5B",
                     "SRPX2", "GEMIN8", "CTPS2", "CLTRN", "NLGN4X", "DUSP21", "ALG13",
                     "SYAP1", "SYTL4", "FUNDC1", "GAB3", "RIBC1", "FAM9C","CA5BP1")

# genes belonging to the male-specific region of chromosome Y (unique genes)
# http://bioinf.wehi.edu.au/software/GenderGenes/index.html
sigs[["Ygenes"]] <-c("AMELY", "DAZ1", "PRKY", "RBMY1A1", "RBMY1HP", "RPS4Y1", "SRY",
                  "TSPY1", "UTY", "ZFY","KDM5D", "USP9Y", "DDX3Y", "PRY", "XKRY",
                  "BPY2", "VCY", "CDY1", "EIF1AY", "TMSB4Y","CDY2A", "NLGN4Y",
                  "PCDH11Y", "HSFY1", "TGIF2LY", "TBL1Y", "RPS4Y2", "HSFY2",
                  "CDY2B", "TXLNGY","CDY1B", "DAZ3", "DAZ2", "DAZ4")
black.list <- sigs
# save object
save(black.list, file = "data/default_black_list.RData")
