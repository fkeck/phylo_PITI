library(ape)
library(adephylo)
library(stringr)
library(phylobase)
library(phylosignal)
library(Rphylopars)
library(diatobc)

# User paths to be edited
# setwd("IF_NECESSARY")
raxml.exe <- "PATH_TO_RAXML"
pathd8.exe <- "PATH_TO_PATHd8"
muscle.exe <- "PATH_TO_MUSCLE"

source("R/utils.R")

rbcl.tree <- read.tree("data/rbcl_ref604_bestML.phy")
rbcl.seq <- read.dna("data/rbcl_ref604_align.fas", format = "fasta")
rbcl.bc <- rbcl.seq[ , 672:983]

trait <- read.csv("data/IPS_values.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(trait) <- trait$species

rbcl.tree.cross <- rbcl.tree
rbcl.tree.cross <- root.phylo(rbcl.tree.cross, outgroup = "Bolidomonas_pacifica", resolve.root = TRUE)
rbcl.tree.cross$tip.label <- gsub("_", " ", rbcl.tree.cross$tip.label)
rbcl.tree.cross <- drop.tip(rbcl.tree.cross, tip = setdiff(rbcl.tree.cross$tip.label, trait$species))
rbcl.tree.cross$node.label <- NULL
rbcl.tree.cross <- pathd8(rbcl.tree.cross, pathd8.exe)
rbcl.tree.cross$tip.label <- gsub("_", " ", rbcl.tree.cross$tip.label)

trait.p4d <- phylo4d(rbcl.tree.cross, tip.data = trait[ , c("IPSS", "IPSV")])

OTU.list <- read.csv("data/OTU_list.csv", header = TRUE,
                     stringsAsFactors = FALSE)
OTU.list.sites <- read.csv("data/OTU_list_sites.csv", header = TRUE,
                           row.names = 1, stringsAsFactors = FALSE)
OTU.list.sites <- t(OTU.list.sites)
OTU.seq <- read.dna("data/OTU_seq.fas", format = "fasta")
IPS.classfied.seq <- read.csv("data/IPS_classified_seq.csv", header = TRUE)

sites.IPS.scores <- read.csv("data/sites_IPS_scores.csv", header = TRUE,
                             row.names = 1, stringsAsFactors = FALSE)


##### BARCODES REPLACEMENT #####

# Barcode insertion in trees
trees.inserted.bc <- list()
for(i in rbcl.tree$tip.label){
  rbcl.tree.no.i <- drop.tip(rbcl.tree, i)
  rbcl.seq.no.i <- rbcl.seq[rownames(rbcl.seq) != i, ]
  rbcl.bc.i <- rbcl.bc[i, ]
  
  rbcl.seq.bc.i.insert <- insertBC(bc = rbcl.bc.i, ref = rbcl.seq.no.i,
                                   region = c(672, 983), align.method = "none")
  
  rbcl.tree.bc.i.insert <- raxmlEPA(rbcl.seq.bc.i.insert, ref.tree = rbcl.tree.no.i,
                                    model = "GTRGAMMAI", nruns = 20, threads = 4,
                                    raxml.exec = raxml.exe)
  
  trees.inserted.bc[[i]] <- rbcl.tree.bc.i.insert$tree
  print(match(i, rbcl.tree$tip.label))
}


trees.inserted.bc <- lapply(trees.inserted.bc,
                            function(x) root.phylo(x, outgroup = "Bolidomonas_pacifica",
                                                   resolve.root = TRUE))

taxa.list <- setdiff(names(trees.inserted.bc), "Bolidomonas_pacifica")
displace.patristic <- vector(mode = "numeric", length = 0)
displace.node <- vector(mode = "numeric", length = 0)

for (i in taxa.list){
  taxa.chk <- i
  tr1 <- rbcl.tree
  tr2 <- trees.inserted.bc[[taxa.chk]]
  
  tr1.mrca <- tr1$edge[tr1$edge[,2] == which(tr1$tip.label == taxa.chk), 1]
  tr1.mrca.dsc <- tr1$tip.label[setdiff(unlist(phangorn::Descendants(tr1, tr1.mrca, "tips")),
                                        which(tr1$tip.label == taxa.chk))]
  if(length(tr1.mrca.dsc) == 1){
    tr2.refmrca <- tr2$edge[tr2$edge[,2] == which(tr2$tip.label == tr1.mrca.dsc), 1]
  } else {
    tr2.refmrca <- getMRCA(tr2, tr1.mrca.dsc)
  }
  tr2.mrca <- tr2$edge[tr2$edge[,2] == which(tr2$tip.label == taxa.chk), 1]
  
  tr2.dnodes <- dist.nodes(tr2)
  displace.patristic[i] <- tr2.dnodes[tr2.refmrca, tr2.mrca]
  
  tr2$edge.length <- rep(1, length(tr2$edge.length))
  tr2.dnodes <- dist.nodes(tr2)
  displace.node[i] <- tr2.dnodes[tr2.refmrca, tr2.mrca]
  
  print(match(i, names(trees.inserted.bc)))
}


pdf(file = "results/Figure_1.pdf", width = 7, height = 7)
plot(table(displace.node), main = "Insertion difference",
     xlab = "Node distance from correct insertion position",
     ylab = "Count", lwd = 5, lend = 1, las = 1)
dev.off()



##### BARCODES ESTIMATION #####

IPS.phylosignal <- phyloSignal(trait.p4d, "Lambda")
IPS.phylosignal <- as.data.frame(IPS.phylosignal)
colnames(IPS.phylosignal) <-  c("stat", "pvalue")
write.table(IPS.phylosignal, "results/phylo_signal_IPS.txt")

mat <- matrix(NA, nrow = nrow(trait), ncol = 2,
              dimnames = list(trait$species, c("IPSS", "IPSV")))
IPS.estim <- list(BM = mat, OU = mat, lambda = mat,
                  kappa = mat, delta = mat, EB = mat, star = mat)
IPS.var <- list(BM = mat, OU = mat, lambda = mat,
                kappa = mat, delta = mat, EB = mat, star = mat)

for (i in trait$species){
  
  tree <- trees.inserted.bc[[gsub(" ", "_", i)]]
  tree <- pathd8(tree, pathd8.exe)
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  tree <- drop.tip(tree, tip = setdiff(tree$tip.label, trait$species))
  trait.NA <- trait[match(tree$tip.label, trait$species), ]
  trait.NA[trait.NA$species == i, c("IPSS", "IPSV")] <- NA
  
  pp <- phylopars(trait_data = trait.NA, tree = tree, model = "BM")
  IPS.estim$BM[i, ] <- pp$anc_recon[i, ]
  IPS.var$BM[i, ] <- pp$anc_var[i, ]
  
  pp <- phylopars(trait_data = trait.NA, tree = force.ultrametric(tree), model = "OU")
  IPS.estim$OU[i, ] <- pp$anc_recon[i, ]
  IPS.var$OU[i, ] <- pp$anc_var[i, ]
  
  pp <- phylopars(trait_data = trait.NA, tree = tree, model = "lambda")
  IPS.estim$lambda[i, ] <- pp$anc_recon[i, ]
  IPS.var$lambda[i, ] <- pp$anc_var[i, ]
  
  pp <- phylopars(trait_data = trait.NA, tree = tree, model = "kappa")
  IPS.estim$kappa[i, ] <- pp$anc_recon[i, ]
  IPS.var$kappa[i, ] <- pp$anc_var[i, ]
  
  pp <- phylopars(trait_data = trait.NA, tree = tree, model = "delta")
  IPS.estim$delta[i, ] <- pp$anc_recon[i, ]
  IPS.var$delta[i, ] <- pp$anc_var[i, ]
  
  pp <- phylopars(trait_data = trait.NA, tree = tree, model = "EB")
  IPS.estim$EB[i, ] <- pp$anc_recon[i, ]
  IPS.var$EB[i, ] <- pp$anc_var[i, ]
  
  pp <- phylopars(trait_data = trait.NA, tree = tree, model = "star")
  IPS.estim$star[i, ] <- pp$anc_recon[i, ]
  IPS.var$star[i, ] <- pp$anc_var[i, ]
  
  print(match(i, trait$species))
}


MSE.IPSS <- sort(sapply(IPS.estim, function(x) sum((trait$IPSS - x[, "IPSS"])^2)/length(trait$IPSS)))
MSE.IPSV <- sort(sapply(IPS.estim, function(x) sum((trait$IPSV - x[, "IPSV"])^2)/length(trait$IPSV)))

SE.IPSS <- sapply(IPS.estim, function(x) (x[, "IPSS"] - trait$IPSS)^2)
SE.IPSV <- sapply(IPS.estim, function(x) (x[, "IPSV"] - trait$IPSV)^2)

SE.IPSS.col <- do.call("rbind", lapply(as.list(as.data.frame(SE.IPSS)), as.matrix))
SE.IPSV.col <- do.call("rbind", lapply(as.list(as.data.frame(SE.IPSV)), as.matrix))
segrp <- rep(colnames(SE.IPSS), each = nrow(SE.IPSS))

wilcox.IPSS.pw.tests <- pairwise.wilcox.test(SE.IPSS.col, segrp, paired = TRUE,
                                             p.adjust.method = "bonferroni")
wilcox.IPSV.pw.tests <- pairwise.wilcox.test(SE.IPSV.col, segrp, paired = TRUE,
                                             p.adjust.method = "bonferroni")

pdf(file = "results/Figure_2.pdf", width = 6, height = 4)
par(mfrow = c(1, 2))
barplot(MSE.IPSS, las = 2, xlab = "Model of evolution", ylab = "Mean squared error", main = "IPSS")
barplot(MSE.IPSV, las = 2, xlab = "Model of evolution", ylab = "Mean squared error", main = "IPSV")
dev.off()


IPSS.e <- matrix(data = IPS.estim[["lambda"]][, "IPSS"] - trait[, "IPSS"], ncol = 1)
rownames(IPSS.e) <- names(IPS.estim[["lambda"]][, "IPSS"])
colnames(IPSS.e) <- "IPSS"

IPSS.col <- IPSS.e
IPSS.col[abs(IPSS.e[, 1]) <= 1, 1] <- "darkgreen"
IPSS.col[abs(IPSS.e[, 1]) > 1 & abs(IPSS.e[, 1]) <= 2, 1] <- "orange"
IPSS.col[abs(IPSS.e[, 1]) > 2, 1] <- "red"

pdf(file = "results/Figure_3.pdf", width = 10, height = 10)
dotplot.phylo4d(trait.p4d, trait = "IPSS", tree.ladderize = TRUE, tree.type = "fan",
                error.bar.sup = IPSS.e, error.bar.col = IPSS.col,
                center = FALSE, scale = FALSE, grid.col = "white", data.xlim = c(0, 5),
                grid.vertical = TRUE, tree.open.angle = 5, 
                tip.cex = 0.5, dot.cex = 0.5, dot.col = IPSS.col, show.trait = FALSE,
                show.data.axis = TRUE, grid.lty = 1)
dev.off()


##### NATURAL SAMPLES OTU ESTIMATION #####

rbcl.seq.OTU.ins <- insertBC(bc = OTU.seq, ref = rbcl.seq,
                             region = c(672, 983),
                             align.method = "muscle",
                             align.exec = muscle.exe)

rbcl.tree.OTU.ins <- raxmlEPA(rbcl.seq.OTU.ins, ref.tree = rbcl.tree,
                              model = "GTRGAMMAI",
                              nruns = 1, threads = 4,
                              raxml.exec = raxml.exe,
                              extra.com = "--epa-keep-placements=10 --epa-prob-threshold=0.1")

rbcl.tree.OTU.ins <- rbcl.tree.OTU.ins$tree
rbcl.tree.OTU.ins <- pathd8(rbcl.tree.OTU.ins, pathd8.exe)
rbcl.tree.OTU.ins$tip.label <- gsub("_", " ", rbcl.tree.OTU.ins$tip.label)

pp.IPSS <- phylopars(trait_data = trait[, c(1, 2)],
                     tree = rbcl.tree.OTU.ins,
                     model = "lambda")
pp.IPSV <- phylopars(trait_data = trait[, c(1, 3)],
                     tree = force.ultrametric(rbcl.tree.OTU.ins),
                     model = "OU")
IPS.phylo <- cbind(pp.IPSS$anc_recon, pp.IPSV$anc_recon)


IPS.phylomoth <- IPS.phylo
IPS.phylomoth[OTU.list$OTU_Label[OTU.list$Species != "unclassified"], "IPSS"] <- IPS.classfied.seq$IPSS
IPS.phylomoth[OTU.list$OTU_Label[OTU.list$Species != "unclassified"], "IPSV"] <- IPS.classfied.seq$IPSV
IPS.phylomoth[is.na(IPS.phylomoth[, "IPSS"]), "IPSS"] <- IPS.phylo[is.na(IPS.phylomoth[, "IPSS"]), "IPSS"]
IPS.phylomoth[is.na(IPS.phylomoth[, "IPSV"]), "IPSV"] <- IPS.phylo[is.na(IPS.phylomoth[, "IPSV"]), "IPSV"]


idx.phylo <- apply(OTU.list.sites, 1, function(x) scoreIPSCustom(x, IPS.custom = IPS.phylo, IPS20 = F))
sites.IPS.scores$IPS.DNAPHYLO <- idx.phylo[rownames(sites.IPS.scores)]

idx.phylomoth <- apply(OTU.list.sites, 1, function(x) scoreIPSCustom(x, IPS.custom = IPS.phylomoth, IPS20 = F))
sites.IPS.scores$IPS.DNAHYBRID <- idx.phylomoth[rownames(sites.IPS.scores)]


pdf(file = "results/Figure_4.pdf", width = 6, height = 6)
pairs(sites.IPS.scores,
      xlim = c(1, 5), ylim = c(1, 5),
      labels = c("IPS-MICROTAXO", "IPS-DNATAXO", "IPS-DNAPHYLO", "IPS-DNAHYBRID"),
      lower.panel = panel.cor,
      upper.panel = panel.txt,
      diag.panel = panel.hist
)
dev.off()


SE.IPS <- apply(sites.IPS.scores[, -1], 2, function(x) (x - sites.IPS.scores[, 1])^2)
SE.IPS.col <- do.call("rbind", lapply(as.list(as.data.frame(SE.IPS)), as.matrix))
segrp.IPS <- rep(colnames(SE.IPS), each = nrow(SE.IPS))

wilcox.IPS.pw.tests <- pairwise.wilcox.test(SE.IPS.col, segrp.IPS, paired = TRUE,
                                            p.adjust.method = "bonferroni")


save.image("results/workspace.RData")
