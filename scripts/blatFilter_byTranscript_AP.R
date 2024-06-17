#############
## chimp ##
#############
source("~/palmer_scratch/XSAnno/bin/Functions_BlatFilter_byTranscript.r")
setwd("/vast/palmer/scratch/noonan/ap2549/XSAnno/BLAT_testing2/blat")

# Updated file names to match the provided data files
sp12sp1 = "blat.exon.hg38Tohg38.filtered.txt"
sp22sp2 = "blat.exon.panTro6TopanTro6.filtered.txt"
sp12sp2 = "blat.exon.hg38TopanTro6.filtered.txt"
sp22sp1 = "blat.exon.panTro6Tohg38.filtered.txt"

sp1 = "hg38"
sp2 = "panTro6"

blat.Sp1ToSp1 <- read.table(sp12sp1, as.is = T, header = T)
blat.Sp1ToSp2 <- read.table(sp12sp2, as.is = T, header = T)
blat.Sp2ToSp1 <- read.table(sp22sp1, as.is = T, header = T)
blat.Sp2ToSp2 <- read.table(sp22sp2, as.is = T, header = T)

# the vectors of PID and PL where the threshold to choose from
IDs <- seq(0.8, 0.999, 0.01)
PLs <- seq(0.8, 0.999, 0.05)

chooseThreshold.inter <- function(blat.Sp1ToSp2, blat.Sp2ToSp1, IDs, PLs) {
    orig.sp1.0 <- do.call(rbind, strsplit(blat.Sp1ToSp2[,1], split = "\\|"))
    orig.sp2.0 <- do.call(rbind, strsplit(blat.Sp2ToSp1[,1], split = "\\|"))

    orig.sp1 <- cbind(apply(orig.sp1.0[,1:5], 1, paste, collapse = "|"), orig.sp1.0[,6:ncol(orig.sp1.0)])
    orig.sp2 <- cbind(apply(orig.sp2.0[,1:5], 1, paste, collapse = "|"), orig.sp2.0[,6:ncol(orig.sp2.0)])

    blat.Sp1ToSp2[,1] <- orig.sp1[,1]
    blat.Sp2ToSp1[,1] <- orig.sp2[,1]
    ########################
    ## filter interSpecies
    # remove low ID low PL and duplicate regions
    exonNums <- sapply(PLs, function(interPL) {
        sapply(IDs, function(interID) {
            sp1_2 <- blat.Sp1ToSp2[blat.Sp1ToSp2$percentID >= interID & blat.Sp1ToSp2$percentLength >= interPL,]
            sp2_1 <- blat.Sp2ToSp1[blat.Sp2ToSp1$percentID >= interID & blat.Sp2ToSp1$percentLength >= interPL,]
            dupGenes <- union(sp1_2[duplicated(sp1_2[,1]),1], sp2_1[duplicated(sp2_1[,1]),1])
            sharedGenes <- intersect(sp1_2[,1], sp2_1[,1])
            uniqGenes <- sharedGenes[!(sharedGenes %in% dupGenes)]

            # the blat region is the same as the liftOver region
            inter.blat0 <- cbind(sp2_1[match(uniqGenes, sp2_1[,1]), 1:4], sp1_2[match(uniqGenes, sp1_2[,1]), 2:4])
            inter.orig0 <- cbind(orig.sp1[match(uniqGenes, orig.sp1[,1]),], orig.sp2[match(uniqGenes, orig.sp2[,1]),-1])

            sum(inter.blat0[,3] == inter.orig0[,2] & inter.blat0[,4] == inter.orig0[,3] & inter.blat0[,6] == inter.orig0[,4] & inter.blat0[,7] == inter.orig0[,5])
        })
    })
    
    # Base R plot
    matplot(IDs, exonNums, type = "n", main = "Exon number vs interspecies ID and PL")
    text(IDs, exonNums, rep(IDs, ncol(exonNums)), col = rep(1:length(PLs), each = length(IDs)), cex = .8)
    legend("bottomleft", legend = paste("PL =", PLs), col = 1:length(PLs), lty = 1, cex = .8)

    colnames(exonNums) <- PLs
    rownames(exonNums) <- IDs
    exonNums
}

# Plot the number of exons against the PID and PL used
pdf("/vast/palmer/scratch/noonan/ap2549/XSAnno/BLAT_testing2/blat/blatfilter_by_transcript/exonNum_vs_blat_thresholds_interSp.pdf", 8, 6)
exonNumTable.inter <- chooseThreshold.inter(blat.Sp1ToSp2, blat.Sp2ToSp1, IDs, PLs)
dev.off()

chooseThreshold.intra <- function(blat.Sp1ToSp1, blat.Sp2ToSp2, IDs, PLs) {
    orig.sp1.0 <- do.call(rbind, strsplit(blat.Sp1ToSp1[,1], split = "\\|"))
    orig.sp2.0 <- do.call(rbind, strsplit(blat.Sp2ToSp2[,1], split = "\\|"))

    orig.sp1 <- cbind(apply(orig.sp1.0[,1:5], 1, paste, collapse = "|"), orig.sp1.0[,6:ncol(orig.sp1.0)])
    orig.sp2 <- cbind(apply(orig.sp2.0[,1:5], 1, paste, collapse = "|"), orig.sp2.0[,6:ncol(orig.sp2.0)])

    blat.Sp1ToSp1[,1] <- orig.sp1[,1]
    blat.Sp2ToSp2[,1] <- orig.sp2[,1]
    ########################
    ## filter intraSpecies
    # remove low ID low PL and duplicate regions
    exonNums <- sapply(PLs, function(PL) {
        sapply(IDs, function(ID) {
            sp1_1 <- blat.Sp1ToSp1[blat.Sp1ToSp1$percentID >= ID & blat.Sp1ToSp1$percentLength >= PL,]
            sp2_2 <- blat.Sp2ToSp2[blat.Sp2ToSp2$percentID >= ID & blat.Sp2ToSp2$percentLength >= PL,]
            dupGenes <- union(sp1_1[duplicated(sp1_1[,1]),1], sp2_2[duplicated(sp2_2[,1]),1])
            sharedGenes <- intersect(sp1_1[,1], sp2_2[,1])
            uniqGenes <- sharedGenes[!(sharedGenes %in% dupGenes)]

            # the blat region is the same as the liftOver region
            inter.blat0 <- cbind(sp1_1[match(uniqGenes, sp1_1[,1]), 1:4], sp2_2[match(uniqGenes, sp2_2[,1]), 2:4])
            inter.orig0 <- cbind(orig.sp1[match(uniqGenes, orig.sp1[,1]),], orig.sp2[match(uniqGenes, orig.sp2[,1]),-1])

            sum(inter.blat0[,3] == inter.orig0[,2] & inter.blat0[,4] == inter.orig0[,3] & inter.blat0[,6] == inter.orig0[,4] & inter.blat0[,7] == inter.orig0[,5])
        })
    })
    
    # Base R plot
    matplot(IDs, exonNums, type = "n", main = "Exon number vs intraspecies ID and PL")
    text(IDs, exonNums, rep(IDs, ncol(exonNums)), col = rep(1:length(PLs), each = length(IDs)), cex = .8)
    legend("bottomleft", legend = paste("PL =", PLs), col = 1:length(PLs), lty = 1, cex = .8)

    colnames(exonNums) <- PLs
    rownames(exonNums) <- IDs
    exonNums
}

pdf("/vast/palmer/scratch/noonan/ap2549/XSAnno/BLAT_testing2/blat/blatfilter_by_transcript/exonNum_vs_blat_thresholds_intraSp.pdf", 8, 6)
exonNumTable.intra <- chooseThreshold.intra(blat.Sp1ToSp1, blat.Sp2ToSp2, IDs, PLs)
dev.off()

# Choose interID, interPL, intraID and intraPL, when maximum exon number reached.
interID <- 0.95
interPL <- 0.95
intraID <- 0.97
intraPL <- 0.97

blatFilter <- function(blat.Sp1ToSp1, blat.Sp1ToSp2, blat.Sp2ToSp1, blat.Sp2ToSp2, interID, interPL, intraID, intraPL, sp1Name, sp2Name) {
    # original coordination
    orig.sp1.0 <- do.call(rbind, strsplit(unique(c(blat.Sp1ToSp1[,1], blat.Sp1ToSp2[,1])), split = "\\|"))
    orig.sp2.0 <- do.call(rbind, strsplit(unique(c(blat.Sp2ToSp1[,1], blat.Sp2ToSp2[,1])), split = "\\|"))

    orig.sp1 <- cbind(apply(orig.sp1.0[,1:5], 1, paste, collapse = "|"), orig.sp1.0[,6:ncol(orig.sp1.0)])
    orig.sp2 <- cbind(apply(orig.sp2.0[,1:5], 1, paste, collapse = "|"), orig.sp2.0[,6:ncol(orig.sp2.0)])

    ########################
    ## filter interSpecies
    # remove low ID low PL and duplicate regions
    sp1_2 <- blat.Sp1ToSp2[blat.Sp1ToSp2$percentID >= interID & blat.Sp1ToSp2$percentLength >= interPL,]
    sp1_2[,1] <- apply(do.call(rbind, strsplit(sp1_2[,1], split = "\\|"))[,1:5], 1, paste, collapse = "|")
    sp2_1 <- blat.Sp2ToSp1[blat.Sp2ToSp1$percentID >= interID & blat.Sp2ToSp1$percentLength >= interPL,]
    sp2_1[,1] <- apply(do.call(rbind, strsplit(sp2_1[,1], split = "\\|"))[,1:5], 1, paste, collapse = "|")
    dupGenes <- union(sp1_2[duplicated(sp1_2[,1]),1], sp2_1[duplicated(sp2_1[,1]),1])
    sharedGenes <- intersect(sp1_2[,1], sp2_1[,1])
    uniqGenes <- sharedGenes[!(sharedGenes %in% dupGenes)]

    # the blat region is the same as the liftOver region
    inter.blat0 <- cbind(sp2_1[match(uniqGenes, sp2_1[,1]), 1:4], sp1_2[match(uniqGenes, sp1_2[,1]), 2:4])
    inter.orig0 <- cbind(orig.sp1[match(uniqGenes, orig.sp1[,1]),], orig.sp2[match(uniqGenes, orig.sp2[,1]),-1])

    inter.blat <- inter.orig0[inter.orig0[,4] == inter.blat0[,3] & inter.orig0[,5] == inter.blat0[,4] & inter.orig0[,8] == inter.blat0[,6] & inter.orig0[,9] == inter.blat0[,7],]

    ########################
    ## filter paralogs
    sp1_1 <- blat.Sp1ToSp1[blat.Sp1ToSp1$percentID >= intraID & blat.Sp1ToSp1$percentLength >= intraPL,]
    sp1_1[,1] <- apply(do.call(rbind, strsplit(sp1_1[,1], split = "\\|"))[,1:5], 1, paste, collapse = "|")
    sp2_2 <- blat.Sp2ToSp2[blat.Sp2ToSp2$percentID >= intraID & blat.Sp2ToSp2$percentLength >= intraPL,]
    sp2_2[,1] <- apply(do.call(rbind, strsplit(sp2_2[,1], split = "\\|"))[,1:5], 1, paste, collapse = "|")
    intra.dupGenes <- union(sp1_1[duplicated(sp1_1[,1]),1], sp2_2[duplicated(sp2_2[,1]),1])
    intra.sharedGenes <- intersect(sp1_1[,1], sp2_2[,1])
    intra.uniqGenes <- intra.sharedGenes[!(intra.sharedGenes %in% intra.dupGenes)]

    ############################
    ## output
    out <- inter.blat[inter.blat[,1] %in% intra.uniqGenes,]
    colnames(out) <- c("ID", paste(rep(c(sp1Name, sp2Name), each = 4), rep(c("chr", "strand", "start", "end"), 2), sep = "."))
    return(out)
}

blatFiltered <- blatFilter(blat.Sp1ToSp1, blat.Sp1ToSp2, blat.Sp2ToSp1, blat.Sp2ToSp2, interID, interPL, intraID, intraPL, sp1Name = sp1, sp2Name = sp2)

# Set the working directory and write the output files
setwd("/vast/palmer/scratch/noonan/ap2549/XSAnno/BLAT_testing2/blat/blatfilter_by_transcript")
write.table(blatFiltered, paste("blatFiltered.", sp1, "To", sp2, ".txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(cbind(blatFiltered[, c(2, 4, 5, 1)], rep(1, nrow(blatFiltered)), blatFiltered[, 3]), paste("blatFiltered.", sp1, "To", sp2, ".", sp1, ".bed", sep = ""), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(cbind(blatFiltered[, c(6, 8, 9, 1)], rep(1, nrow(blatFiltered)), blatFiltered[, 7]), paste("blatFiltered.", sp1, "To", sp2, ".", sp2, ".bed", sep = ""), quote = F, sep = "\t", col.names = F, row.names = F)

q(save = "no")
