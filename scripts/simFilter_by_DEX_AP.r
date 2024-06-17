# Load the necessary library
library(DESeq2)

# Read the combined count table
countTable <- read.delim("exon_count_matrix_hs_pt.txt", header=TRUE, stringsAsFactors=FALSE)

# Read the blatFiltered file
blatFiltered <- read.delim("../../BLAT_testing2/blat/blatfilter_by_transcript/blatFiltered.hg38TopanTro6.txt", header=TRUE, stringsAsFactors=FALSE)

# Filter the countTable to only include rows present in blatFiltered
exonCount.blat <- countTable[countTable$ID %in% blatFiltered$ID,]

# Filter the countTable based on counts
countTable <- exonCount.blat[apply(exonCount.blat[, 2:21], 1, max) > 10, 2:21]

# Create the coldata dataframe
coldata <- data.frame(
  name = paste(rep(c("hs", "pt"), each=10), rep(1:10, 2), sep=""),
  species = rep(c("hs", "pt"), each=10)
)

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countTable, colData=coldata, design=~species)
rownames(dds) <- exonCount.blat$ID[apply(exonCount.blat[, 2:21], 1, max) > 10]

# Run DESeq
dds <- DESeq(dds)

# Extract results
res <- results(dds, name="species_pt_vs_hs")

# Write the DESeq2 results to a file
write.table(data.frame(ID=rownames(res), res), "simDESeq_blat.exon.results.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# Keep exons that are not significant (padj >= 0.01)
inExons.blat <- exonCount.blat$ID[!(exonCount.blat$ID %in% rownames(res)[res$padj < 0.01])]
simFilteredExon <- blatFiltered[match(inExons.blat, blatFiltered$ID),]

# Write the filtered exons to a file
write.table(simFilteredExon, "simFiltered.hg38TopanTro6.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# Write the BED files for hg38 and panTro6
write.table(cbind(simFilteredExon[, c(2, 4, 5, 1)], rep(1, nrow(simFilteredExon)), simFilteredExon[, 3]),
            "simFiltered.hg38TopanTro6.hg38.bed", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
write.table(cbind(simFilteredExon[, c(6, 8, 9, 1)], rep(1, nrow(simFilteredExon)), simFilteredExon[, 7]),
            "simFiltered.hg38TopanTro6.panTro6.bed", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
