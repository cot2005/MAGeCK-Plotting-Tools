# author: Colin
# function to make a generic rank plot from a mageck gene summary output file.
#

library("ggbiplot")
library("dplyr")
library("ggplot2")
library("ggrepel")


mageck.rankplot<-function(mageckGeneSum, hits = 3, width = 6, height = 6) {
  graphname <- unlist(strsplit(mageckGeneSum, split = "\\."))[1]
  outputName <- paste(graphname, "_rankplot.pdf", sep = "")
  datatable <- read.table(mageckGeneSum, sep = "\t", header = T, stringsAsFactors = F)
  datatable <- subset(datatable, select = c(1,14))
  colnames(datatable) <- c("gene", "lfc")
  datatable <- datatable[order(datatable$lfc),]
  datatable <- na.omit(datatable)
  genenum <- nrow(datatable)
  datatable$rank <- 1:genenum
  labeltext <- rep("", genenum)
  if (hits > 0) {
    rowIndex <- c(seq(1,hits), seq(genenum - hits + 1, genenum))
    generows <- data.frame(gene = datatable$gene[rowIndex], index = rowIndex)
    labeltext[generows[,2]] <- as.character(generows[,1])
  }
  ylimit <- max(abs(min(datatable$lfc)), abs(max(datatable$lfc)))
  ggdata <- ggplot(datatable, aes(x = rank, y = lfc, color=lfc))
  graph <- ggdata + geom_hline(yintercept = 0, color="black", linetype="dashed", size=0.5) +
    geom_point(size = 2) + 
    scale_color_gradientn(colors = c("steelblue","grey","red"), 
                          values = scales::rescale(c(min(datatable$lfc), -0.5, 0, 0.5, max(datatable$lfc)))) + 
    ylab("Log2FC") + xlab("sgRNA Rank") + theme_bw() + 
    theme(axis.text = element_text(size=18), axis.title = element_text(size=18, face="bold"), legend.position = "none") +
    ggtitle(graphname) + 
    ylim(c(min(datatable$lfc), max(datatable$lfc))) + #for automatic lims
    geom_text_repel(label = labeltext, box.padding = 0.4, min.segment.length = 0.25, size = 5.5, max.overlaps = 20, segment.size = 0.3)
  ggsave(outputName, width = width, height = height)
}
