# author: Colin
#
# function to make a volcano plot for MAGeCK gene summary output file. x and y lim percentile parameters set the percentage of data points for each axis to plot either as
# a circle or triangle and placed at the graph limits.
#

library("ggbiplot")
library("dplyr")
library("ggplot2")
library("ggrepel")


mageck.volcano.plot<-function(mageckGeneSum, alphavalue = 0.05, lfcLim = log2(1.5), sigLabel = TRUE, genesetFile = NULL,
                                   xlimsPercentile = 0.00, ylimsPercentile = 0.00, width = 6, height = 6) {
  datadf <- read.table(mageckGeneSum, sep = "\t", header = T, stringsAsFactors = F)
  datadf <- na.omit(datadf)
  datadf$log2FoldChange <- datadf[,14]
  datadf$padj <- apply(subset(datadf, select=c(5,11)), 1, min)
  datadf$GeneSymbol <- datadf[,1]
  graphname <- unlist(strsplit(mageckGeneSum, split = "\\."))[1]
  outputName <- paste(graphname, "_volcanoPlot.pdf", sep = "")
  
  # manages gene labels
  input <- mutate(datadf, sig = ifelse(datadf$padj < alphavalue & abs(datadf$log2FoldChange) >= lfcLim, "Sig", "Not_Sig")) #Will have different colors depending on significance
  input$label <- ""
  if (is.null(genesetFile) == F) {
    geneset <- read.table(genesetFile, stringsAsFactors = F)
    geneset$rownum <- match(geneset[,1], datadf$GeneSymbol)
    geneset <- na.omit(geneset)
    input[geneset[,2],length(input)] <- as.character(geneset[,1])
  } else if (sigLabel == TRUE && is.null(genesetFile) == T) {
    input$label[which(input$sig == "Sig")] <- input$GeneSymbol[which(input$sig == "Sig")]
  }
  # determines x and y limits
  upperFDR <- -log10(quantile(input$padj, ylimsPercentile))
  # for equal x axis both sides
  #upperlfc <- max(abs(quantile(input$log2FoldChange, c(xlimsPercentile, 1-xlimsPercentile))))
  #xlims <- c(-upperlfc, upperlfc)
  lfclims <- quantile(input$log2FoldChange, c(xlimsPercentile, 1-xlimsPercentile))
  xlims <- unname(lfclims)
  ylims <- c(0, max(ceiling(upperFDR), 2))
  
  # changes outliers to triangles
  input$shape <- ifelse(-log10(input$padj) > ylims[2] | input$log2FoldChange > xlims[2] | input$log2FoldChange < xlims[1], "triangle", "circle")
  input$padj[-log10(input$padj) > ylims[2]] <- 10^-(ylims[2])
  input$log2FoldChange[input$log2FoldChange > xlims[2]] <- xlims[2]
  input$log2FoldChange[input$log2FoldChange < xlims[1]] <- xlims[1]
  
  # begins plotting
  v <- ggplot(input, aes(log2FoldChange, -log10(padj))) #volcanoplot with log2Foldchange versus pvalue
  v + geom_point(aes(col = sig, shape=shape), alpha = 0.8) + #add points colored by significance
    scale_color_manual(values = c("grey72","red")) + xlim(xlims) + ylim(ylims) +
    xlab("log2(fold change)") + ylab("-log10(FDR)") + ggtitle(graphname) + theme_bw() + 
    geom_hline(yintercept = -log10(alphavalue), color="black", linetype="dotted", size=0.5) + geom_vline(xintercept = c(-lfcLim, lfcLim), color="black", linetype="dotted", size=0.5) +
    theme(legend.position = "none", axis.text = element_text(size=14), axis.title = element_text(size=16,face="bold")) + ggtitle(graphname) +
    geom_text_repel(label = input$label, box.padding = 0.5, size = 5, force = 10, segment.size = 0.3, min.segment.length = 0.01, segment.alpha = 0.5, max.overlaps = 100000) #adding text for the top 20 genes
  ggsave(outputName, width = width, height = height)
}
