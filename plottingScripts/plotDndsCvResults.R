#written by noah friedman
library(ggplot2)
library(ggrepel)

plot_cancer_type_dndscv_results <- function(df, cancerType, title){
  df <- df[df$cancerType == cancerType,]
  df <- df[(df$pancanQVal < .2) | (df$qglobal_cv < .2),]
  
  capThresh <- 1e-10
  minCoord = .9
  maxCoord = 2-log10(capThresh)
  
  plt <- ggplot(df, aes(x=1-log10(pancanQVal), y=1-log10(qglobal_cv)))+
    geom_point()+
    geom_text_repel(aes(label=displayName, colour=KnownCancerGene))+
    
    geom_segment(aes(x=minCoord, xend=maxCoord, y=1- log10(.01), yend= 1- log10(.01)), colour='black', linetype=2)+
    geom_segment(aes(x=1- log10(.01), xend=1- log10(.01), y=minCoord, yend= maxCoord), colour='black', linetype=2)+
    
    geom_segment(aes(x=minCoord, xend=maxCoord, y=1- log10(.1), yend= 1- log10(.1)), colour='black')+
    geom_segment(aes(x=1- log10(.1), xend=1- log10(.1), y=minCoord, yend= maxCoord), colour='black')+
    
    xlab('1 - log10(pan-dog cancer q value)')+
    ylab('1 - log10(cancer type specific q value)')+
    ggtitle(title)
  
  return(plt)
    
}

df <- read.table('/Users/friedman/Desktop/dog/dndsCvPlottingData.tsv', sep = '\t', header=TRUE)
mammary <- plot_cancer_type_dndscv_results(df, 'mammary')
lymphoma <- plot_cancer_type_dndscv_results(df, 'lymphoma')
oralMelanoma <- plot_cancer_type_dndscv_results(df, 'oralMelanoma')
sarcoma <- plot_cancer_type_dndscv_results(df, 'sarcoma')
other <- plot_cancer_type_dndscv_results(df, 'other')

arrangedPlot <- plot_grid(mammary, lymphoma, oralMelanoma, sarcoma, other, nrow=2, ncol=3)
arrangedPlotWithTitle <- plot_grid(ggplot()+ggtitle('DNDS-CV analysis results'),arrangedPlot, nrow=2, rel_heights = c(.2, 1))
ggsave('~/Desktop/plot.pdf', arrangedPlotWithTitle,  width = 18, height = 10, units = c("in"))





#
#####
#########
##########

#TEMPORARY please move e: dndscv results


df <- read.table('/Users/friedman/Desktop/dog/geneMutationsForLolipop.tsv', sep = '\t', header=TRUE)

plot_exon_mutation_counts(df, 'KRAS')


ggsave('~/Desktop/plotKDM6A.pdf', p,  width = 100, height = 5, units = c("in"), limitsize = FALSE)

is.integer(df$orderingVal)

dim(df[(df$gene == 'TP53'),])
dim(df[(df$gene == 'TP53') & (df$exonNumber != 'intron'),])


