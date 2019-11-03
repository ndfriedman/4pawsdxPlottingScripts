#written by noah friedman

#written by noah friedman
library(ggplot2)
library(ggrepel)


plot_exon_mutation_counts <- function(df, gene, hideLegend=TRUE){
  
  pal <- rep(c('black', 'gray'), 200) #alternating black and gray color
  
  dfRed <- df[(df$gene == gene) &  (df$exonNumber != 'intron'),]
  
   p <- ggplot(dfRed, aes(x=exon))+
    
    geom_segment(aes(x=exonCumSum, xend=exonCumSumPrev, y=0, yend=0, colour=exon))+
    geom_bar(aes(x=(exonCumSum + exonCumSumPrev)/2,
                 fill=consequenceSimplified), width =10)+ 
    
    geom_text(aes(x=(exonCumSum + exonCumSumPrev)/2, y=-2, label=exonNumber), angle = 90)+
    
    scale_fill_brewer(palette = 'Set2', drop=FALSE)+
    scale_color_manual(values = pal)+
    theme(axis.text.x = element_text(angle=60))+
    ggtitle(gene)+
    ylab('N mutations at exon')+
    expand_limits(y=-5)+ #MAKE THE CHART GO NEGATIVE
    
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    xlab('cds number')
  
  #LEGEND info
  if(hideLegend){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}


#version corrected for ilya with mutation counts corrected by exon size
plot_exon_mutation_counts_v2 <- function(df, gene, hideLegend=TRUE){
  
  pal <- rep(c('black', 'gray'), 200) #alternating black and gray color
  
  dfRed <- df[(df$gene == gene) &  (df$exonNumber != 'intron'),]
  
  #ggplot(dfRed, aes(x=reorder(exon, orderingVal)))+
  p <- ggplot(dfRed, aes(x=exon))+
    
    geom_segment(aes(x=exonCumSum, xend=exonCumSumPrev, y=0, yend=0, colour=exon))+
    geom_bar(aes(x=(exonCumSum + exonCumSumPrev)/2, y=1.0/exonSize,
                 fill=consequenceSimplified), width =10, stat='identity')+ 
    
    ylim(-.1,1)+
    #geom_text(aes(x=(exonCumSum + exonCumSumPrev)/2, y=-2, label=exonNumber), angle = 90)+
    geom_text(aes(x=(exonCumSum + exonCumSumPrev)/2, y=-.05, label=exonNumber), angle = 90)+
    
    
    scale_fill_brewer(palette = 'Set2', drop=FALSE)+
    scale_color_manual(values = pal)+
    theme(axis.text.x = element_text(angle=60))+
    ggtitle(gene)+
    ylab('N mutations at exon')+
    #expand_limits(y=-5)+ #MAKE THE CHART GO NEGATIVE
    
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank())+
    xlab('cds number')
    
    #LEGEND info
  if(hideLegend){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}

make_big_plot <- function(df){
  
  maxWidth <- my.max(df$exonCumSum)
  spicyP <- plot_grid(
    plot_grid(plot_exon_mutation_counts(df, 'TP53'),
              rel_widths = c(my.max(df[df['gene'] == 'TP53',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'TP53',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'KRAS'),
              rel_widths = c(my.max(df[df['gene'] == 'KRAS',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'KRAS',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'PIK3CA'),
              rel_widths = c(my.max(df[df['gene'] == 'PIK3CA',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'PIK3CA',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'ARID1A'),
              rel_widths = c(my.max(df[df['gene'] == 'ARID1A',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'ARID1A',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'PTEN'),
              rel_widths = c(my.max(df[df['gene'] == 'PTEN',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'PTEN',]$exonCumSum)), ncol=2),
    
    plot_grid(plot_exon_mutation_counts(df, 'EGFR'),
              rel_widths = c(my.max(df[df['gene'] == 'EGFR',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'EGFR',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'NF1'),
              rel_widths = c(my.max(df[df['gene'] == 'NF1',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'NF1',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'ATM'),
              rel_widths = c(my.max(df[df['gene'] == 'ATM',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'ATM',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'BRAF'),
              rel_widths = c(my.max(df[df['gene'] == 'BRAF',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'BRAF',]$exonCumSum)), ncol=2),
    
    plot_grid(plot_exon_mutation_counts(df, 'IDH1'),
              rel_widths = c(my.max(df[df['gene'] == 'IDH1',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'IDH1',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'NRAS'),
              rel_widths = c(my.max(df[df['gene'] == 'NRAS',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'NRAS',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'AKT1'),
              rel_widths = c(my.max(df[df['gene'] == 'AKT1',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'AKT1',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'CTNNB1'),
              rel_widths = c(my.max(df[df['gene'] == 'CTNNB1',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'CTNNB1',]$exonCumSum)), ncol=2),
    #plot_grid(plot_exon_mutation_counts(df, 'APC'),
    #         rel_widths = c(my.max(df[df['gene'] == 'APC',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'APC',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'ERBB2'),
              rel_widths = c(my.max(df[df['gene'] == 'ERBB2',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'ERBB2',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'ERBB3'),
              rel_widths = c(my.max(df[df['gene'] == 'ERBB3',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'ERBB3',]$exonCumSum)), ncol=2),
    
    #lymphoma/leukemia/other
    plot_grid(plot_exon_mutation_counts(df, 'KMT2D'),
              rel_widths = c(my.max(df[df['gene'] == 'KMT2D',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'KMT2D',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'CREBBP'),
              rel_widths = c(my.max(df[df['gene'] == 'CREBBP',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'CREBBP',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'TET2'),
              rel_widths = c(my.max(df[df['gene'] == 'TET2',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'TET2',]$exonCumSum)), ncol=2),
    #plot_grid(plot_exon_mutation_counts(df, 'DNMT3A'),
    #          rel_widths = c(my.max(df[df['gene'] == 'DNMT3A',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'DNMT3A',]$exonCumSum)), ncol=2),
    #plot_grid(plot_exon_mutation_counts(df, 'SETD2'),
    #          rel_widths = c(my.max(df[df['gene'] == 'SETD2',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'SETD2',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'KDM5C'),
              rel_widths = c(my.max(df[df['gene'] == 'KDM5C',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'KDM5C',]$exonCumSum)), ncol=2),
    plot_grid(plot_exon_mutation_counts(df, 'KDM6A'),
              rel_widths = c(my.max(df[df['gene'] == 'KDM6A',]$exonCumSum), maxWidth - my.max(df[df['gene'] == 'KDM6A',]$exonCumSum)), ncol=2),
    
    nrow=22
  )
  return(spicyP)
}

plot_exon_specific_mut_sprectrum <- function(df, gene, exonNum){
  dfRed <- df[(df$gene == gene) & (df$exonNumber == exonNum),]
  ggplot(dfRed, aes(x=POS))+
    geom_bar(aes(fill=consequenceSimplified))+
    ggtitle(gene)
    #geom_bar(aes(fill=ALT))
}

my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

#MAKE THE BIG PLOT WHICH IS ABOUT MUTATION RECURRENCE
df <- read.table('/Users/friedman/Desktop/dog/geneMutationsForLolipop.tsv', sep = '\t', header=TRUE)
spicyP <- make_big_plot(df)
legend <- get_legend(plot_exon_mutation_counts(df, 'KRAS', hideLegend=FALSE)
                              )
finalPlot <- plot_grid(spicyP, legend, ncol = 2, rel_widths = c(1,.2))
ggsave('~/Desktop/plot.pdf', plot=finalPlot,  width = 40, height = 80, units = c("in"), limitsize = FALSE)


#
###
########
#############
##################
##############
########
###
#
#
#



#Plots of exons of interest:
p <- plot_grid(
plot_exon_specific_mut_sprectrum(df, 'TP53', 3), plot_exon_specific_mut_sprectrum(df, 'TP53', 4), plot_exon_specific_mut_sprectrum(df, 'TP53', 5), plot_exon_specific_mut_sprectrum(df, 'TP53', 8),
plot_exon_specific_mut_sprectrum(df, 'NRAS', 1), plot_exon_specific_mut_sprectrum(df, 'NRAS', 2), plot_exon_specific_mut_sprectrum(df, 'NRAS', 3),plot_exon_specific_mut_sprectrum(df, 'NRAS', 4),
plot_exon_specific_mut_sprectrum(df, 'KRAS', 1), plot_exon_specific_mut_sprectrum(df, 'KRAS', 3), plot_exon_specific_mut_sprectrum(df, 'KRAS', 4), ggplot(),
#plot_exon_specific_mut_sprectrum(df, 'PIK3CA', 4), plot_exon_specific_mut_sprectrum(df, 'PIK3CA', 5), plot_exon_specific_mut_sprectrum(df, 'PIK3CA', 6), ggplot(df, 'PIK3CA', 12),
plot_exon_specific_mut_sprectrum(df, 'ARID1A', 5), plot_exon_specific_mut_sprectrum(df, 'ARID1A', 19), ggplot(), ggplot(),
plot_exon_specific_mut_sprectrum(df, 'EGFR', 12), ggplot(), ggplot(), ggplot(),
plot_exon_specific_mut_sprectrum(df, 'NF1', 17), ggplot(), ggplot(), ggplot(),
plot_exon_specific_mut_sprectrum(df, 'BRAF', 6), ggplot(), ggplot(), ggplot(),
plot_exon_specific_mut_sprectrum(df, 'IDH1', 7), ggplot(), ggplot(), ggplot(),
plot_exon_specific_mut_sprectrum(df, 'CTNNB1', 2), plot_exon_specific_mut_sprectrum(df, 'CTNNB1', 12), ggplot(), ggplot(),
plot_exon_specific_mut_sprectrum(df, 'KMT2D', 26), plot_exon_specific_mut_sprectrum(df, 'KMT2D', 28), plot_exon_specific_mut_sprectrum(df, 'KMT2D', 32), plot_exon_specific_mut_sprectrum(df, 'KMT2D', 48),
plot_exon_specific_mut_sprectrum(df, 'CREBBP', 31), ggplot(), ggplot(), ggplot(),
plot_exon_specific_mut_sprectrum(df, 'TET2', 2), plot_exon_specific_mut_sprectrum(df, 'TET2', 9), ggplot(), ggplot(),
plot_exon_specific_mut_sprectrum(df, 'KDM5C', 3), plot_exon_specific_mut_sprectrum(df, 'KDM5C', 18), plot_exon_specific_mut_sprectrum(df, 'KDM5C', 29), plot_exon_specific_mut_sprectrum(df, 'KDM5C', 30),
plot_exon_specific_mut_sprectrum(df, 'KDM6A', 1), plot_exon_specific_mut_sprectrum(df, 'KDM6A', 3), plot_exon_specific_mut_sprectrum(df, 'KDM6A', 6), plot_exon_specific_mut_sprectrum(df, 'KDM6A', 20),
ncol=4, nrow=14)

#pT <- plot_grid(ggplot()+ggtitle('NRAS mutations'), p, rel_heights = c(.5,1), nrow=2)

ggsave('~/Desktop/plot.pdf', plot=p,  width = 20, height = 50, units = c("in"), limitsize = FALSE)

