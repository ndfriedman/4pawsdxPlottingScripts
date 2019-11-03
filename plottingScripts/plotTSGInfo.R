#written by noah friedman

library(ggrepel)

df <- read.table('~/Desktop/dog/truncatingMutationsInTSGs.tsv',  sep = '\t', header=TRUE)

ggplot(df[df$size < 1e5,], aes(x=size, y=nOccurences))+
  geom_text_repel(aes(label=gene))+
  scale_x_log10()

ggplot(df[df$size < 1e5,], aes(x=reorder(gene, nOccurences/size), y=nOccurences/size))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=60, hjust=.5))+
  ylab('N occurences in pancan cohort/gene size')+
  xlab('gene')+
  ggtitle('Truncating Mutations in Tumor Suppressors in Dog Cancer')
  #geom_smooth(method='lm')+
  #scale_x_log10()


plot_recurrence_of_muts <- function(df, title, lim = 2){
  plotDf <- df[df$nmut > lim,]
  plt <- ggplot(plotDf, aes(x=reorder(gene, nmut/geneSize), y=nmut/geneSize))+
    geom_bar(stat = 'identity')+
    theme(axis.text.x = element_text(angle=60, hjust=.5))+
    ggtitle(title)+
    ylab('nmut in cancer type\n normed by gene size')+
    xlab('gene')
  return(plt)
}

df2 <- read.table('~/Desktop/dog/tsgMutCountsByCancerType.tsv',  sep = '\t', header=TRUE)
p1 <- plot_recurrence_of_muts(df2[df2$cancerType == 'Lymphoma',], 'Lymphoma')
p2 <- plot_recurrence_of_muts(df2[df2$cancerType == 'Oral Melanoma',], 'Oral Melanoma')
p3 <- plot_recurrence_of_muts(df2[df2$cancerType == 'MammaryCancer',], 'MammaryCancer')
p4 <- plot_recurrence_of_muts(df2[df2$cancerType == 'Soft Tissue Sarcoma',], 'Soft Tissue Sarcoma', lim=10)
p5 <- plot_recurrence_of_muts(df2[df2$cancerType == 'Osteosarcoma',], 'Osteosarcoma', lim=0)
p6 <- plot_recurrence_of_muts(df2[df2$cancerType == 'Hemangiosarcoma',], 'Hemangiosarcoma', lim=0)

alignedP <- plot_grid(p1, p2, p3, p4, p5, p6, nrow=3, ncol=2)

fP <- plot_grid(ggplot() + ggtitle('Truncating mutations in Tumor Suppressors in Dog Cancer'),
                       alignedP, nrow=2, rel_heights=c(.2,1))

ggsave('~/Desktop/plot.pdf', plot=fP,  width = 10, height = 10, units = c("in"))

#
####
#########
################
#######################
################################
######################
#################
##########
#####
#

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/tsgFracsPOLE_MMR.tsv', sep = '\t', header=TRUE)

p <- ggplot(df, aes(x=class, y=fracClass))+
  geom_bar(stat='identity', aes(fill=geneType))+
  scale_fill_manual(values = c('black', 'gray'))+
  emptyTheme+
  xlab('Mutation Type')+
  theme(axis.text.x = element_text(angle=90))+
  geom_segment(aes(x=0, xend=3, y=.558, yend=.558), colour='white')+
  geom_text(aes(x=1.5, y=.5, label='Fraction of panel\nCDS in TSGs'), size=2)+
  ylab('Frac mutations')+
  ggtitle('All truncating mutations\nin hypermutators')+
  labs(caption='plotTSGInfo.R\nsummarize_types_of_genes_mutated_hypermutators.ipynb')+
  theme(legend.position = 'none')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 3, height = 5, units = c("in"))

#####
df2 <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/tsgFracsRecurrentAllelesPOLE_MMR.tsv', sep = '\t', header=TRUE)

p2 <- ggplot(df2, aes(x=class, y=fracClass))+
  geom_bar(stat='identity', aes(fill=geneType))+
  scale_fill_manual(values = c('black', 'gray'))+
  emptyTheme+
  xlab('Mutation Type')+
  theme(axis.text.x = element_text(angle=90))+
  geom_segment(aes(x=0, xend=3, y=.558, yend=.558), colour='white')+
  geom_text(aes(x=1.5, y=.5, label='Fraction of panel\nCDS in TSGs'), size=2)+
  ylab('Frac mutations')+
  ggtitle('Recurrent truncating alleles')+
  labs(caption='plotTSGInfo.R\nsummarize_mutations_seen_and_never_seen.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p2,  width = 3, height = 5, units = c("in"))



#
#####
########
#############
###################
########################
##################
############
#########
#####
#

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/alleleRecurrences.tsv', sep = '\t', header=TRUE)

p <- ggplot(df, aes(x=reorder(alleleClass, orderingVal), y=fracOccurencesMultiplet))+
  stat_summary()+
  xlab('Type of allele')+
  ylab('Fraction of times allele is observed\nin conjunction with another oncogenic mutation')+
  ggtitle('Double oncogenic mutation\n by allele type')+
  theme(axis.text.x = element_text(angle=90))+
  emptyTheme+
  labs(caption = 'plotTSGInfo.R\nmultiple_mutation_analyses.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 2.5, height = 5, units = c("in"))



