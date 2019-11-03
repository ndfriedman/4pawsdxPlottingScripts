#written by Noah Friedman
#an R script for plotting dog mutational signatures by cancer type

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)


plot_signature_waterfall <- function(df, title,  noLegend = TRUE){
  
  plottingLevels <- c('Signature.1','Signature.2','Signature.3', 'Signature.4', 'Signature.5',
                      'Signature.6','Signature.7','Signature.8', 'Signature.9', 'Signature.10',
                      'Signature.11','Signature.12','Signature.13', 'Signature.14', 'Signature.15',
                      'Signature.16','Signature.17','Signature.18', 'Signature.19', 'Signature.20',             
                      'Signature.21','Signature.22','Signature.23', 'Signature.24', 'Signature.25',
                      'Signature.26','Signature.27','Signature.28', 'Signature.29', 'Signature.30' 
                      )
  barColorPalette <- c('#2A52BE','red','#FF1493', 'orange', 'purple',
  '#008000','#00FA9A','#20B2AA', '#6B8E23', '#7CFC00',
  '#F4A460','#A52A2A','#D2B48C', '#CD853F', '#FFEBCD',
  '#FA8072','#DB7093','#B22222', '#FF6347', '#DC143C',             
  '#FFFACD','#FF7F50','#F0E68C', '#FF6347', '#CCCC00',
  '#F0F8FF','#87CEEB','#483D8B', '#FFD700', 'gray' 
  )
  
  
  #c(
  #  "#00DFFF", "#FF0000", "#FF1493", "#FFA500", 
  #  "#FFB6C1", "#267574", "#FFF600", "#ADFF2F",
  #  "#2A52BE","#551A8B","#D3D3D3"
  #)
  
  plt <- ggplot(df)+
    
    #The first bar of the predominant signature in the positive direction
    geom_bar(aes(x = reorder(sampleID, -orderingVal), y=signatureOfInterestMagnitude, 
                 #fill = signatureOfInterestName), stat="identity")+
                 
                  fill = factor(signatureOfInterestName, levels=plottingLevels)), stat="identity")+
    
    #TODO FIX THE THIRD BAR SO IT 
    #THIRD BAR (plot it first so it is covered by the second bar)
    #geom_bar(aes(x = reorder(sampleID, -orderingVal), y=-secondPredominantSigMagnitude - thirdPredominantSigMagnitude, 
    #            #fill = thirdPredominantSigName), stat = "identity")+
    #             fill = factor(thirdPredominantSigName, levels=plottingLevels)), stat = "identity")+ #color the lower signature column by which signature it is                                                             
    
    #The second bar of the second predominant signature in the negative direction
    geom_bar(aes(x = reorder(sampleID, -orderingVal), y=-secondPredominantSigMagnitude, 
                 #fill = secondPredominantSigName), stat = "identity")+
                 fill = factor(secondPredominantSigName, levels=plottingLevels)), stat = "identity")+ #color the lower signature column by which signature it is                                                                                
    
    scale_fill_manual(values=barColorPalette, drop = FALSE)+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ylab('Signature Fraction')+
    scale_y_continuous(limits = c(-1,1), breaks=c(-1, -.5, 0, .5, 1),
                       labels=c('-1'='1', '-5.'='.5', '0'='0', '.5'='.5', '1'='1'))+
    guides(fill=guide_legend(title="Signature"))+
    ggtitle(title)+
    xlab('Sample ID')
  if(noLegend == TRUE){
    plt <- plt + theme(legend.position = 'none')
  }
  return(plt)
}


df <- read.table('/Users/friedman/Desktop/dog/signaturesPlottingFormat.tsv',sep = '\t', header=TRUE)

lymphomaDf = df[df$cancerType == 'Lymphoma',]
osteosarcomaDf = df[df$cancerType == 'Osteosarcoma',]
oralMelanomaDf = df[df$cancerType == 'Oral Melanoma',]
mammaryCancerDf = df[df$cancerType == 'MammaryCancer',]
hemangiosarcomaDf = df[df$cancerType == 'Hemangiosarcoma',]
otherDf = df[df$cancerType == 'Other',]

totalNCases <- dim(mammaryCancerDf)[1]

#CORRECT PLOTS BY the number of cases 
pltLymphoma <- plot_grid(plot_signature_waterfall(lymphomaDf, title=paste('Lymphoma; n cases: ', dim(lymphomaDf)[1])),
                         ggplot(), rel_widths = c(dim(lymphomaDf)[1], totalNCases -  dim(lymphomaDf)[1]), ncol=2)#set to false to just return the waterfall plot
pltOsteosarcoma <-  plot_grid(plot_signature_waterfall(osteosarcomaDf, title=paste('Osteosarcoma; n cases: ', dim(osteosarcomaDf)[1])),
                         ggplot(), rel_widths = c(dim(osteosarcomaDf)[1], totalNCases -  dim(osteosarcomaDf)[1]), ncol=2)#set to false to just return the waterfall plot
pltOralMelanoma <-  plot_grid(plot_signature_waterfall(oralMelanomaDf, title=paste('Oral Melanoma; n cases: ', dim(oralMelanomaDf)[1])),
                         ggplot(), rel_widths = c(dim(oralMelanomaDf)[1], totalNCases -  dim(oralMelanomaDf)[1]), ncol=2) #set to false to just return the waterfall plot
pltMammaryCancer <-  plot_grid(plot_signature_waterfall(mammaryCancerDf, title=paste('Mammary Cancer; n cases: ', dim(mammaryCancerDf)[1])),
                         ggplot(), rel_widths = c(dim(mammaryCancerDf)[1], totalNCases -  dim(mammaryCancerDf)[1]), ncol=2)#set to false to just return the waterfall plot
pltHemangiosarcoma <-  plot_grid(plot_signature_waterfall(hemangiosarcomaDf, title=paste('Hemangiosarcoma; n cases: ', dim(hemangiosarcomaDf)[1])),
                         ggplot(), rel_widths = c(dim(hemangiosarcomaDf)[1], totalNCases - dim(hemangiosarcomaDf)[1], ncol=2))#set to false to just return the waterfall plot
pltOther <-  plot_grid(plot_signature_waterfall(otherDf, title=paste('Other; n cases: ', dim(otherDf)[1])),
                                 ggplot(), rel_widths = c(dim(otherDf)[1], totalNCases - dim(otherDf)[1], ncol=2))#set to false to just return the waterfall plot


legend <- get_legend(plot_signature_waterfall(hemangiosarcomaDf, title=paste('Hemangiosarcoma; n cases: ', dim(hemangiosarcomaDf)[1]), noLegend = FALSE))

fullPlot <- plot_grid(pltLymphoma, pltOsteosarcoma, pltOralMelanoma, pltMammaryCancer, pltHemangiosarcoma, pltOther,
                      nrow=6)

fullPlotWithLegend <- plot_grid(fullPlot, legend, ncol=2, rel_widths = c(1,.1))

fullPlotWithLegendAndTitle <- plot_grid(ggplot()+ggtitle('Landscape of mutational signatures in canine cancer')+theme(title =element_text(size=40, face='bold')),
      fullPlotWithLegend, nrow=2, rel_heights = c(.1,1))

ggsave('~/Desktop/plot.pdf', fullPlotWithLegendAndTitle,  width = 25, height = 30, units = c("in"))

#TODO add cases with cancer type equals other!


