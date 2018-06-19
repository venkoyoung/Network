 #finalMatrixCOmparison
library(data.table)
library(stringr)
forPlot<-dfOrigEvents
#row.names(forPlot)<-rownames(dfOrigEvents)
row.names(forPlot)<-1:31
colnames(forPlot)<-colnames(dfOrigEvents)
length(which(is.na(forPlot)))#681
forPlot[is.na(forPlot)] <-4
forPlot[1:5,1:9]
#SD controles
forPlotMelted<-melt(forPlot)
head(forPlotMelted)
table(forPlotMelted$value)
#xlabels
eventName<-row.names(dfOrigEvents)
head(dfOrigEvents)
head(forPlotMelted)
table(forPlotMelted$value)
####################################
#pdf('withNs_hor_large.pdf', width = 80, height = 10)
pdf("withNs_hor_A4.pdf", width = 11.69,  height = 8.27)
  ggplot(forPlotMelted, 
         aes(x = Var2, y = Var1, fill = factor(value)))  +
    theme_minimal() +
    geom_tile(position = "identity",  
              color = "grey") +
    scale_y_continuous(breaks=c(1:31), labels=eventName,
                       limits=c(0,32) , expand = c(0,-0.5)) +
    scale_x_discrete(position = "top") +
    scale_fill_manual(values=c("white","blue","green","pink","red"))+
    theme(
      axis.text.x = element_text(angle = 90, hjust=0, vjust=-5,size=3),
      axis.ticks.length=unit(0,"cm"),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border=element_blank(),
      legend.position="none",
      axis.line=element_blank()
    )
  dev.off()

  
  