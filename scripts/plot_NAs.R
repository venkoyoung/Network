#desde aca vengo con el script renames
plotNAs<-function(matrixI)
{
forPlot<-matrixI
forPlot[!is.na(forPlot)] <-1
forPlot[is.na(forPlot)] <-4
forPlotMelted<-data.table::melt(data=forPlot)
print(dim(forPlotMelted))
head(forPlotMelted)
table(forPlotMelted$value)
####################################
#only 2 colors
pdf("withNs_hor_A4.pdf", width = 11.69,  height = 8.27)
p<-ggplot(forPlotMelted, 
         aes(x = Var2, y = Var1, fill = factor(value)))  +
    theme_minimal() +
    geom_tile(position = "identity",  
              color = "grey")     +
    scale_fill_manual(values=c("white","red"))  +
    scale_x_discrete(position = "top") +
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
print(p)
  dev.off()
}

  
  