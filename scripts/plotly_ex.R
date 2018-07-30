p<-
  plot_ly(dataA, 
          x = ~x, 
          y = ~y, 
          z = ~cut, 
          color = ~cut, 
          colors = heat.colors(5),
          marker = list(symbol = 'dot', size = 0.5)) 

toy<- matrix(rep(c(0,1),each=5), ncol=5, byrow = T)
toy<- matrix(0, ncol=2,nrow=2)
toy2<-toy+1
plot_ly(showscale = FALSE) %>%
  add_surface(z = toy)%>%
  add_surface(z = toy2)
toy<-matrix(rep(1:2,2), nrow=2, ncol=2, byrow = T)
plot_ly(showscale = FALSE) %>%
  add_surface(z = toy)
