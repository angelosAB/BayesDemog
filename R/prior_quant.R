prior_quant<-function(low,up){
  s<- (up-low)/4.66
  m<- low+ 2.33*s
  return(c(m,s))
}
