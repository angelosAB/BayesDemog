#'@importFrom stats dbinom

HPloglik<-function(theta,y,n,len_ages,len_years,myages){

  ages<-len_ages
  yrs = len_years
  x=myages

  t<-rep(1:yrs,rep(ages,yrs))

  a<-theta[1];b<-theta[2];c<-theta[3];d<-theta[4];e<-theta[5];f<-theta[6];g<-theta[7];h<-theta[8]

  a1<-exp(a)/(1+exp(a));b1<-exp(b)/(1+exp(b));c1<-exp(c)/(1+exp(c));d1<-exp(d)/(1+exp(d));e1<-exp(e);f1<-(40*exp(f)+10)/(1+exp(f));g1<-exp(g)/(1+exp(g));h1<-exp(h)

  hp<-a1^((x+b1)^c1)+d1*exp(-e1*(log(x/f1))^2)+g1*(h1^x)
  p<-hp/(1+hp)
  loglik<-sum(dbinom(y,n,p,log=T))
  return(loglik)
}
