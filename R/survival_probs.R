#'@title survival_probs
#'@description compute and plot predictions of survival probabilities
#'@param year_ahead numeric, the year in which the projection of the survival probabilities is conducted
#'@param pred_dist  a matrix with samples (columns) from the predictive distribution for the year_ahead of mortality rates at each age (rows).
#'@param ages numeric vector with ages at the year_ahead.
#'@param s numeric, the prediction of survival probability is defined as the probability that a person
#'aged z at the year_ahead will survive up to age z+s.
#'@param plot if TRUE then the boxplots of the predictive distribution of the survival probabilities are plotted
#'
#'@details
#'See also in \insertCite{alexopoulos2019bayesian;textual}{BayesDemog} for details.
#'
#'
#'@return The value returned is a list holding
#'#'\itemize{
#'\item{\code{survival_probs}}  {a matrix with predictive distribution of survival probabilities}
#'\item{\code{survival_plot}}  {boxplots of the predictive distribution of survival probabilities}
#'}
#'
#'@references
#'  \insertAllCited{}
#'
#'
#'
#'@examples
#'
#'\dontrun{
#'ages = c(0,1:89)
#'#we will use the years 1991:2000 to fit the models and
#'#the years 2001:2010 for out-of-sample tests
#'years =  1991:2010
#'country = "GBRTENW"
#' #please provide your HMD user otherwise do not continue
#'HMDuser = "user@..."
#'#please provide the password that corresponds to HMDuser otherwise do not continue
#'HMDpass = ""
#'data_label = "England & Wales Total Population"
#'data_series = "female"
#'
#'mydata =load_demog_data(ages,years,country,HMDuser,HMDpass,data_label,data_series,TRUE,FALSE)
#'
#'its = 5*10^5
#'thin = 200
#'burn = 10^5
#'adapt = burn-10000
#'monitor =1000
#'
#'myD = mydata$deaths
#'myE = mydata$exposures
#'
#'HPmcmc = HP(myD,myE,its,thin,burn,df=2,A=10^10,adapt,monitor,ages,years)
#'HPforecast = HPpred(c(0,1:89),10,10,2001:2010,HPmcmc,c(0.025,0.975))
#'
#'
#'
#'x = survival_probs(2005,exp(HPforecast$pred_dist[5,,]),seq(0,80,by=10),5,TRUE)
#'print(x$survival_plot)
#'}
#'
#'
#'@importFrom utils stack
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 ggtitle
#'@importFrom ggplot2 aes
#'@importFrom ggplot2 geom_boxplot
#'@importFrom ggplot2 coord_cartesian
#'@importFrom ggplot2 theme
#'@importFrom ggplot2 element_text
#'@importFrom ggplot2 labs
#'
#'
#'
#'
#'@export

survival_probs = function(year_ahead,pred_dist,ages,s,plot=T){

  ind = values = NULL
  num_mcmc_samples = dim(pred_dist)[2]
  surv<-matrix(NA,length(ages), num_mcmc_samples)
  for(i in 1:num_mcmc_samples){
    jj<-0
    for(j in (ages+1)){
      jj<-jj+1
      surv[jj,i]<-prod( 1-pred_dist[j:(j+(s-1)),i]      )
    }
  }
  surv<-t(surv)
  colnames(surv)<-as.character(ages)
  if(plot==T){
  dat <- stack(as.data.frame(surv))
  P1<-ggplot(dat) + geom_boxplot(aes(x = ind, y = values))+coord_cartesian(ylim = c(0.5, 1.0))
  P1<-P1+ ggtitle(as.character(year_ahead)) +theme(plot.title = element_text(hjust = 0.5)) +labs(x = "Age",y="Probability")
  }
  list(survival_probs=surv,survival_plot=P1)
}
