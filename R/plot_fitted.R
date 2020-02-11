#'@title plot_fitted
#'@description plots fitted mortality rates by using either the Heligman-Pollard or the Gaussian Markov random field model.
#'@param output a DemogMCMC object
#'@param years_fitted years used to fit the model
#'@param years_plot years to plot
#'@param quant numeric vector with quantiles that form the credible interval of the fitted mortality rates
#'@param ages numeric vector with the ages
#'@param deaths  (number of ages) times (number of years) matrix with deaths
#'@param exposures  (number of ages) times (number of years) matrix with exposures to risk
#'@param log  indicates if the rates should be plotted in the log-scale, default is TRUE
#
#'
#'@examples
#'
#'
#'\dontrun{
#'ages = c(0,1:89)
#'years =  1991:2000
#'country = "GBRTENW"
#'#please provide your HMD user otherwise do not continue
#'HMDuser = "user@..."
#'#please provide the password that corresponds to HMDuser otherwise do not continue
#'HMDpass = " "
#'data_label = "England & Wales Total Population"
#'data_series = "female"
#'
#'mydata =load_demog_data(ages,years,country,HMDuser, HMDpass,data_label,data_series,TRUE,FALSE)
#'
#'
#'myD = mydata$deaths
#'myE = mydata$exposures
#'
#'#Example 1: Fitted mortality rates by using the Heligman-Pollard model
#'
#'its = 5*10^5
#'thin = 200
#'burn = 10^5
#'adapt = burn-10000
#'monitor =1000
#'
#'HPmcmc = HP(myD,myE,its,thin,burn,df=2,A=10^10,adapt,monitor,ages,years)
#'#plots for all the years used to fit the model
#'plot_fitted(HPmcmc,1991:2000,1991:2000,c(0.025,0.975),c(0,1:89),myD,myE,TRUE)
#'#plots for 1991 and 1995
#'plot_fitted(HPmcmc,1991:2000,c(1991,1995),c(0.025,0.975),c(0,1:89),myD,myE,TRUE)
#'
#'
#'#Example 2: Fitted mortality rates by using the Gaussian Markov random field model
#'its=15000
#'thin=10
#'burn=5000
#'adapt=3000
#'monitor=1000
#'
#'GMRFmcmc = GMRF(myD,myE,its,thin,burn,adapt,monitor,ages,years)
#'
#'
#'#plots for all the years used to fit the model
#'plot_fitted(GMRFmcmc,1991:2000,1991:2000,c(0.025,0.975),c(0,1:89),myD,myE,TRUE)
#'#plots for 1991 and 1995
#'plot_fitted(GMRFmcmc,1991:2000,c(1991,1995),c(0.025,0.975),c(0,1:89),myD,myE,TRUE)
#'}
#'
#'
#'@importFrom methods is
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 ggtitle
#'@importFrom ggplot2 theme
#'@importFrom ggplot2 element_text
#'@importFrom ggplot2 aes
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 labs
#'
#'
#'
#'
#'
#'@export
#'
#'
#'
#'
plot_fitted<-function(output,years_fitted,years_plot,quant,ages,deaths,exposures,log=T){

  age = p.mean=p.quant.up=p.quant.low=trues = quantile = NULL
  if (!is(output, "DemogMCMC")) stop("This function expects a 'DemogMCMC' object.")

  colnames(output$mortality_rates) = as.character(years_fitted)
  for(j in years_plot){
    i = as.character(j)
    if(log==T){
      true_rate = log(deaths[,i]/exposures[,i])
      rate_dist = log(output$mortality_rates[,i,])
      ytitle = "Probability of death (log-scale)"
    }else{
      true_rate = deaths[,i]/exposures[,i]
      rate_dist = output$mortality_rates[,i,]
      ytitle = "Probability of death"
    }

      df<-data.frame(age=rep(ages,1),p.mean=apply(rate_dist,1,mean),p.quant.up=apply(rate_dist,1,quantile,quant[2]),p.quant.low=apply(rate_dist,1,quantile,quant[1]),trues = true_rate)
      P2<-ggplot()+ ggtitle(i) +
        theme(plot.title = element_text(hjust = 0.5)) +labs(x = "Age",y=ytitle) +geom_line(data = df, aes(x = age, y = p.quant.up), color = "#F8766D",linetype=22)+geom_line(data = df, aes(x = age, y = p.quant.low), color = "#F8766D",linetype=22)+
        geom_line(data = df, aes(x = age, y = p.mean), color = "#F8766D")
      P2 = P2 + geom_point(data=df, aes(x = age, y = trues))
    print(P2)
  }
}
