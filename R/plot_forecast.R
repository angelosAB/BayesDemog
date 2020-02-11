#'@title plot_forecast
#'@description plots of future mortality rates predicted by using either the Heligman-Pollard or the Gaussian Markov random field model.
#'@param predictive a \code{Pred} object
#'@param quant numeric vector with quantiles that form the credible interval of the fitted mortality rates
#'@param years_ahead numeric vector with years in which to plot the predicted mortality rates
#'@param ages numeric vector with the ages
#'@param deaths  (number of ages) times (number of years) matrix with deaths
#'@param exposures  (number of ages) times (number of years) matrix with exposures to risk
#'@param true_rates logical, if TRUE points with true mortality rates are added in the plot
#'@param log  indicates if the rates should be plotted in the log-scale, default is TRUE
#
#'
#'@examples
#'#Example 1:
#'
#'data(UKandWales)
#'
#'years_fit = 1991:2000
#'myD = mydata$deaths[,as.character(years_fit)]
#'myE = mydata$exposures[,as.character(years_fit)]
#'
#'years_fit =1991:2000
#'myD = mydata$deaths[,as.character(years_fit)]
#'myE = mydata$exposures[,as.character(years_fit)]
#'
#'#Predictions by using the Heligman-Pollard model
#'
#'its = 5*10^5
#'thin = 200
#'burn = 10^5
#'adapt = burn-10000
#'monitor =1000
#'HPmcmc = HP(myD,myE,its,thin,burn,df=2,A=10^10,adapt,monitor,ages,years_fit)
#'HPforecast = HPpred(c(0,1:89),10,10,2001:2010,HPmcmc,c(0.025,0.975))
#'#10 plots should be produced
#'plot_forecast(HPforecast,c(0.025,0.975),2001:2010,c(0,1:89),myD,myE,TRUE,TRUE)
#'#2 plots should be produced
#'plot_forecast(HPforecast,c(0.025,0.975),2001:2002,c(0,1:89),myD,myE,TRUE,TRUE)
#'
#'#Predictions by using the Gaussian Markov random field model
#'its=15000
#'thin=10
#'burn=5000
#'adapt=3000
#'monitor=1000
#'
#'GMRFmcmc = GMRF(myD,myE,its,thin,burn,adapt,monitor,ages,years_fit)
#'
#'GMRFforecast = GMRFpred(c(0,1:89),10,10,2001:2010,GMRFmcmc,c(0.025,0.975))
#'#10 plots should be produced
#'plot_forecast(GMRFforecast,c(0.025,0.975),2001:2010,c(0,1:89),myD,myE,TRUE,TRUE)
#'#2 plots should be produced
#'plot_forecast(GMRFforecast,c(0.025,0.975),2001:2002,c(0,1:89),myD,myE,TRUE,TRUE)
#'
#'#Example 2:
#'
#'\dontrun{
#'
#'ages = c(0,1:89)
#'#we will use the years 1991:2000 to fit the models and
#'#the years 2001:2010 for out-of-sample tests
#'years =  1991:2010
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
#'years_fit =1991:2000
#'myD = mydata$deaths[,as.character(years_fit)]
#'myE = mydata$exposures[,as.character(years_fit)]
#'
#'#Predictions by using the Heligman-Pollard model
#'
#'its = 5*10^5
#'thin = 200
#'burn = 10^5
#'adapt = burn-10000
#'monitor =1000
#'HPmcmc = HP(myD,myE,its,thin,burn,df=2,A=10^10,adapt,monitor,ages,years_fit)
#'HPforecast = HPpred(c(0,1:89),10,10,2001:2010,HPmcmc,c(0.025,0.975))
#'#10 plots should be produced
#'plot_forecast(HPforecast,c(0.025,0.975),2001:2010,c(0,1:89),myD,myE,TRUE,TRUE)
#'#2 plots should be produced
#'plot_forecast(HPforecast,c(0.025,0.975),2001:2002,c(0,1:89),myD,myE,TRUE,TRUE)
#'
#'#Predictions by using the Gaussian Markov random field model
#'its=15000
#'thin=10
#'burn=5000
#'adapt=3000
#'monitor=1000
#'
#'GMRFmcmc = GMRF(myD,myE,its,thin,burn,adapt,monitor,ages,years_fit)
#'
#'GMRFforecast = GMRFpred(c(0,1:89),10,10,2001:2010,GMRFmcmc,c(0.025,0.975))
#'#10 plots should be produced
#'plot_forecast(GMRFforecast,c(0.025,0.975),2001:2010,c(0,1:89),myD,myE,TRUE,TRUE)
#'#2 plots should be produced
#'plot_forecast(GMRFforecast,c(0.025,0.975),2001:2002,c(0,1:89),myD,myE,TRUE,TRUE)
#'}
#'
#'
#'
#'@importFrom methods is
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 ggtitle
#'@importFrom ggplot2 aes
#'@importFrom ggplot2 theme
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 element_text
#'
#'
#'@export
#'
#'
#'
#'
plot_forecast<-function(predictive,quant,years_ahead,ages,deaths,exposures,true_rates=T,log=T){

  age = p.mean=p.quant.up=p.quant.low=true = quantile = NULL

  if (!is(predictive, "Pred")) stop("This function expects a 'Pred' object.")


  for(j in years_ahead){
    i = as.character(j)
  if(log==T){
    true_rate = log(deaths[,i]/exposures[,i])
    pred_dist = predictive$pred_dist
    ytitle = "Probability of death (log-scale)"
  }else{
    true_rate = deaths[,i]/exposures[,i]
    pred_dist = exp(predictive$pred_dist)
    ytitle = "Probability of death"
  }


    df<-data.frame(age=rep(ages,1),true=true_rate,p.mean=apply(pred_dist,c(1,2),mean)[i,],p.quant.up=apply(pred_dist,c(1,2),quantile,quant[2])[i,],p.quant.low=apply(pred_dist,c(1,2),quantile,quant[1])[i,])
    P2<-ggplot()+ ggtitle(i) +
      theme(plot.title = element_text(hjust = 0.5)) +labs(x = "Age",y=ytitle) +geom_line(data = df, aes(x = age, y = p.quant.up), color = "#F8766D",linetype=22)+geom_line(data = df, aes(x = age, y = p.quant.low), color = "#F8766D",linetype=22)+
      geom_line(data = df, aes(x = age, y = p.mean), color = "#F8766D")
     if(true_rates){
       P2 = P2 +geom_point(data=df, aes(x = age, y = true))
     }


    print(P2)
}
}
