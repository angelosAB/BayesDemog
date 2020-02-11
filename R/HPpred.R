#'@title HPpred
#'@description returns the predictive distribution of mortality rates obtained by fitting the dynamic Heligman-Pollard model
#'@param ages numeric vector with the ages used to fit the model, default to the vector c(0,1:89)
#'@param num_years_fitted number of years used to fit the model, default is 10
#'@param num_years_ahead  number of years to forecast the mortality rates for each age in ages
#'@param years_ahead years in which predictions are required
#'@param HP_output  output of the MCMC sampler for the dynamic Heligman-Pollard model
#'@param pred_quants a numeric vector with two quantiles that will form the credible interval of the predictive distributions
#'
#'@details
#'For details concerning the method used to draw samples from the predictive distribution of mortality rates please see the paper by \insertCite{alexopoulos2019bayesian;textual}{BayesDemog}
#'
#'
#'
#'@return The value returned is a list object of class \code{Pred} holding
#'\itemize{
#'\item{\code{pred_quant_up}}  {a matrix consisted of the upper quantile of the prediction credible interval for the log-mortality rate of each age in the input vector \code{ages}}
#'\item{\code{pred_quant_low}}  {a matrix consisted of the lower quantile of the prediction credible interval for the log-mortality rate of each age in the input vector \code{ages}}
#'\item{\code{pred_mean}}  {a matrix consisted of the mean of the posterior predictive distribution of the log-mortality rate for each age in the input vector \code{ages}}
#'\item{\code{pred_dist}}  {a matrix with samples drawn from the predictive distribution of the log-mortality rate for each age in the input vector \code{ages}}
#'}
#'
#'
#'@references
#'  \insertAllCited{}
#'
#'@examples
#'#Example 1:
#'data(UKandWales)
#'
#'its = 5*10^5
#'thin = 200
#'burn = 10^5
#'adapt = burn-10000
#'monitor =1000
#'
#'years_fit = 1991:2000
#'myD = mydata$deaths[,as.character(years_fit)]
#'myE = mydata$exposures[,as.character(years_fit)]
#'
#'mcmc = HP(myD,myE,its,thin,burn,adapt,monitor,df=2,A=10^10,ages,years_fit)
#'HPforecast = HPpred(c(0,1:89),10,10,2001:2010,mcmc,c(0.025,0.975))
#'
#'#Example 2:
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
#'years_fit = 1991:2000
#'myD = mydata$deaths[,as.character(years_fit)]
#'myE = mydata$exposures[,as.character(years_fit)]
#'
#'HPmcmc = HP(myD,myE,its,thin,burn,df=2,A=10^10,adapt,monitor,ages,years_fit)
#'HPforecast = HPpred(c(0,1:89),10,10,2001:2010,HPmcmc,c(0.025,0.975))
#'}
#'
#'
#'@importFrom methods is
#'@importFrom TruncatedNormal rtmvnorm
#'
#'
#'
#'@export
#'
#'
#'
#'
HPpred<-function(ages,num_years_fitted,num_years_ahead,years_ahead,HP_output,pred_quants){

  quantile = NULL
  if (!is(HP_output,'DemogMCMC')) stop("This function expects a 'DemogMCMC' object.")


  len_ages = length(ages)
  priorA<-c(log(0.25*0.0001/(1-0.25*0.0001)),log(3*0.02/(1-3*0.02)))
  priorB<-c(log(0.25*0.0001/(1-0.25*0.0001)),log(3*0.15/(1-3*0.15)))
  priorC<-c(log(0.25*0.01/(1-0.25*0.01)),log(3*0.3/(1-3*0.3)))
  priorD<-c(log(0.25*0.00005/(1-0.25*0.00005)),log(3*0.01/(1-3*0.01)))
  priorE<-c(log(0.00000000001),log(3*20))
  priorG<-c(log(0.25*0.0000001/(1-0.25*0.0000001)),log(3*0.01/(1-3*0.01)))
  priorH<-c(log(0.25*1),log(1*1.2))
  priorF<-prior_quant(log((15-10)/(40-15)),log((40-10)/(40-33)))


  cov_mat = HP_output$S
  thetas = HP_output$theta

  num_mcmc_samples = dim(thetas)[1]

  pHP= array(NA,c(num_years_ahead,len_ages,num_mcmc_samples))
  wnew<-matrix(NA,num_years_ahead,num_mcmc_samples)
  SS = array(NA,c(8,8,num_mcmc_samples))
  xnew = array(NA,c(num_years_ahead,8,num_mcmc_samples))

  x = array(NA,c(num_years_fitted,8,num_mcmc_samples))
  mu = HP_output$mu
  for(i in 1:num_mcmc_samples){
    SS[,,i] = matrix(cov_mat[i,],8,8)
    x[,,i] =  matrix(thetas[i,],num_years_fitted,8,byrow=F)

    for(t in 1:num_years_ahead){
      if(t==1){
        pred = rtmvnorm(1,mu= (x[num_years_fitted,,i] + mu[i,]), sigma=SS[,,i],lb=c(priorA[1],priorB[1],priorC[1],priorD[1],priorE[1],-Inf,priorG[1],priorH[1]), ub=c(priorA[2],priorB[2],priorC[2],priorD[2],priorE[2],priorF[2],priorG[2],priorH[2]))
      }else{
        pred = rtmvnorm(1,mu= (xnew[t-1,,i] + mu[i,]), sigma=SS[,,i],lb=c(priorA[1],priorB[1],priorC[1],priorD[1],priorE[1],-Inf,priorG[1],priorH[1]), ub=c(priorA[2],priorB[2],priorC[2],priorD[2],priorE[2],priorF[2],priorG[2],priorH[2]))
      }

      xnew[t,,i] =   pred
      pHP[t,,i] = prob(xnew[t,,i],ages)
    }
  }


  p.quant_up = apply(log(pHP),c(1,2),quantile,pred_quants[2])
  p.quant_low= apply(log(pHP),c(1,2),quantile,pred_quants[1])
  p.mean = apply(log(pHP),c(1,2),mean)

  rownames(pHP) <- as.character(years_ahead)
  colnames(pHP) <- as.character(ages)
  res = list(pred_quant_up = p.quant_up,pred_quant_low =p.quant_low,pred_mean =p.mean,pred_dist=log(pHP))
  class(res) <- "Pred"
  res
  }






































