#'@title GMRFpred
#'@description returns the predictive distribution of the mortality rates obtained by fitting the Gaussian Markov random field model
#'@param ages numeric vector with the ages used to fit the model, default to the vector c(0,1:89)
#'@param num_years_fitted number of years used to fit the model, default is 10
#'@param num_years_ahead  number of years to forecast the mortality rates for each age in ages
#'@param years_ahead years in which predictions are required
#'@param GMRF_output  output of the MCMC sampler for the Gaussian Markov random field model
#'@param pred_quants a numeric vector with two quantiles that will form the credible interval of the predictive distributions
#'
#'@details
#'For details concerning the method used to draw samples from the predictive distribution of the mortality rates please see the paper by \insertCite{alexopoulos2019bayesian;textual}{BayesDemog}
#'
#'
#'
#'@return The value returned is a list object of class \code{Pred} holding
#'\itemize{
#'\item{\code{pred_quant_up}}  {a matrix consisted of the upper quantile of the prediction credible interval for the log-mortality rate of each age in the input vector \code{ages}}
#'\item{\code{pred_quant_low}}  {a matrix consisted of the lower quantile of the prediction credible interval for the log-mortality rate of each age in the input vector \code{ages}}
#'\item{\code{pred_mean}}  {a matrix consisted of the mean of the posterior predictive distribution of the log-mortality rate for each age in the input vector \code{ages}}
#'\item{\code{pred_dist}}  {a matrix with samples drawn from the predictive distributions of the mortality rates for each age in the input vector \code{ages}}
#'}
#'
#'@references
#'  \insertAllCited{}
#'
#'
#'@examples
#'#Example 1:
#'
#'data(UKandWales)
#'
#'its=15000
#'thin=10
#'burnin=5000
#'adapt=3000
#'monitor=1000
#'
#'years_fit = 1991:2000
#'myD = mydata$deaths[,as.character(years_fit)]
#'myE = mydata$exposures[,as.character(years_fit)]
#'
#'mcmc = GMRF(myD,myE,its,thin,burnin,adapt,monitor,ages,years_fit)
#'GMRFforecast = GMRFpred(c(0,1:89),10,10,2001:2010,GMRFmcmc,c(0.025,0.975))
#'
#'
#'#Example 2:
#'
#'\dontrun{
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
#'mydata =load_demog_data(ages,years,country,HMDuser,HMDpass,data_label,data_series,TRUE,FALSE)
#'
#'its=15000
#'thin=10
#'burn=5000
#'adapt=3000
#'monitor=1000
#'
#'years_fit = 1991:2000
#'myD = mydata$deaths[,as.character(years_fit)]
#'myE = mydata$exposures[,as.character(years_fit)]
#'
#'GMRFmcmc = GMRF(myD,myE,its,thin,burn,adapt,monitor,ages,years_fit)
#'GMRFforecast = GMRFpred(c(0,1:89),10,10,2001:2010,GMRFmcmc,c(0.025,0.975))
#'}
#'
#'@importFrom methods is
#'@importFrom spam as.spam
#'@importFrom spam solve.spam
#'@importFrom spam as.dgCMatrix.spam
#'@importFrom sparseMVN rmvn.sparse
#'@importFrom Matrix Cholesky
#'
#'
#'@export
#'
#'
#'
#'
GMRFpred<-function(ages,num_years_fitted,num_years_ahead,years_ahead,GMRF_output,pred_quants){

  quantile = NULL

  if (!is(GMRF_output,'DemogMCMC')) stop("This function expects a 'DemogMCMC' object.")


  alpha = GMRF_output$alpha
  b = GMRF_output$b
  tau=GMRF_output$tau2
  x=t(GMRF_output$z)
  len_ages = length(ages)
  num_mcmc_samples = length(alpha)

  yrs = num_years_ahead
  newyrs = yrs+num_years_fitted
  RT =  tridiag( upper = rep(-1,newyrs-1),lower =  rep(-1,newyrs-1),main = c(1,rep(2,newyrs-2),1))
  Rn =  tridiag(upper =  rep(-1,len_ages-1), lower = rep(-1,len_ages-1),main= c(1,rep(2,len_ages-2),1))
  xnew = matrix(NA,num_years_ahead*len_ages,num_mcmc_samples)

  for(i in 1:num_mcmc_samples){
    QQ = tau[i]*alpha[i]*kronecker(Rn,diag(rep(1,newyrs))) + tau[i]*(2-alpha[i])*kronecker(diag(rep(1,len_ages)),RT)



    indsA = c(1:num_years_fitted,(num_years_fitted+num_years_ahead+1):(num_years_fitted+num_years_ahead+num_years_fitted));TT=num_years_fitted
    for(ii in 3:len_ages){
      TT=TT+num_years_fitted+num_years_ahead
      indsA = append(indsA,(TT+(num_years_ahead+1)):(TT+(num_years_fitted+num_years_ahead)))
    }
    indsB = c(1:((num_years_fitted+num_years_ahead)*len_ages))[-indsA]


    QBB = QQ[indsB,indsB]
    QBA = QQ[indsB,indsA]

    rm(QQ)
    gc()

    bb = rep(b[i]*c(1:num_years_fitted),len_ages)
    bnew = rep(b[i]*((num_years_fitted+1):(num_years_fitted+num_years_ahead)),len_ages)

    spamQBB<-as.spam( QBB)
    QBBinv<-solve.spam(spamQBB)
    M<-QBBinv%*%as.spam(QBA)%*%(x[,i]-bb)
    xnew[,i] = bnew -M[,1]+ rmvn.sparse(1, rep(0,num_years_ahead*len_ages), Cholesky(as.dgCMatrix.spam(spamQBB)), prec = TRUE)

  }

  p= array(NA,c(num_years_ahead,len_ages,num_mcmc_samples))
  for(i in 1:num_mcmc_samples){
    p[,,i]= matrix(exp(xnew[,i])/(1+exp(xnew[,i])),yrs,len_ages,byrow=F)
  }

  p.quant_up = apply(log(p),c(1,2),quantile,pred_quants[2])
  p.quant_low= apply(log(p),c(1,2),quantile,pred_quants[1])
  p.mean = apply(log(p),c(1,2),mean)

  rownames(p) <- as.character(years_ahead)
  colnames(p) <- as.character(ages)
  res = list(pred_quant_up = p.quant_up,pred_quant_low =p.quant_low,pred_mean =p.mean,pred_dist=log(p))
  class(res) <- "Pred"
  res
  }


