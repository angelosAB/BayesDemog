#'@title computeESS
#'@description computes the effective sample size (ESS) of the
#'samples drawn from the joint posterior of parameters and
#'latent states of either the dynamic Heligman-Polard model or the Gaussian Markov random field model
#'
#'@param x a \code{DemogMCMC} object
#'@param model either "Dynamic_HP" for the dynamic Heligman-Polard model or "GMRF" for the Gaussian Markov random field model
#'
#'
#'@return The value returned is a list holding matrices with the effective sample sizes of the parameters and latent states of the model specified by the argument \code{model}. See also Details.
#'
#'
#'@details
#'In the case of the Dynamic Heligman-Pollard model the function returns a list consisted of the matrices
#'\code{SigmaESS}, \code{muESS} and \code{thetaESS}
#'with the ESS for the parameters \code{Sigma} and \code{mu} and the latent states \code{theta}.
#'In the case of the Gaussian Markov random field model the returned list contains the vector \code{parsESS} consisted of
#'the ESS for the parameters \code{alpha}, \code{tau}^2, \code{b} and the matrix zESS with the ESS for the latent states \code{z}.
#'Please see in \insertCite{alexopoulos2019bayesian;textual}{BayesDemog} for details.
#'
#'
#'
#'@examples
#'
#'\dontrun{
#'ages = c(0,1:89)
#'years =  1991:2000
#'country = "GBRTENW"
#'#please provide your HMD user otherwise do not continue
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
#'mcmc = HP(mydata$deaths,mydata$exposures,its,thin,burn,adapt,monitor,df=2,A=10^10,ages,years)
#'ESS = computeESS(mcmc,"Dynamic_HP")
#'}
#'
#'@references
#'  \insertAllCited{}
#'
#'@importFrom methods is
#'@importFrom coda effectiveSize
#'
#'
#'
#'
#'@export
#'
#'
#'

computeESS = function(x,model){

  if (!is(x, "DemogMCMC")) stop("This function expects a 'DemogMCMC' object.")

  if(model=="Dynamic_HP"){

    SigmaESS = effectiveSize(x$S)
    thetaESS = effectiveSize(x$theta)
    muESS = effectiveSize(x$mu)
    res = list(SigmaESS=SigmaESS,muESS=muESS,thetaESS=thetaESS)
  }
  if(model=="GMRF"){
    alphaESS =effectiveSize(x$alpha)
    tau2ESS =effectiveSize(x$tau2)
    bESS = effectiveSize(x$b)
    zESS = effectiveSize(x$z)
    res = list(parsESS =c(alphaESS,tau2ESS,bESS),zESS=zESS)
  }
  return(res)
}
