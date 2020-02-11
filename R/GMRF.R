#'@title GMRF
#'@description MCMC algorithm to conduct Bayesian inference for the parameters and the latent states of the Gaussian Markov random field (GMRF) model
#'@param deaths  (number of ages) times (number of years) matrix with deaths
#'@param exposures  (number of ages) times (number of years) matrix with exposures to risk
#'@param its  number of MCMC iterations.
#'@param thin  choose desired thinning for the MCMC samples.
#'@param burnin  choose burn-in period for the sampler.
#'@param adapt  number of iterations within burn-in to learn step-size of MH
#'@param monitor  display the (monitor)th iteration
#'@param ages  vector with ages used to fit the model
#'@param years  vector with years used to fit the model
#'
#'@details
#'For details concerning the model and the algorithm please see the paper by \insertCite{alexopoulos2019bayesian;textual}{BayesDemog}
#'
#'
#'
#'@return The value returned is a list object of class \code{DemogMCMC} holding
#'\itemize{
#'\item{\code{alpha}}  {a vector with the thinned MCMC samples from the posterior of the parameter \code{alpha} of the model}
#'\item{\code{tau2}}  {a vector with the thinned MCMC samples from the posterior of the parameter \code{tau^2} of the model}
#'\item{\code{b}}  {a vector with the thinned MCMC samples from the posterior of the parameter \code{b} of the model}
#'\item{\code{z}}  {a matrix with the thinned MCMC samples from the posterior of the latent states \code{z} of the model}
#'\item{\code{mortality rates}}  {a 3D array with thinned MCMC samples from the posterior of the mortality rates}
#'\item{\code{accept_rate_z}}  {the acceptance rate of the Metropolis-Hastings step used to draw the latent states \code{z}}
#'\item{\code{accept_rate_alpha}}  {the acceptance rate of the Metropolis-Hastings step used to draw the parameter \code{alpha}}
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
#'
#'
#'
#'#Example 2:
#'\dontrun{
#'ages = c(0,1:89)
#'years = 1991:2000
#'country = "GBRTENW"
#'#please provide your HMD user otherwise do not continue
#'HMDuser = "user@..."
#'#please provide the password that corresponds to HMDuser otherwise do not continue
#'HMDpass = ""
#'data_label = "England & Wales Total Population"
#'data_series = "female"
#'
#'mydata =load_demog_data(ages ,years,country,HMDuser, HMDpass,data_label,data_series,TRUE,FALSE)
#'
#'its=15000
#'thin=10
#'burnin=5000
#'adapt=3000
#'monitor=1000
#'
#'mcmc = GMRF(mydata$deaths,mydata$exposures,its,thin,burnin,adapt,monitor,ages,years)
#'}
#'
#'
#'
#'@importFrom stats rgamma
#'@importFrom stats runif
#'@importFrom stats dnorm
#'@importFrom stats rnorm
#'@importFrom stats dbinom
#'@importFrom spam as.spam
#'@importFrom spam `diag.spam<-`
#'@importFrom spam kronecker.spam
#'@importFrom spam diag.spam
#'@importFrom spam chol.spam
#'@importFrom spam forwardsolve.spam
#'@importFrom spam backsolve.spam
#'
#'
#'
#'@export
#'
#'
#'
#'
GMRF<-function(deaths,exposures,its,thin,burnin,adapt,monitor,ages=c(0,1:89),years=1991:2000){

  y = as.vector(deaths)
  n = as.vector(exposures)
  yrs<-length(years)
  len_ages = length(ages)

  p<-y/n
  logitp<-log(p/(1-p))

  ymat = round(matrix(y,yrs,len_ages,byrow=T))
  nmat = round(matrix(n,yrs,len_ages,byrow=T))
  t<-rep(1:yrs,rep(length(ages),yrs))

  number<-(its-burnin)/thin
  indsave  = 0

  accrho<-0
  acc2 = 0

  accY =rep(0,burnin)
  alphasave=bsave=tausave =rep(NA,number)
  zsave=matrix(NA,yrs*len_ages,number)

  D = rep(1:yrs,len_ages)
  alpha = 0.5
  b = 0
  z = logitp


  eigs = eigsprop = matrix(NA,yrs,len_ages)

  RT =  tridiag( upper = rep(-1,yrs-1),lower =  rep(-1,yrs-1),main = c(1,rep(2,yrs-2),1))
  RT = as.spam(RT)
  Rn =  tridiag(upper =  rep(-1,len_ages-1), lower = rep(-1,len_ages-1),main= c(1,rep(2,len_ages-2),1))
  Rn = as.spam(Rn)
  srw = 0.001

  startings = z
  rate_save = array(NA,c(length(ages),yrs,(its-burnin)/thin))
  #start mcmc iterations
  for(i in 1:its){

    Q = alpha*kronecker.spam(Rn,diag(rep(1,yrs))) + (2-alpha)*kronecker.spam(diag(rep(1,len_ages)),RT)
    tau = rgamma(1, 0.001 + 0.5*len_ages*yrs, rate=0.001 + 0.5*(t(z-D*b)%*%Q%*%(z-D*b))     )

    b = rnorm(1,(solve(tau*(t(D)%*%Q%*%D)+0.001))%*%t(D)%*%(tau*Q)%*%z, solve(tau*(t(D)%*%Q%*%D)+0.001)) #gibbs step for the beta's

    alphapropN =log(alpha/(2-alpha)) + rnorm(1,0,0.5)
    alphaprop = 2*exp(alphapropN)/(1+exp(alphapropN))

    Qprop = alphaprop*kronecker.spam(Rn,diag(rep(1,yrs))) + (2-alphaprop)*kronecker.spam(diag(rep(1,len_ages)),RT)

    logratioa =dnorm(alphapropN,0,sqrt(100),log=T)-dnorm(log(alpha),0,sqrt(100),log=T)+ 0.5*logdet(alphaprop,yrs,len_ages) - 0.5*logdet(alpha,yrs,len_ages) -0.5*tau*(t(z-D*b)%*%Qprop%*%(z-D*b)) +0.5*tau*(t(z-D*b)%*%Q%*%(z-D*b))
    if(log(runif(1))<logratioa){
      alpha=alphaprop
      accrho = accrho+1
    }

    Q = alpha*kronecker.spam(Rn,diag(rep(1,yrs))) + (2-alpha)*kronecker.spam(diag(rep(1,len_ages)),RT)

    derivative = (  as.vector(ymat)*((exp(-z))^(1))  )/( 1 + exp(-z)   ) -(as.vector(nmat-ymat)*exp(-z)*( ((1+exp(-z))^(-2))   ))/(1-  ((1+exp(-z))^(-1))  )

    zz1 = rnorm(yrs*len_ages,0,sqrt(0.5*srw)) + z+ 0.5*srw*derivative

    cc = 2*zz1/srw + (tau*Q)%*%(D*b)

    QQ=tau*Q

    diag.spam(QQ) = diag.spam(QQ) + 2/srw
    U = chol.spam(QQ) # assuming Q is already a spam object
    bb = rnorm(nrow(QQ))
    meanMVN = forwardsolve.spam(U, cc)
    zprop = backsolve.spam(U, meanMVN + bb)

    derivprop= (  as.vector(ymat)*((exp(-zprop))^(1))  )/( 1 + exp(-zprop)   ) -(as.vector(nmat-ymat)*exp(-zprop)*( ((1+exp(-zprop))^(-2))   ))/(1-  ((1+exp(-zprop))^(-1))  )

    logratioz =  sum(( dbinom( as.vector(ymat),as.vector(nmat), exp(zprop)/(1+exp(zprop)),log=T)  )   )-sum(( dbinom( as.vector(ymat),as.vector(nmat), exp(z)/(1+exp(z)),log=T)  )   ) +sum( dnorm(zz1,zprop+0.5*srw*derivprop,sqrt(rep(srw/2,yrs*len_ages) ),log=T  )    )-sum( dnorm(zz1,z+0.5*srw*derivative,sqrt(rep(srw/2,yrs*len_ages) ),log=T  )    )+sum( dnorm(zz1,z,sqrt(rep(srw/2,yrs*len_ages) ),log=T  )    )-sum( dnorm(zz1,zprop,sqrt(rep(srw/2,yrs*len_ages) ),log=T  )    )


    if( log(runif(1))<logratioz   )
    {
      z = zprop
      acc2<-acc2+1
      if(i<burnin)accY[i]<-1
    }


    if(i>20 & i<adapt){
      if (sum( accY[(i-19):i]   ) < 8){
        srw<-srw-0.15*srw
      }
      if (sum( accY[(i-19):i]   ) > 12){
        srw<-srw+0.15*srw
      }
    }
    mortality_rate= matrix(exp(z)/(1+exp(z)),yrs,len_ages,byrow=F)

    if(i>burnin){
      if(i%%thin==0){
        indsave<-indsave+1
        rate_save[,,indsave] = t(mortality_rate)
        alphasave[indsave]= alpha
        tausave[indsave] = tau
        bsave[indsave] = b
        zsave[,indsave] = z

      }
    }
    if ((i%%monitor)==0) cat(paste("iteration",i),"\n")
  }

  res = list(alpha =alphasave,tau2 = tausave,b =bsave,z = t(zsave),mortality_rates = rate_save,accept_rate_z=acc2/its,accept_rate_alpha=accrho/its)
  class(res) <- "DemogMCMC"
  return(res)
  }







