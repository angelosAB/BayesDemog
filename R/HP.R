#'@title HP
#'@description MCMC algorithm to conduct Bayesian inference for the parameters and the latent states of the dynamic Heligman-Pollard model
#'@param deaths  (number of ages) times (number of years) matrix with deaths
#'@param exposures  (number of ages) times (number of years) matrix with exposures to risk
#'@param its  number of MCMC iterations.
#'@param thin  choose desired thinning for the MCMC samples.
#'@param burn  choose burn-in period for the MCMC samples.
#'@param adapt  number of iterations within burn-in to learn step-size of MH
#'@param monitor  display the (monitor)th iteration
#'@param df  degrees of freedom for the IW prior for the covariance of the random walk.
#'@param A scale parameter for the prior on the covariance of the random walk
#'@param ages  vector with ages used to fit the model
#'@param years  vector with years used to fit the model
#'
#'@details
#'For details concerning the model and the algorithm please see the paper by \insertCite{alexopoulos2019bayesian;textual}{BayesDemog}
#'
#'
#'
#'
#'@return The value returned is a list object of class \code{DemogMCMC} holding
#'\itemize{
#'\item{\code{S}}  {a matrix with the thinned MCMC samples from the posterior of the parameter \code{Sigma} of the model}
#'\item{\code{theta}}  {a matrix with the thinned MCMC samples from the posterior of the latent states \code{theta} of the model}
#'\item{\code{mu}}  {a matrix with the thinned MCMC samples from the posterior of the parameter \code{mu} of the model}
#'\item{\code{mortality rates}}  {a 3D array with thinned MCMC samples from the posterior of the mortality rates.}
#'\item{\code{accept_rate_theta}}  {the acceptance rate of the Metropolis-Hastings step used to draw the latent states \code{theta}}
#'}
#'
#'
#'@references
#'  \insertAllCited{}
#'
#'@examples
#'
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
#'
#'#Example 2:
#'\dontrun{
#'ages = c(0,1:89)
#'years = 1991:2000
#'country = "GBRTENW"
#'HMDuser = "user@..." #please provide your HMD user otherwise do not continue
#'HMDpass = ""  #please provide the password that corresponds to HMDuser otherwise do not continue
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
#'mcmc = HP(myD,myE,its,thin,burn,adapt,monitor,df=2,A=10^10,ages,years)
#'}
#'
#'
#'@importFrom stats rgamma
#'@importFrom stats runif
#'@importFrom stats nls
#'@importFrom stats coef
#'@importFrom stats vcov
#'@importFrom MASS mvrnorm
#'@importFrom mvtnorm dmvnorm
#'@importFrom MCMCpack riwish
#'
#'
#'
#'@export
#'
#'
#'
#'
HP<-function(deaths,exposures,its=5*10^5,thin=200,burn=10^5,adapt,monitor,df=2,A=10^10,ages=c(0,1:89),years=1991:2000){


  y = as.vector(deaths)
  n = as.vector(exposures)
  x<-ages
  yrs<-length(years)

  #truncation point of the random walk
  priorA<-c(log(0.25*0.0001/(1-0.25*0.0001)),log(3*0.02/(1-3*0.02)))
  priorB<-c(log(0.25*0.0001/(1-0.25*0.0001)),log(3*0.15/(1-3*0.15)))
  priorC<-c(log(0.25*0.01/(1-0.25*0.01)),log(3*0.3/(1-3*0.3)))
  priorD<-c(log(0.25*0.00005/(1-0.25*0.00005)),log(3*0.01/(1-3*0.01)))
  priorE<-c(log(0.00000000001),log(3*20))
  priorG<-c(log(0.25*0.0000001/(1-0.25*0.0000001)),log(3*0.01/(1-3*0.01)))
  priorH<-c(log(0.25*1),log(1*1.2))
  priorF<-prior_quant(log((15-10)/(40-15)),log((40-10)/(40-33)))


  propvar<-array(NA,c(8,8,yrs))
  starthp<-matrix(0,yrs,8)
  j=0
  t<-rep(1:yrs,rep(length(ages),yrs))
  for(jj in 1:yrs){
    j=j+1
    ps<-y[t==jj]/n[t==jj]
    opt<-nls(ps~((exp(a)/(1+exp(a)))^((x+(exp(b)/(1+exp(b))))^((exp(c)/(1+exp(c)))))+(exp(d)/(1+exp(d)))*exp((-exp(e))*(log(x)-log((40*exp(f)+10)/(1+exp(f))))^2)+(exp(g)/(1+exp(g)))*(exp(h*x)))/(1+((exp(a)/(1+exp(a)))^((x+(exp(b)/(1+exp(b))))^((exp(c)/(1+exp(c)))))+(exp(d)/(1+exp(d)))*exp((-exp(e))*(log(x)-log((40*exp(f)+10)/(1+exp(f))))^2)+(exp(g)/(1+exp(g)))*(exp(h*x)))),start=list(a=-7,b=-4.5,c=-2.2,d=-9,e=2.5,f=log((18-10)/(40-18)),g=-10.5,h=0.098),weights=1/(ps^2),control=list(maxiter=400,warnOnly = T))
    starthp[j,]<-coef(opt)
    propvar[,,j]<-vcov(opt)
  }
  startmu=rep(0,8)

  muprior<-diag(0.001,8)
  Pprior<-diag(0.0001,8)
  a<-matrix(NA,(its-burn)/thin,8)
  aa<-1/rgamma(8,0.5,scale=A)

  S<-riwish(  8+1,  diag(rep(1,8))  )
  jj<-0
  tun<-rep(0.1,yrs)
  accY =matrix(0,burn,yrs)
  musave<-matrix(0,(its-burn)/thin,8)
  Ssave<-matrix(0,(its-burn)/thin,64)
  thsave<-matrix(0,(its-burn)/thin,yrs*8)
  mu<-startmu
  th<-starthp
  acc<-rep(0,yrs)
  P<-solve(S)
  accept = rep(0,yrs)
  mortality_rate = matrix(NA,length(ages),yrs)
  rate_save = array(NA,c(length(ages),yrs,(its-burn)/thin))
  colnames(rate_save) <- as.character(years)
  rownames(rate_save) <- as.character(ages)
  #start mcmc iterations
  for (i in 1:its){
    aa<-1/rgamma(8,(df+8)/2,scale=1/((1/A)+df*diag(P)))
    diff<-th[2:yrs,]-th[1:(yrs-1),]
    V<-2*df*diag(1/aa)

    for (j in 1:(yrs-1)) V<-V+outer( diff[j,]-mu,diff[j,]-mu)

    S<-riwish(yrs+df+8-2,V)

    P<-solve(S)
    mu<-mvrnorm(1,solve(muprior+(yrs-1)*P)%*%P%*%(th[yrs,]-th[1,]),solve(muprior+(yrs-1)*P))


    for (k in 1:yrs){

      thprop<-th[k,]+mvrnorm(1,rep(0,8),tun[k]*propvar[,,k])
      if( (thprop[1]<=priorA[1]) || (thprop[1]>=priorA[2]) ||
          (thprop[2]<=priorB[1]) || (thprop[2]>=priorB[2]) ||
          (thprop[3]<=priorC[1]) || (thprop[3]>=priorC[2]) ||
          (thprop[4]<=priorD[1]) || (thprop[4]>=priorD[2]) ||
          (thprop[5]<=priorE[1]) || (thprop[5]>=priorE[2]) ||
          (thprop[6]<=-Inf) || (thprop[6]>=2.64) ||
          (thprop[7]<=priorG[1]) || (thprop[7]>=priorG[2]) ||
          (thprop[8]<=priorH[1]) || (thprop[8]>=priorH[2]) )
      {
        alpha<- -Inf
      }
      else{
        alpha<-HPloglik(thprop,y[t==k],n[t==k],length(ages),yrs,ages)-HPloglik(th[k,],y[t==k],n[t==k],length(ages),yrs,ages)

        if (k==1) {alpha<-alpha+dmvnorm(thprop,th[2,]-mu,S,log=T)-dmvnorm(th[k,],th[2,]-mu,S,log=T)
        }
        else {if (k==yrs)
        {alpha<-alpha+dmvnorm(thprop,th[yrs-1,]+mu,S,log=T)-dmvnorm(th[k,],th[yrs-1,]+mu,S,log=T)}
          else
          {alpha<-alpha+dmvnorm(thprop,(th[k-1,]+th[k+1,])/2,S/2,log=T)-dmvnorm(th[k,],(th[k-1,]+th[k+1,])/2,S/2,log=T)}}
      }


      if (log(runif(1))<alpha){
        th[k,]<-thprop
        accept[k]<-accept[k]+1
        if(i<=burn){
          accY[i,k]<-1
        }
      }

      mortality_rate[,k] = prob(th[k,],x)

      if(i>20 & i<adapt){
        if (sum( accY[(i-19):i,k]   ) < 6){
          tun[k]<-tun[k]-0.15*tun[k]
        }
        if (sum( accY[(i-19):i,k]   ) > 9){
          tun[k]<-tun[k]+0.15*tun[k]
        }
      }
    }

    if(i>burn){
      if(i%%thin==0){
        jj<-jj+1
        rate_save[,,jj] = mortality_rate
        Ssave[jj,]<-as.vector(S)
        musave[jj,]<-mu
        thsave[jj,]<-as.vector(th)
        a[jj,]<-aa
      }
    }
    if ((i%%monitor)==0) cat(paste("iteration",i),"\n")
  }
  res = list(S=Ssave,theta=thsave,mu=musave,mortality_rates=rate_save,accept_rate_theta=accept/its)
  class(res) <- "DemogMCMC"
  return(res)
}


