#'@title load_demog_data
#'@description Read raw data from the website of the Human Mortality Database (HMD) and return matrices with number of deaths and exposures to risk for ages and years
#'specified by the user.
#'
#'
#'@param ages numeric vector with the ages to be analyzed, default is c(0,1:89).
#'@param years numeric vector with years to be used for the estimation of the models, default to 1991:2000
#'@param country directory abbreviation from the HMD. For instance, England & Wales Total Population = "GBRTENW". See help page of \code{hmd.mx} function in \code{demography} package for details.
#'@param HMDuser username of HMD account (case-sensitive).
#'@param HMDpass password of HMD account (case-sensitive).
#'@param data_label character string giving name of country from which the data are taken.
#'@param data_series character string giving the name of the extracted series, default is "female".
#'@param Central2Initial {logical, if \code{TRUE} we transform exposures to risk from central to initial}. See also \code{Details}
#'@param Initial2Central {logical, if \code{TRUE} we transform exposures to risk from initial to central}. See also \code{Details}
#'
#'@details
#'This function utilizes routines from the R packages \code{demography} \insertCite{demogR}{BayesDemog}  and \code{StMoMo} \insertCite{StMoMo}{BayesDemog}
#'in order to read raw data from the website of the Human Mortality Database and to convert them
#'in suitable form for Bayesian analysis by using the MCMC algorithms developed by \insertCite{alexopoulos2019bayesian;textual}{BayesDemog}.
#'By using the functions \code{central2initial} and \code{initial2central} from \code{StMoMo} package
#'we are able to transform central exposures to risk to initial and initial exposures to central respectively if needed.
#'
#'
#'@return The value returned is a list object of class \code{Pred} holding
#'
#'
#'\itemize{
#'\item{\code{deaths}}  {...}
#'\item{\code{exposures}}  {...}
#'}
#'
#'@examples
#'
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
#'}
#'
#'
#'@references
#'  \insertAllCited{}
#'
#'@importFrom demography hmd.mx
#'@importFrom StMoMo StMoMoData
#'@importFrom StMoMo central2initial
#'@importFrom StMoMo initial2central
#'
#'
#'
#'
#'
#'@export




load_demog_data = function(ages = c(0,1:89),years = 1991:2000,country ,
                           HMDuser, HMDpass,
                           data_label ,data_series,
                           Central2Initial=TRUE,Initial2Central = FALSE){
  ################read raw data from HMD################################
  mydata<-hmd.mx(country ,HMDuser, HMDpass , data_label)

  ################choose gender (univariate modelling)################################
  mydataStMoMo <- StMoMoData(mydata, data_series)

  ################make suitable transformations################################
  if(Central2Initial){
    data<-central2initial(mydataStMoMo)
    data$Dxt<-round(data$Dxt)
    data$Ext<-round(data$Ext)
  }
  if(Initial2Central){
    data<-initial2central(mydataStMoMo)
    data$Dxt<-round(data$Dxt)
    data$Ext<-round(data$Ext)
  }
  if((Central2Initial==FALSE) & (Initial2Central==FALSE) ){
    data<- mydataStMoMo
    data$Dxt<-round(data$Dxt)
    data$Ext<-round(data$Ext)
  }

  myages = as.character(ages)
  insample_years = as.character(years)

  myD =  data$Dxt[myages,insample_years]
  myE =  data$Ext[myages,insample_years]

  res <- list(deaths =myD,exposures = myE)

  class(res) <- "DemogData"
  res
}



