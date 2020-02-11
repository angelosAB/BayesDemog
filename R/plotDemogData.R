#'@title plotDemogData
#'@description plots raw mortality rates read from the Human Mortality Database
#'@param x a DemogData object
#'@param ages numeric vector with the ages
#'@param years numeric vector with years
#'@param log  indicates if the rates should be plotted in the log-scale, default is TRUE
#
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
#'mydata =load_demog_data(ages,years,country,HMDuser, HMDpass,data_label,data_series,TRUE,FALSE)
#'
#'#plot in log-scale
#'plotDemogData(mydata,c(0,1:89),c(1991,1995,2000),T)
#'#plot in the actual scale
#'plotDemogData(mydata,c(0,1:89),c(1991,1995,2000),F)
#'}
#'
#'
#'@importFrom methods is
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 aes
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 xlab
#'@importFrom ggplot2 ylab
#'@importFrom ggplot2 labs
#'
#'
#'@export
#'
#'
#'
#'
plotDemogData <- function(x,ages,years,log=TRUE) {

   Age = prop = Year = NULL

   if (!is(x, "DemogData")) stop("This function expects a 'DemogData' object.")

   myages = as.character(ages)
   myYears = as.character(years)
   jj=0
   df=NULL
   if(log==TRUE){
      ytitle ="Probability of death (log-scale)"
   }else{
      ytitle ="Probability of death"

   }
   for(j in myYears){
      jj=jj+1
      if(log==TRUE){
         y = log(x$deaths[myages,j]/x$exposures[myages,j])
      }else{
         y = x$deaths[myages,j]/x$exposures[myages,j]
      }
      df = rbind(df,data.frame(Age = ages,Year=rep(years[jj],length(ages)), prop=y ))
   }
   gg = ggplot(df,mapping = aes(x=Age,y=prop,color=as.factor(Year)))+geom_line()
   gg = gg+xlab("Age") +ylab(ytitle)
   gg = gg+ labs(color = "Year")
   gg
}



