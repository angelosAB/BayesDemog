
logdet <- function(a,yrs,n){

  ll = matrix(NA,n,yrs)

  for(i in 1:n)for(j in 1:yrs) ll[i,j] = 2*(  a*(   1- cos( pi*(i-1)/n   ))  +(2-a)*(   1- cos( pi*(j-1)/yrs   ))          )

  lll<-log(ll[which(ll!=0)])
  ld = sum(lll)
  return(ld)
}
