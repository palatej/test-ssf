library(rjdssf)

ar_emp<-c(0.6710,0.1590)
sigma2u_e<-0.000004
ar_emp_act<-c(0.6687,0.1762)

ar_u<-c(0.4281,0.1880)
sigma2u_u<-0.000348
#ar_u<-c(0.4522,0.2239)
#sigma2u_u<-0.007645
#ar_u<-c(0.4239,0.1739)
#sigma2u_u<-0.001405

un<-read.csv("./Data/abs_un.csv")
un_sig<-read.csv("./Data/abs_un_sig2.csv")
emp<-read.csv("./Data/abs_emp.csv")
emp_sig<-read.csv("./Data/abs_emp_sig2.csv")

# y=series, ar=ar coefficients, uvar=u variance (either fixed or time varying)
# seasonal =seasonal model, noise component in the structural model
bsm_sae<-function(y, ar, uvar, seasonal="Trigonometric", noise=T){
  # create the model
  bsm<-rjdssf::model()
  
  # create the components and add them to the model
  # trend component
  rjdssf::add(bsm, rjdssf::locallineartrend("ll"))
  
  # seasonal component. Several specifcations available
  rjdssf::add(bsm, rjdssf::seasonal("s", 12, type=seasonal))
  
  #noise
  if (noise)
    rjdssf::add(bsm, rjdssf::noise("n"))
  
  #sample error
  if (length(uvar) == 1){
    rjdssf::add(bsm, rjdssf::ar("u", ar, fixedar = T, variance = 1, fixedvariance = T))
  }else{
    rjdssf::add(bsm, rjdssf::ar("u", ar, fixedar = T, variance = 1, fixedvariance = T))
  }
  # create the equation (no measurement error)
  eq<-rjdssf::equation("eq")
  
  rjdssf::add(eq, "ll")
  rjdssf::add(eq, "s")
  if (noise){
    rjdssf::add(eq, "n")
  }
  if (length(uvar) == 1){
    rjdssf::add(eq, "u", sqrt(uvar), T)
  }else{
    rjdssf::add(eq, "u", loading=rjdssf::varloading(0, sqrt(uvar)))
  }
  rjdssf::add(bsm, eq)
  
  return (rjdssf::estimate(bsm, y, optimizer="BFGS", marginal=T, concentrated=F))
}

printsae<-function(rslt, uvar, filtering=T){
  p<-result(rslt, "parameters")
  ll<-result(rslt, "likelihood.ll")
  cat("\nloglikelihood =", ll, "\n\n")
  cat("variances:\n")
  cat("level: ", p[1], "\n")
  cat("slope: ", p[2], "\n")
  cat("seasonal: ", p[3], "\n")
  cat("noise: ", p[4], "\n")
  cmp<-result(rslt, "ssf.cmppos")
  if (filtering){
    code<-paste("ssf.filtering.cmp(", length(cmp)-1, ")", sep="")
  }else{
    code<-paste("ssf.smoothing.cmp(", length(cmp)-1, ")", sep="")
  }
  sae<-sqrt(uvar)*result(rslt, code)
  return (sae)
}

plotsae<-function(rslt, uvar, filtering=T){
  cmp<-result(rslt, "ssf.cmppos")
  if (filtering){
    code<-paste("ssf.filtered.cmp(", length(cmp)-1, ")", sep="")
  }else{
    code<-paste("ssf.smoothing.cmp(", length(cmp)-1, ")", sep="")
  }
  sae<-exp(sqrt(uvar)*result(rslt, code))
  plot(sae, type="l")
}




# Survey variances (corrected for ar coefficients)
# if you want to you uncorrected survey variances, you should you rjdssf::sae instead of rjdssf::ar in the function bsm_sae
# (the correction is then automatically computed, which could allow direct (poor) estimation of the AR coefficients)
uvar<-un_sig[,13]
s<-log(un[,163])
rslt<-bsm_sae(s, ar_u, sigma2u_u)

zu<-printsae(rslt, sigma2u_u)
plotsae(rslt, uvar)

# Survey variances (corrected for ar coefficients)
uvar<-emp_sig[,4]
s<-log(emp[,118])
rslt<-bsm_sae(s, ar_emp_act, uvar)

ze<-printsae(rslt, uvar)
plotsae(rslt, uvar)
