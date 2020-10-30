library(rjdssf)

load("./Data/vat.rda")

# create the model
bsm_td<-function(s, tdgroups=NULL, nweights=NULL){
  bsm<-rjdssf::model()
  freq<-frequency(s)
  start<-start(s)
  
  # create the components and add them to the model
  rjdssf::add(bsm, rjdssf::locallineartrend("ll"))
  rjdssf::add(bsm, rjdssf::seasonal("s", freq, type="HarrisonStevens"))
    rjdssf::add(bsm, rjdssf::noise("n"))
  if (! is.null(tdgroups)){
    rjdssf::add(bsm, rjdssf::td("td", freq, start, length(s), tdgroups, T, variance = 0, fixed=T))
  }
  # create the equation (fix the variance of the noise to 1)
  eq<-rjdssf::equation("eq")
  rjdssf::add(eq, "ll")
  rjdssf::add(eq, "s")
  if (is.null(nweights)){
    rjdssf::add(eq, "n")
  }else{             
    rjdssf::add(eq, "n", 1, T, rjdssf::varloading(0, nweights))
  }
  if (! is.null(tdgroups)){
    rjdssf::add(eq, "td")
  }
  rjdssf::add(bsm, eq)
  #estimate the model
  rslt<-rjdssf::estimate(bsm, s, initialization="SqrtDiffuse")
  
  # variances
  p<-result(rslt, "parameters")*result(rslt, "scalingfactor")
  model<-list(levelVariance=p[1],
              slopeVariance=p[2],
              seasVariance=p[3],
              noiseVariance=p[4])
  ss<-result(rslt, "ssf.smoothing.states")
  n<-ss[,2+freq]
  if (! is.null(nweights)){

    n<-n*nweights
  }
  seas<-ts(ss[,3], frequency = freq, start=start)
  decomposition<-list(
    series=s,
    sa=s-seas,
    level=ts(ss[,1], frequency = freq, start=start),
    slope=ts(ss[,2], frequency = freq, start=start),
    seas=seas,
    noise=ts(n, frequency = freq, start=start)
  )
  return(structure(list(
    model=model,
    decomposition=decomposition),
    class="JDBSM"))
}

plot.JDBSM<-function(q){
  par(mfrow=c(2,1))
  ts.plot(ts.union(q$decomposition$series, q$decomposition$level, q$decomposition$sa),
          col=c("gray", "red", "blue"))
  ts.plot(ts.union(q$decomposition$seas, q$decomposition$noise),
          col=c("magenta", "green"))
  par(mfrow=c(1,1))
}

test<-function(s, o){
  n<-length(s)
  w<-array(1, dim=n)
  w[(n-length(0)):n]<-o
  return (bsm_td(s, nweights = w))
}

s<-log(vat$I.H.bergement.et.restauration)

q<-test(s, c(1,1))
plot(q)
a<-test(s, c(2,2))
plot(a)
b<-test(s, c(3,3))
plot(b)
c<-test(s, c(5,5))
plot(c)
d<-test(s, c(10,10))
plot(d)

