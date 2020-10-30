library(rjdssf)

#Sugar production in Belgium. Break in 1987 (pos=156)
data<-c(34.8,21,32.2,22.4,18.2,15.9,13.3,16.7,19,91.5,246.5,258.3,126.5,19.7,24.8,16.7,12.2,14.3,10.9,14.9,17.6,211.7,358.3,201,
        4,34.6,31.3,17,16.1,16.8,9.6,18.1,18.5,215.4,329.6,244.3,21.4,16.9,37.5,16.8,14.4,30,12,14.3,17.4,242.3,323.4,257.3,16,
        13.6,18.9,28.8,20.4,43.1,7.8,17.2,23.3,240.2,371.5,272.5,24.5,22.1,25.6,29.8,47.3,16.4,10,17.8,24.3,322.1,380,318.9,31,
        30.8,31.5,42.3,22.6,18.2,9.3,21,27.3,305.6,352.2,251.2,17.2,17.9,34.4,31.7,28.7,21.2,9,16.1,57.5,351.1,379.1,358.5,108.2,
        20,34.2,28.2,27.1,17.5,9.2,16.3,107.8,373.5,399.3,381,26.7,15.4,18.8,39.7,21.6,21.1,6.6,19.1,18.8,257.9,419.3,234,16.9,
        19.1,19.5,28.6,23.5,21.5,17.6,19,10,279.3,406.3,291.2,19,16.1,24.6,30.1,18.9,19.2,5,12.3,34.6,367.8,423.5,229.5,17,12.6,
        14.7,48.4,24.9,23.2,18.6,17.1,33.9,310.1,365.1,177.1,64.1,74.1,82.5,75.7,70.7,78.2,25.6,78.3,96.5,138.1,137.2,93.6,65.9,71.6,
        79.4,67.8,64,75.4,16.1,74.2,90.7,134.3,132.1,82,56.7,61,63.7,63.2,57.7,68.5,16.1,64,80.9,115.7,115.8,91.3,56,60.3,65.6,53.5,
        59.5,61.4,17.3,64.1,82.8,112.1,108.9,79.6,60.1,56.8,62.2,59.1,54.4,59.1,17.2,56.6,70.7,111.2,107.1,74.2,53.8,53.1,55.3,50.1,
        47.3,53.1,19.5,47.8,72.6,99.6,95.5,64.2,48.6,50.3,57.6,46.1,43.2,52.8,13.9,48.2,71.1,91.9,89.3,73,37.1,42.9,51.1,39.8,37,
        49.1,17.3,45.2,52.5,72.6,73.5,61.3)
sugar<-ts(data, freq=12, start=c(1974,1))

splitsa<-function(s, breakPos, seasonal=c("HarrisonStevens", "Crude", "Trigonometric", "Dummy"), common.trend=T, common.irregular=F, common.seasonal=F){
  freq=frequency(s)
  seasonal<-match.arg(seasonal)
  n<-length(s)
  m<-array(dim=c(n, 2))
  idx1<-1:breakPos
  idx2<-(breakPos+1):n
  m[idx1,1]<-s[idx1]
  m[idx2, 2]<-s[idx2]
  
  cmp1<-array(dim=3)
  cmp2<-array(dim=3)
  m1<-array(dim=4)
  m2<-array(dim=4)
  pos<-1
  mpos<-1
  # Create the model 
  model<-rjdssf::model()
  # create the equations (fix the variance to 0)
  eq1<-rjdssf::equation("eq1")
  eq2<-rjdssf::equation("eq2")
  # create the components and add them to the model
  if (common.trend){
    rjdssf::add(model, rjdssf::locallineartrend("l"))
    rjdssf::add(eq1, "l")
    rjdssf::add(eq2, "l")
    cmp1[1]<-1
    cmp2[1]<-1
    pos<-3
    m1[1]<-1
    m2[1]<-1
    m1[2]<-2
    m2[2]<-2
    mpos<-3
  }else{
    rjdssf::add(model, rjdssf::locallineartrend("l1"))
    rjdssf::add(eq1, "l1")
    rjdssf::add(model, rjdssf::locallineartrend("l2"))
    rjdssf::add(eq2, "l2")
    cmp1[1]<-1
    cmp2[1]<-3
    pos<-5
    m1[1]<-1
    m2[1]<-2
    m1[2]<-3
    m2[2]<-4
    mpos<-5
  }
  if (common.seasonal){
    rjdssf::add(model, rjdssf::seasonal("s", freq, type=seasonal))
    rjdssf::add(eq1, "s")
    rjdssf::add(eq2, "s")
    cmp1[2]<-pos
    cmp2[2]<-pos
    pos<-pos+freq-1
    m1[3]<-mpos
    m2[3]<-mpos
    mpos<-mpos+1
  }else{
    rjdssf::add(model, rjdssf::seasonal("s1", freq, type=seasonal))
    rjdssf::add(eq1, "s1")
    rjdssf::add(model, rjdssf::seasonal("s2", freq, type=seasonal))
    rjdssf::add(eq2, "s2")
    cmp1[2]<-pos
    cmp2[2]<-pos+freq-1
    pos<-pos+2*(freq-1)
    m1[3]<-mpos
    m2[3]<-mpos+1
    mpos<-mpos+2
  }
  if (common.irregular){
    rjdssf::add(model, rjdssf::noise("n"))
    rjdssf::add(eq1, "n")
    rjdssf::add(eq2, "n")
    cmp1[3]<-pos
    cmp2[3]<-pos
    m1[4]<-mpos
    m2[4]<-mpos
  }else{
    rjdssf::add(model, rjdssf::noise("n1"))
    rjdssf::add(eq1, "n1")
    rjdssf::add(model, rjdssf::noise("n2"))
    rjdssf::add(eq2, "n2")
    cmp1[3]<-pos
    cmp2[3]<-pos+1
    m1[4]<-mpos
    m2[4]<-mpos+1
  }
  
  rjdssf::add(model, eq1)
  rjdssf::add(model, eq2)
  rslt<-rjdssf::estimate(model, m, marginal=F, initialization="SqrtDiffuse", optimizer="LevenbergMarquardt", concentrated=T, precision=1e-15)
  
  p<-result(rslt, "parameters")
  factor<-result(rslt, "scalingfactor")
  names<-c("level", "slope", "seasonal", "noise")
  model<-data.frame(var1=c(factor*p[m1[1]],factor*p[m1[2]], factor*p[m1[3]], factor*p[m1[4]]), 
                    var2=c(factor*p[m2[1]],factor*p[m2[2]], factor*p[m2[3]], factor*p[m2[4]]),
                    row.names = names)

    ll<-result(rslt, "likelihood.ll")
  ser<-result(rslt, "likelihood.ser")
  res<-result(rslt, "likelihood.residuals")
  
  likelihood<-list(loglikelihood=ll, ser=ser, residuals=res)
  sm<-rjdssf::smoothedstates(rslt)
  start<-start(s)
  decomposition<-list(
    series=s)
  start=start(s)
  if (common.trend){
    level=ts(sm[,cmp1[1]], frequency = freq, start = start)
    decomposition[["level"]]<-level
  }else{
    level1=ts(sm[,cmp1[1]], frequency = freq, start = start)
    level2=ts(sm[,cmp2[1]], frequency = freq, start = start)
    decomposition[["level1"]]<-level1
    decomposition[["level2"]]<-level2
    level<-level1
    level[idx2]<-level2[idx2]
    decomposition[["level"]]<-level
  }
  if (common.seasonal){
    seas=ts(sm[,cmp1[2]], frequency = freq, start = start)
    decomposition[["seas"]]<-seas
    decomposition[["sa"]]<-s-seas
  }else{
    seas1=ts(sm[,cmp1[2]], frequency = freq, start = start)
    seas2=ts(sm[,cmp2[2]], frequency = freq, start = start)
    decomposition[["seas1"]]<-seas1
    decomposition[["seas2"]]<-seas2
    seas<-seas1
    seas[idx2]<-seas2[idx2]
    decomposition[["seas"]]<-seas
    decomposition[["sa"]]<-s-seas
  }
  if (common.irregular ){
    noise=ts(sm[,cmp1[3]], frequency = freq, start = start)
    decomposition[["noise"]]<-noise
  }else{
    noise1=ts(sm[,cmp1[3]], frequency = freq, start = start)
    noise2=ts(sm[,cmp2[3]], frequency = freq, start = start)
    decomposition[["noise1"]]<-noise1
    decomposition[["noise2"]]<-noise2
    noise<-noise1
    noise[idx2]<-noise2[idx2]
    decomposition[["noise"]]<-noise
  }
  
  return (structure(
    list(
         model=model,
         decomposition=decomposition,
         likelihood=likelihood,
         internal=rslt)
    , class="JD3SplitSA"))
}

q<-splitsa(sugar, 156, common.trend = T, common.seasonal = F, common.irregular = F)
par(mfrow=c(2,1))
plot(q$decomposition$series, col="gray")
lines(q$decomposition$level, col="red")
lines(q$decomposition$sa, col="blue")
plot(q$decomposition$seas, col="magenta")
lines(q$decomposition$noise, col="green")
par(mfrow=c(1,1))

q<-splitsa(log(sugar), 156, common.trend = T, common.seasonal = F, common.irregular = F)
par(mfrow=c(2,1))
plot(q$decomposition$series, col="gray")
lines(q$decomposition$level, col="red")
lines(q$decomposition$sa, col="blue")
plot(q$decomposition$seas, col="magenta")
lines(q$decomposition$noise, col="green")
par(mfrow=c(1,1))
