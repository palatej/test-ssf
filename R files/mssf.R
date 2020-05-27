library(rjdssf)

### Multivariate model 
a<-read.table("./Data/bematrix.txt", sep='\t')
#a<-scale(a)

### normalize a. The results are much more stable

mdata<-cbind(a[,1], a[,10], a[,3], a[,4])

# create the model
model<-rjdssf::model()

# create the components and add them to the model
rjdssf::add(model, rjdssf::locallineartrend("tu", levelVariance = 0, fixedLevelVariance = T ))
rjdssf::add(model, rjdssf::locallineartrend("ty", levelVariance = 0, fixedLevelVariance = T))
rjdssf::add(model, rjdssf::locallevel("tpicore"))
rjdssf::add(model, rjdssf::locallevel("tpi"))
rjdssf::add(model, rjdssf::ar("cycle", c(1, -.5), fixedar = FALSE, variance= 1, fixedvariance=TRUE, nlags= 5))

# create the equations 
eq1<-rjdssf::equation("eq1", 1, F)
rjdssf::add(eq1, "tu")
rjdssf::add(eq1, "cycle", .1, F)
rjdssf::add(model, eq1)
eq2<-rjdssf::equation("eq2", 1, F)
rjdssf::add(eq2, "ty")
rjdssf::add(eq2, "cycle", .1, F)
rjdssf::add(model, eq2)
eq3<-rjdssf::equation("eq3", 1, F)
rjdssf::add(eq3, "tpicore")
rjdssf::add(eq3, "cycle", .1, F, rjdssf::loading(4))
rjdssf::add(model, eq3)
eq4<-rjdssf::equation("eq4", 1, F)
rjdssf::add(eq4, "tpi")
rjdssf::add(eq4, "cycle", .1, F)
rjdssf::add(model, eq4)

#estimate the model
rslt2<-rjdssf::estimate(model, mdata, marginal=T, concentrated=T)

print(result(rslt2, "parameters"))
print(result(rslt2, "likelihood.ll"))
factor=sqrt(result(rslt2, "scalingfactor"))
print(factor)
ts.plot(-result(rslt2, "ssf.smoothing.cmp(4)"))
