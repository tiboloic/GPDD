# LT 26/08/2019

# fits random effects model and Hierarchical model to restricted dataset

load("TaxonomyData.Rdata")


# load fits on Knape's dataset
load("Knape.RData")

#keep only relevant time-series
allmodelfits$full = allmodelfits$full[order.sel]
allmodelfits$nodd = allmodelfits$full[order.sel]

library("TMB")

#########
# random effect model fit
########
#compile("RandombRepar.cpp")
dyn.load(dynlib("RandombRepar"))

init = allmodelfits.regul$randomb$randomb$f$env$parList()
init$mu = init$mu[order.sel]
init$lsig = init$lsig[order.sel]
init$ltau = init$ltau[order.sel]
init$bdevs = init$bdevs[order.sel]

f = MakeADFun(
  data=list(obs= obs, lengths= lengths, nas= nas),
  parameters=init,
  random=c("bdevs"),
  DLL="RandombRepar", silent=T)
fit = optim(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))
sdrep = sdreport(f)
allmodelfits[["randomb"]]=list("randomb" = list(name="randomb", n=n, fit=fit, f=f, prof=NA, sd=sdrep))


### hierarchical model

#compile("RandomOrder.cpp")
dyn.load(dynlib("RandomOrder"))


init = list(b=0.5, mu=mus, lsig = rep(-1,n),
            ltau = rep(-1,n),
            lnu0=-3, lnu1= -3,
            bdevs0=rep(0,14), bdevs1=rep(0,n))
f = MakeADFun(
  data=list(obs= obs, lengths= lengths, nas= nas, orders=order.fact),
  parameters=init,
  random=c("bdevs0", "bdevs1"),
  DLL="RandomOrder", silent=F, inner.control = list(maxit=1000))

fit = optim(g$par,g$fn,g$gr, method="BFGS", control=list(maxit=3000), hessian=T)
sdrep = sdreport(f)

allmodelfits[["hierarchical"]]=list("hierarchical" = list(name="hierarchical", n=n, fit=fit, f=f, prof=NA, sd=sdrep))

### model selection
nLL = sapply(allmodelfits, function(model) sum(sapply(model, function(fit) fit$fit$value)))
df = sapply(allmodelfits, function(model) sum(sapply(model, function(fit) length(fit$fit$par))))
nobs = sum(!nas)
AIC = 2 * df + 2 * nLL
AICc = AIC + 2*df*(df+1)/(nobs-df-1)
BIC = log(nobs) *df + 2 * nLL

save(allmodelfits, file="Hierarchical.Rdata")
