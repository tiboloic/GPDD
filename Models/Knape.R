# LT 26/08/2019

# fits full model, no density dependence model and randomb model to Knape's dataset 

load("TaxonomyData.Rdata")

library("TMB")
#compile("commonb.cpp")
dyn.load(dynlib("commonb"))

# remove family
models = c("full")
modelfactors = list(full = factor(1:627)) 

allmodelfits=c()

# fit full model
mus = sapply(seq_len(nrow(obs)), function(i) mean(obs[i,!nas[i,]]))

for (model in models) {
  
  modelfactor = modelfactors[[model]]
  
  # to store the fits
  onemodelfits = c()
  
  for (group in unique(modelfactor)) {
    
    cat(paste("Now fitting model ", model," group ", group, "\n"))
    
    # select a group
    sel = modelfactor==group
    
    # # of time series in that group
    n=sum(sel)
    
    # get intit values from regularized fits
    init = list(b=  .8, mu=mus[sel], lsig=  rep(0.1,n), ltau=  rep(0.1,n))
    #fit group
    f = MakeADFun(
      data=list(obs= obs[sel,,drop=F], lengths= lengths[sel], nas= nas[sel,,drop=F]),
      parameters=init,
      DLL="commonbRepar", silent=T)
    
#    fit = optim(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))
    fit = nlminb(f$par,f$fn,f$gr, control=list(eval.max=5000, iter.max=3000))
    
    # get likelihood profile for b
    prof = tmbprofile(f,"b", trace=F)
    
    # get standard erroes
    sdrep = sdreport(f)
    
    # save fit, likelihod profile
    onemodelfits[[group]] = list(name=group, n=n, fit=fit, f=f, prof=prof, sd=sdrep)
    
  }
  allmodelfits[[model]] = onemodelfits
}

# now do the density independent stuff
n = 627
sel = rep(TRUE,n)

init = allmodelfits.regul$nodd$nodd$fit$par
init = list(b=1,
            a = init[names(init)=="a"],
            lsig = init[names(init)=="lsig"],
            ltau = init[names(init)=="ltau"])

# allows for geometric growth/decline
f = MakeADFun(
  data=list(obs= obs, lengths= lengths, nas= nas),
  parameters=init,
  map = list(b=factor(NA)),
  DLL="commonb", silent=T)
fit = optim(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))

# get standard errors
sdrep = sdreport(f)

onemodelfits=c()
onemodelfits[["nodd"]] = list(name="nodd", n=n, fit=fit, f=f, prof=NA, sd=sdrep)
allmodelfits[["nodd"]] = onemodelfits

#########
# random effect model fit
########
#compile("RandombRepar.cpp")
dyn.load(dynlib("RandombRepar"))

n = 627
init = allmodelfits.regul$randomb$randomb$f$env$parList()

f = MakeADFun(
  data=list(obs= obs, lengths= lengths, nas= nas),
  parameters=init,
  random=c("bdevs"),
  DLL="RandombRepar", silent=T)
fit = optim(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))
sdrep = sdreport(f)
allmodelfits[["randomb"]]=list("randomb" = list(name="randomb", n=n, fit=fit, f=f, prof=NA, sd=sdrep))

####
# model selection
nLL = sapply(allmodelfits, function(model) sum(sapply(model, function(fit) fit$fit$objective)))
df = sapply(allmodelfits, function(model) sum(sapply(model, function(fit) length(fit$fit$par))))
nobs = sum(!nas)
AIC = 2 * df + 2 * nLL
AICc = AIC + 2*df*(df+1)/(nobs-df-1)
BIC = log(nobs) *df + 2 * nLL

