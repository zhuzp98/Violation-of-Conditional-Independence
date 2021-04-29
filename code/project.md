# Settings

    require(rjags)
    require(MCMCvis)
    ### global reproducibility
    set.seed(1234)

    logit = function(p) { log(p)-log(1-p) }
    expit = function(z) { 1/(1+exp(-z)) }

# Data

    ### data size
    n = 5000

    ### balanced case-control study
    y = sample(c(rep(0, n/2), rep(1, n/2)))

    ### true OR of 2
    gamma.0.TR = 0.1
    gamma.1.TR = expit(logit(gamma.0.TR)+log(2))

    x = rbinom(n, size=1, prob=(1-y)*gamma.0.TR + y*gamma.1.TR)

    ### two independent surrogates M1
    xstr01 = rbinom(n, size=1, prob=(1-(x))*(1-0.85)+(x)*(0.75))
    xstr02 = rbinom(n, size=1, prob=(1-(x))*(1-0.95)+(x)*(0.60))

    ### two semi_independent surrogates M2
    xstr11 = rbinom(n, size=1, prob=(1-(0.9*x+0.1*y))*(1-0.85)+(0.9*x+0.1*y)*(0.75))
    xstr12 = rbinom(n, size=1, prob=(1-(0.9*x+0.1*y))*(1-0.95)+(0.9*x+0.1*y)*(0.60))

    ### two semi_independent surrogates M3
    xstr21 = rbinom(n, size=1, prob=(1-(0.4*x+0.6*y))*(1-0.85)+(0.4*x+0.6*y)*(0.75))
    xstr22 = rbinom(n, size=1, prob=(1-(0.4*x+0.6*y))*(1-0.95)+(0.4*x+0.6*y)*(0.60))

    ### Mother Nature WITHHOLDS the actual x
    rm(x)

    dta0 = data.frame(xstr1=xstr01, xstr2=xstr02, y=y) 
    dta1 = data.frame(xstr1=xstr11, xstr2=xstr12, y=y)
    dta2 = data.frame(xstr1=xstr21, xstr2=xstr22, y=y)

# Generative Model

    genmod.string <- "model{
      
    ### prior distribution
    gamma.0 ~ dunif(0,1)
    gamma.1 ~ dunif(0,1)
    sn1 ~ dunif(0.5, 1)
    sp1 ~ dunif(0.5, 1)
    sn2 ~ dunif(0.5, 1)
    sp2 ~ dunif(0.5, 1)

    trgt <- logit(gamma.1)-logit(gamma.0)

    for (i in 1:n) {
      x[i] ~ dbern((1-y[i])*gamma.0+y[i]*gamma.1)
      xstr1[i] ~ dbern((1-x[i])*(1-sp1)+x[i]*sn1)
      xstr2[i] ~ dbern((1-x[i])*(1-sp2)+x[i]*sn2)
    }  

    }" 

    GEN_opt = function(data, modstring) {
      dta = data
      ### generative model, data go in
      mod <- jags.model(textConnection(modstring),
        data=list(y=dta$y, xstr1=dta$xstr1, xstr2=dta$xsrt2,
                  n=dim(dta)[1]),
        n.chains=3)

      ###  MC output comes out
      opt.JAGS <- coda.samples(mod, n.iter=20000, thin=10, 
        variable.names=c("gamma.0","gamma.1","sn1","sp1",
                        "sn2","sp2","trgt"))

      return(opt.JAGS)
    }

    dt0 = GEN_opt(dta0, genmod.string)
    dt1 = GEN_opt(dta1, genmod.string)
    dt2 = GEN_opt(dta2, genmod.string)

# Model Summary

    MCMCsummary(dt0)

    MCMCsummary(dt1)

    MCMCsummary(dt2)

# Plots

    MCMCtrace(dt0, params="trgt", pdf=F)

    MCMCtrace(dt1, params="trgt", pdf=F)

    MCMCtrace(dt2, params="trgt", pdf=F)

    MCMCplot(dt0, dt1)

    MCMCplot(dt0, dt2)
