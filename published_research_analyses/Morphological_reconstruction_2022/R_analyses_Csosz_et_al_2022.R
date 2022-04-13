
###----- loading data -----

# load data frame with morphological (trait) variables
tt = read.csv("Csosz_et_al_Temnothorax_morpho_data.csv", sep=";")

# variables to work with
vars = c("CL","POC","CW","FR","SL","ML","PEH","PEL","NOH","NOL","PPH","PPL","SPST","SPL","MW","PEW","PPW","SPBA","SPWI","SPTI")
spp = unique(tt$species)

rsc = function(x){
  x.rsc = (x-mean(x, na.rm=T))/sd(x, na.rm=T)
  return(x.rsc)
}

###----- PCA -----

pca.1 = stats::prcomp(x = tt[,vars], center = TRUE, scale = TRUE)
summary(pca.1)

tt = cbind(tt, pca.1$x[,1:2])

###----- preliminary analytical step -----

require("tidyr")

tt.long = gather(tt, trait, measurement, CL:SPTI)

m1 = lm(measurement ~ as.factor(infected)*species*trait, data = tt.long)

stats::anova(m1)


###----- fit models and acquire marginal mean estimates -----

df.t = tt
model.list = list()
anova.list = list()
emm.list = list()
for(t in vars){
  
  df.t$y = df.t[,t]
  for(s in spp){
    df.t$y[df.t$species==s] = rsc(df.t$y[df.t$species==s])
  }
  
  # save model fit
  model.list[[t]] = lm(y ~ as.factor(infected)*species, data = df.t)
  
  par(mfrow=c(2,1))
  plot(fitted(model.list[[t]]), resid(model.list[[t]]), xlab="Fitted values", ylab="Residuals")
  hist(resid(model.list[[t]]), xlab="Residuals", main = paste0("Residuals for ", t))
  par(mfrow=c(1,1))
  
  # save anova table for the fitted model
  anova.list[[t]] = stats::anova(model.list[[t]])
  
  # save estimated marginal means and contrasts from the fitted model
  emm.list[[t]] = emmeans(model.list[[t]], pairwise ~ as.factor(infected)|species)
  
  # ask for input before stepping to next trait
  readline(paste0("Current trait: ",t,". Press ENTER for next trait!"))
  
}

# apply Bonferroni correction on P-values of interaction terms from ANOVA tables
anova.Pvals = unlist(lapply(anova.list, function(x){return(x[3,5])}))
stats::p.adjust(anova.Pvals, method = "bonferroni")

# apply Bonferroni correction on P-values marginal contrast terms from EMMs
emm.Pvals = data.frame(trait="", species="", P=0)[-1,]
for(t in vars){
  for(l in 1:3){
    contr.tl = emm.list[[t]]$contrasts[l]
    emm.Pvals[nrow(emm.Pvals)+1,"trait"] = t
    emm.Pvals[nrow(emm.Pvals),"species"] = as.character(as.data.frame(contr.tl)$species)
    emm.Pvals[nrow(emm.Pvals),"P"] = as.data.frame(contr.tl)$p.value
  }
}
stats::p.adjust(emm.Pvals$P, method = "bonferroni")


###----- estimate putative uninfected trait values -----

loglik = function(par, y) {
  sum(dnorm(y, par[1], sqrt(par[2]), log = TRUE))
}

# new table in which "restored" trait values will be saved
tt.2 = tt[tt$infected==1,1:24]
for(s in spp){
  
  tt.s = tt[tt$species==s,]
  
  for(t in vars){
    
    # estimate mean and variance of value distributions, for given species and trait, separately for infested and non-infested
    pars.wild = optim(c(mean = 0, var = 1), fn = loglik, y = tt.s[tt.s$infected==0,t], control = list(fnscale = -1, reltol = 1e-16))$par
    pars.inf = optim(c(mean = 0, var = 1), fn = loglik, y = tt.s[tt.s$infected==1,t], control = list(fnscale = -1, reltol = 1e-16))$par
    
    # calculate the SD-standardized difference between infested and non-infested
    d.inf = (tt.s[tt.s$infected==1,t] - pars.inf[1])/sqrt(pars.inf[2])
    S = lm(tt.s[,t] ~ as.factor(tt.s$infected))$coefficients[2] * (-1)
    
    # calculate predicted non-infested trait values
    trait.inf.shifted = d.inf*sqrt(pars.wild[2]) + pars.inf[1] + S
    
    # save new values
    tt.2[tt.2$species==s,t] = trait.inf.shifted
    
  }
  
}
