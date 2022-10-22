
###----- FUNCTIONS -----

rsc = function(x, na.rm=T){
  x.rsc = (x-mean(x, na.rm=na.rm))/sd(x, na.rm=na.rm)
  return(x.rsc)
}

tbeta = function(x){
  x.t = (x*(length(x)-1)+0.5)/length(x)
  return(x.t)
}

range2 = function(x){
  range2.x = c(
    min(x, na.rm = T),
    mean(x, na.rm = T),
    max(x, na.rm = T)
  )
  return(range2.x)
}

#-------------------------------------------------------------------------------

###----- DATA -----

evar.labels = data.frame(
  label = c("Error rate", "Sequencing depth", "PCR duplicate ratio", "Optical duplicate ratio"),
  evar = c("err", "cov", "pdup", "odup")
)

gvar.labels = data.frame(
  label = c("Genome size (Mbp)", "GC%", "Unique ratio"),
  gvar = c("size.mb.full", "GC", "complexity")
)

bacs = read.csv("bacteria.csv", sep = ";")
asq = read.csv("genome_error_simulations.csv", sep = ";")
for(v in evar.labels$evar){
  asq[,paste(v,"rsc", sep=".")] = rsc(asq[,v])
}

qvar.labels = read.csv("quality_metrics_labels.csv", sep = ";")
qvar.labels = qvar.labels[!(qvar.labels$qvar %in% vars.to.omit),]
qvars = qvar.labels$qvar

N.nominal = length(unique(asq$cov)) * length(unique(asq$err)) * length(unique(asq$odup)) * length(unique(asq$pdup)) * length(unique(asq$rep))

vals = data.frame(
  evar = "",
  evar.label = "",
  original.value = 0,
  rescaled.value = 0,
  original.mean = 0,
  original.sd = 0
)[-1,]
for(e in evar.labels$evar){
  for(v in unique(asq[,e])){
    vals[nrow(vals)+1,"evar"] = e
    vals[nrow(vals),"evar.label"] = evar.labels$label[evar.labels$evar==e]
    vals[nrow(vals),"original.value"] = v
    vals[nrow(vals),"rescaled.value"] = unique(asq[asq[,e]==v,paste0(e,".rsc")])
    vals[nrow(vals),"original.mean"] = mean(asq[,e])
    vals[nrow(vals),"original.sd"] = sd(asq[,e])
  }
}

#--------------------------------------------------------------------------------

###----- Additive models -----

models.add = list()
for(b in unique(bacs$working.ID)){
  print(b)
  for(y in qvar.labels$qvar){
    
    df.by = asq[asq$bact.ID==b,]
    df.by$y = df.by[,y]
    df.by = df.by[!is.na(df.by$y),]
    df.by$y = tbeta(df.by$y/max(df.by$y))
    
    mod.by = betareg(y ~ cov.rsc + err.rsc + odup.rsc + pdup.rsc, data = df.by,
                     control = betareg.control(maxit = 1e4, fsmaxit = 2e3, fstol = 1e-8))
    
    models.add[[paste0(b,"_",y)]] = mod.by
    
  }
}

###----- Multiplicative models -----

models.mult = list()
for(b in unique(bacs$working.ID)){
  print(b)
  for(y in qvar.labels$qvar){
    
    df.by = asq[asq$bact.ID==b,]
    df.by$y = df.by[,y]
    df.by = df.by[!is.na(df.by$y),]
    df.by$y = tbeta(df.by$y/max(df.by$y))
    
    mod.by = betareg(y ~ cov.rsc * err.rsc * odup.rsc * pdup.rsc, data = df.by,
                     control = betareg.control(maxit = 1e4, fsmaxit = 2e3, fstol = 1e-8))
    
    models.mult[[paste0(b,"_",y)]] = mod.by
    
  }
}


###----- create EMT tables -----

# additive models
emt.add = data.frame(bact.ID = "", qvar = "", evar = "", estimate = 0, SE = 0, CI.low = 0, CI.up = 0, z = 0, P = 0, N = 0)[-1,]
for(i in 1:length(models.add)){
  
  cat(i," / ",length(models.add)," (",round(i/length(models.add)*100),"%)\n", sep="")
  
  nm.i = base::strsplit(x = names(models.add)[i], split = "_")[[1]]
  bac.i = nm.i[1]
  y.i = nm.i[2]
  
  coefs.i = as.data.frame(summary(models.add[[i]])$coefficients)
  ci.i = stats::confint(models.add[[i]], level = 0.997)
  for(j in paste0(evar.labels$evar, ".rsc")){
    emt.add[nrow(emt.add)+1,"bact.ID"] = bac.i
    emt.add[nrow(emt.add),"qvar"] = y.i
    emt.add[nrow(emt.add),"evar"] = j
    emt.add[nrow(emt.add),"estimate"] = coefs.i[j,"mean.Estimate"]
    emt.add[nrow(emt.add),"SE"] = coefs.i[j,"mean.Std..Error"]
    emt.add[nrow(emt.add),"CI.low"] = ci.i[j,1]
    emt.add[nrow(emt.add),"CI.up"] = ci.i[j,2]
    emt.add[nrow(emt.add),"z"] = coefs.i[j,"mean.z.value"]
    emt.add[nrow(emt.add),"P"] = coefs.i[j,"mean.Pr...z.."]
    emt.add[nrow(emt.add),"N"] = nrow(asq[!is.na(asq[,y.i]),])
  }
  
}
for(b in unique(emt.add$bact.ID)){
  for(g in gvar.labels$gvar){
    emt.add[emt.add$bact.ID==b,g] = bacs[bacs$working.ID==b,g]
  }
}

# multiplicative models
emt.mult = data.frame(bact.ID = "", qvar = "", focal.evar = "", estimate = 0, SE = 0, CI.low = 0, CI.up = 0, z = 0, P = 0, N = 0)
focal.evar = "err.rsc"
marginal.vars = paste0(evar.labels$evar, ".rsc")
marginal.vars = marginal.vars[!grepl(pattern = focal.evar, x = marginal.vars)]
evar.levels = list()
for(v in marginal.vars){
  emt.mult[1,v] = 0
  evar.levels[[v]] = range2(asq[,v])
}
emt.mult = emt.mult[-1,]
coefs.mult = data.frame(bact.ID = "", qvar = "", coef = "", estimate = 0, SE = 0, CI.low = 0, CI.up = 0, z = 0, P = 0, N = 0)[-1,]
for(i in 1:length(models.mult)){
  
  cat(i," / ",length(models.mult)," (",round(i/length(models.mult)*100),"%)\n", sep="")
  
  nm.i = base::strsplit(x = names(models.mult)[i], split = "_")[[1]]
  bac.i = nm.i[1]
  y.i = nm.i[2]
  
  df.by = asq[asq$bact.ID==bac.i,]
  df.by$y = df.by[,y.i]
  df.by = df.by[!is.na(df.by$y),]
  df.by$y = tbeta(df.by$y/max(df.by$y))
  
  # coefs
  coefs.i = as.data.frame(summary(models.mult[[i]])$coefficients)
  ci.i = stats::confint(models.mult[[i]], level = 0.997)
  for(j in rownames(coefs.i)[-1]){
    coefs.mult[nrow(coefs.mult)+1,"bact.ID"] = bac.i
    coefs.mult[nrow(coefs.mult),"qvar"] = y.i
    coefs.mult[nrow(coefs.mult),"coef"] = j
    coefs.mult[nrow(coefs.mult),"estimate"] = coefs.i[j,"mean.Estimate"]
    coefs.mult[nrow(coefs.mult),"SE"] = coefs.i[j,"mean.Std..Error"]
    coefs.mult[nrow(coefs.mult),"CI.low"] = ci.i[j,1]
    coefs.mult[nrow(coefs.mult),"CI.up"] = ci.i[j,2]
    coefs.mult[nrow(coefs.mult),"z"] = coefs.i[j,"mean.z.value"]
    coefs.mult[nrow(coefs.mult),"P"] = coefs.i[j,"mean.Pr...z.."]
    coefs.mult[nrow(coefs.mult),"N"] = nrow(df.by)
  }
  
  # EMTs
  emt.i = as.data.frame(emtrends(models.mult[[i]], pairwise ~ cov.rsc+odup.rsc+pdup.rsc, var=focal.evar, at=evar.levels, mode="link", level=0.997, data = df.by)$emtrends)
  emt.i.p = as.data.frame(emmeans::test(emtrends(models.mult[[i]], pairwise ~ cov.rsc+odup.rsc+pdup.rsc, var=focal.evar, at=evar.levels, mode="link", level=0.997, data = df.by)$emtrends))
  for(j in unique(emt.i[,marginal.vars[1]])){
    for(k in unique(emt.i[,marginal.vars[2]])){
      for(l in unique(emt.i[,marginal.vars[3]])){
        
        logic.index = emt.i[,marginal.vars[1]]==j & emt.i[,marginal.vars[2]]==k & emt.i[,marginal.vars[3]]==l
        
        emt.mult[nrow(emt.mult)+1,"bact.ID"] = bac.i
        emt.mult[nrow(emt.mult),"qvar"] = y.i
        emt.mult[nrow(emt.mult),"focal.evar"] = focal.evar
        emt.mult[nrow(emt.mult),"estimate"] = emt.i[logic.index, grepl(pattern = focal.evar, x = colnames(emt.i))]
        emt.mult[nrow(emt.mult),"SE"] = emt.i[logic.index, "SE"]
        emt.mult[nrow(emt.mult),"CI.low"] = emt.i[logic.index, "asymp.LCL"]
        emt.mult[nrow(emt.mult),"CI.up"] = emt.i[logic.index, "asymp.UCL"]
        emt.mult[nrow(emt.mult),"z"] = emt.i.p[logic.index, "z.ratio"]
        emt.mult[nrow(emt.mult),"P"] = emt.i.p[logic.index, "p.value"]
        emt.mult[nrow(emt.mult),"N"] = nrow(df.by)
        
        marg.values = c(j,k,l)
        names(marg.values) = c(marginal.vars)
        for(w in marginal.vars){
          w0 = gsub(pattern = ".rsc", replacement = "", x = w)
          original.mu = unique(vals$original.mean[vals$evar==w0])
          original.sd = unique(vals$original.sd[vals$evar==w0])
          original.val = marg.values[w]*original.sd + original.mu
          emt.mult[nrow(emt.mult),w0] = original.val
          emt.mult[nrow(emt.mult),w] = emt.i[logic.index, w]
        }
        
      }
    }
  }
}
for(b in unique(emt.add$bact.ID)){
  for(g in gvar.labels$gvar){
    emt.mult[emt.mult$bact.ID==b,g] = bacs[bacs$working.ID==b,g]
  }
}


###----- meta-regressions -----

MA.method = "REML"

#--- additive ---

# no moderators
meta.add = data.frame(qvar = "", evar = "", estimate = 0, SE = 0, CI.low = 0, CI.up = 0, z = 0, P = 0, tau2 = 0, I2 = 0, H2 = 0, N = 0)[-1,]
for(q in qvar.labels$qvar){
  
  print(q)
  
  for(e in evar.labels$evar){
    
    emt.qe = emt.add[emt.add$qvar==q & emt.add$evar==paste(e,"rsc", sep="."),]
    mr.qe = rma(yi = emt.qe$estimate, sei = emt.qe$SE, slab = emt.qe$bact.ID, method = MA.method)
    mr.ci = confint(mr.qe, level = 0.997, fixed=T, random=F)
    
    meta.add[nrow(meta.add)+1,"qvar"] = q
    meta.add[nrow(meta.add),"evar"] = paste(e,"rsc", sep=".")
    meta.add[nrow(meta.add),"estimate"] = mr.qe$beta
    meta.add[nrow(meta.add),"SE"] = mr.qe$se
    meta.add[nrow(meta.add),"CI.low"] = mr.ci$fixed[1,"ci.lb"]
    meta.add[nrow(meta.add),"CI.up"] = mr.ci$fixed[1,"ci.ub"]
    meta.add[nrow(meta.add),"z"] = mr.qe$zval
    meta.add[nrow(meta.add),"P"] = mr.qe$pval
    meta.add[nrow(meta.add),"tau2"] = mr.qe$tau2
    meta.add[nrow(meta.add),"I2"] = mr.qe$I2
    meta.add[nrow(meta.add),"H2"] = mr.qe$H2
    meta.add[nrow(meta.add),"N"] = nrow(emt.qe)
    
  }
}

# WITH moderators
meta.add.WithMod = data.frame(qvar = "", evar = "", gvar = "", mod.var = "", estimate = 0, SE = 0, CI.low = 0, CI.up = 0, z = 0, P = 0, tau2 = 0, I2 = 0, H2 = 0, N = 0)[-1,]
for(q in qvar.labels$qvar){
  
  print(q)
  
  for(e in evar.labels$evar){
    
    emt.qe = emt.add[emt.add$qvar==q & emt.add$evar==paste(e,"rsc", sep="."),]
    
    for(g in gvar.labels$gvar){
      
      moderator.form = as.formula(paste0("~",gsub(pattern = "GC", replacement = "poly(GC,2)", x = g)))
      mr.qe = rma(yi = emt.qe$estimate, sei = emt.qe$SE, slab = emt.qe$bact.ID, method = MA.method, mods = moderator.form, data = emt.qe)
      mr.ci = confint(mr.qe, level = 0.997, fixed = T, random = F)
      
      for(m in rownames(mr.qe$beta)){
        
        meta.add.WithMod[nrow(meta.add.WithMod)+1,"qvar"] = q
        meta.add.WithMod[nrow(meta.add.WithMod),"evar"] = e
        meta.add.WithMod[nrow(meta.add.WithMod),"gvar"] = g
        meta.add.WithMod[nrow(meta.add.WithMod),"mod.var"] = m
        meta.add.WithMod[nrow(meta.add.WithMod),"estimate"] = mr.qe$beta[m,1]
        meta.add.WithMod[nrow(meta.add.WithMod),"SE"] = mr.qe$se[grep(pattern = m, x = rownames(mr.qe$beta), fixed = T)]
        meta.add.WithMod[nrow(meta.add.WithMod),"CI.low"] = mr.ci$fixed[m,2]
        meta.add.WithMod[nrow(meta.add.WithMod),"CI.up"] = mr.ci$fixed[m,3]
        meta.add.WithMod[nrow(meta.add.WithMod),"z"] = mr.qe$zval[grep(pattern = m, x = rownames(mr.qe$beta), fixed = T)]
        meta.add.WithMod[nrow(meta.add.WithMod),"P"] = mr.qe$pval[grep(pattern = m, x = rownames(mr.qe$beta), fixed = T)]
        meta.add.WithMod[nrow(meta.add.WithMod),"tau2"] = mr.qe$tau2
        meta.add.WithMod[nrow(meta.add.WithMod),"I2"] = mr.qe$I2
        meta.add.WithMod[nrow(meta.add.WithMod),"H2"] = mr.qe$H2
        meta.add.WithMod[nrow(meta.add.WithMod),"N"] = nrow(emt.qe)
        
      }
    }
  }
}


#--- multiplicative ---

# no moderators
meta.mult = data.frame(qvar = "", focal.evar = "", estimate = 0, SE = 0, CI.low = 0, CI.up = 0, z = 0, P = 0, tau2 = 0, I2 = 0, H2 = 0, N = 0)
for(v in marginal.vars){
  meta.mult[1,v] = 0
}
meta.mult = meta.mult[-1,]
for(q in qvar.labels$qvar){
  
  print(q)
  
  for(i in unique(emt.mult[,marginal.vars[1]])){
    for(j in unique(emt.mult[,marginal.vars[2]])){
      for(k in unique(emt.mult[,marginal.vars[3]])){
        
        emt.qe = emt.mult[emt.mult$qvar==q & emt.mult[,marginal.vars[1]]==i & emt.mult[,marginal.vars[2]]==j & emt.mult[,marginal.vars[3]]==k,]
        
        mr.qe = rma(yi = emt.qe$estimate, sei = emt.qe$SE, slab = emt.qe$bact.ID, method = MA.method)
        mr.ci = confint(mr.qe, level = 0.997, fixed=T, random=F)
        
        meta.mult[nrow(meta.mult)+1,"qvar"] = q
        meta.mult[nrow(meta.mult),"focal.evar"] = focal.evar
        meta.mult[nrow(meta.mult),"estimate"] = mr.qe$beta
        meta.mult[nrow(meta.mult),"SE"] = mr.qe$se
        meta.mult[nrow(meta.mult),"CI.low"] = mr.ci$fixed[1,"ci.lb"]
        meta.mult[nrow(meta.mult),"CI.up"] = mr.ci$fixed[1,"ci.ub"]
        meta.mult[nrow(meta.mult),"z"] = mr.qe$zval
        meta.mult[nrow(meta.mult),"P"] = mr.qe$pval
        meta.mult[nrow(meta.mult),"tau2"] = mr.qe$tau2
        meta.mult[nrow(meta.mult),"I2"] = mr.qe$I2
        meta.mult[nrow(meta.mult),"H2"] = mr.qe$H2
        meta.mult[nrow(meta.mult),"N"] = nrow(emt.qe)
        
        marg.values = c(i,j,k)
        names(marg.values) = c(marginal.vars)
        for(w in marginal.vars){
          w0 = gsub(pattern = ".rsc", replacement = "", x = w)
          original.mu = unique(vals$original.mean[vals$evar==w0])
          original.sd = unique(vals$original.sd[vals$evar==w0])
          original.val = marg.values[w]*original.sd + original.mu
          meta.mult[nrow(meta.mult),w0] = original.val
          meta.mult[nrow(meta.mult),w] = marg.values[w]
        }
        
      }
    }
  }
}


# WITH moderators
meta.mult.WithMod = data.frame(qvar = "", focal.evar = "", gvar = "", mod.var = "", estimate = 0, SE = 0, CI.low = 0, CI.up = 0, z = 0, P = 0, tau2 = 0, I2 = 0, H2 = 0, N = 0)
for(v in marginal.vars){
  meta.mult.WithMod[1,v] = 0
}
meta.mult.WithMod = meta.mult.WithMod[-1,]
for(q in qvar.labels$qvar){
  
  print(q)
  
  for(i in unique(emt.mult[,marginal.vars[1]])){
    for(j in unique(emt.mult[,marginal.vars[2]])){
      for(k in unique(emt.mult[,marginal.vars[3]])){
        
        emt.qe = emt.mult[emt.mult$qvar==q & emt.mult[,marginal.vars[1]]==i & emt.mult[,marginal.vars[2]]==j & emt.mult[,marginal.vars[3]]==k,]
        
        for(g in gvar.labels$gvar){
          
          moderator.form = as.formula(paste0("~",gsub(pattern = "GC", replacement = "poly(GC,2)", x = g)))
          mr.qe = rma(yi = emt.qe$estimate, sei = emt.qe$SE, slab = emt.qe$bact.ID, method = MA.method, mods = moderator.form, data = emt.qe)
          mr.ci = confint(mr.qe, level = 0.997, fixed = T, random = F)
          
          for(m in rownames(mr.qe$beta)){
            meta.mult.WithMod[nrow(meta.mult.WithMod)+1,"qvar"] = q
            meta.mult.WithMod[nrow(meta.mult.WithMod),"focal.evar"] = focal.evar
            meta.mult.WithMod[nrow(meta.mult.WithMod),"gvar"] = g
            meta.mult.WithMod[nrow(meta.mult.WithMod),"mod.var"] = m
            meta.mult.WithMod[nrow(meta.mult.WithMod),"estimate"] = mr.qe$beta[m,1]
            meta.mult.WithMod[nrow(meta.mult.WithMod),"SE"] = mr.qe$se[grep(pattern = m, x = rownames(mr.qe$beta), fixed = T)]
            meta.mult.WithMod[nrow(meta.mult.WithMod),"CI.low"] = mr.ci$fixed[m,2]
            meta.mult.WithMod[nrow(meta.mult.WithMod),"CI.up"] = mr.ci$fixed[m,3]
            meta.mult.WithMod[nrow(meta.mult.WithMod),"z"] = mr.qe$zval[grep(pattern = m, x = rownames(mr.qe$beta), fixed = T)]
            meta.mult.WithMod[nrow(meta.mult.WithMod),"P"] = mr.qe$pval[grep(pattern = m, x = rownames(mr.qe$beta), fixed = T)]
            meta.mult.WithMod[nrow(meta.mult.WithMod),"tau2"] = mr.qe$tau2
            meta.mult.WithMod[nrow(meta.mult.WithMod),"I2"] = mr.qe$I2
            meta.mult.WithMod[nrow(meta.mult.WithMod),"H2"] = mr.qe$H2
            meta.mult.WithMod[nrow(meta.mult.WithMod),"N"] = nrow(emt.qe)
            
            marg.values = c(i,j,k)
            names(marg.values) = c(marginal.vars)
            for(w in marginal.vars){
              w0 = gsub(pattern = ".rsc", replacement = "", x = w)
              original.mu = unique(vals$original.mean[vals$evar==w0])
              original.sd = unique(vals$original.sd[vals$evar==w0])
              original.val = marg.values[w]*original.sd + original.mu
              meta.mult.WithMod[nrow(meta.mult.WithMod),w0] = original.val
              meta.mult.WithMod[nrow(meta.mult.WithMod),w] = marg.values[w]
            }
            
          }
          
        }
        
      }
    }
  }
}

#------------------------------------------------------------------------------------------------

save(rsc, tbeta, range2, evar.labels, bacs, asq, qvar.labels, gvar.labels, N.nominal, vals,
     models.add, models.mult, emt.add, emt.mult, coefs.mult,
     meta.add, meta.add.WithMod, meta.mult, meta.mult.WithMod,
     file = "load_vars.RData")
