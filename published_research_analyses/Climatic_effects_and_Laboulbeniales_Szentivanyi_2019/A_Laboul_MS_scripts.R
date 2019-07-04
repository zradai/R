
lb<-read.csv("/media/zozoo/3191D18644760ED0/Tan/PhD/Others/Laboulbeniales_clim/A_newest_results_MERRAclim/labinf_data_V8_wFungiNames_merraclim.csv", sep=";")
pr<-read.csv("/media/zozoo/3191D18644760ED0/Tan/PhD/Others/Laboulbeniales_clim/A_newest_results_MERRAclim/labinf_PREVALENCE_V3_merraclim.csv", sep=";")

require("corrplot")
require("lme4")
require("lmerTest")

################################################################################################

### Presence / absence ###

###----- Check which variables seem to have significant effect on Laboulbeniales presence -----

# Mann-Whitney's tests
asd<-data.frame(variable=0, W=0, p=0)
for(v in 1:19){
  kwt.v<-wilcox.test(lb[lb$Labinf==0,(v+32)], lb[lb$Labinf==1,(v+32)])
  asd[v,"variable"]<-colnames(lb)[(v+32)]
  asd[v,"W"]<-as.numeric(kwt.v$statistic)
  asd[v,"p"]<-as.numeric(kwt.v$p.value)
}
asd$p.adj<-p.adjust(asd$p, method="holm")

###----- Principal component analyses (PCA) -----

# build matrix from relevant variables 
pre.prds<-c()
s.vars<-asd$variable[asd$p.adj<=0.05]
for(v in s.vars){
  pre.prds<-c(pre.prds, lb[,v])
}
preds<-matrix(pre.prds, ncol=nrow(asd[asd$p.adj<=0.05,]))

# PCA (on mean+SD-scaled variables to bring them to the same scale)
pc.preds<-prcomp(preds, scale=TRUE, center=TRUE)

# save principel components with Eigen-values higher than 1
lb$pc.preds1<-pc.preds$x[,1]
lb$pc.preds2<-pc.preds$x[,2]

###----- Visualise correlations between raw variables and principal components -----

var.names = asd$variable[asd$p.adj<=0.05]
vars.matrix = cbind(lb[,var.names], pc.preds$x[,1:2])

cor.a<-cor(vars.matrix)
colnames(cor.a)<-c(substr(var.names, 1, nchar(var.names)-2), "PC1-A", "PC2-A")
corrplot.mixed(cor.a, tl.col="black", upper.col=colorRampPalette(c("blue","white", "red"))(100), lower.col=colorRampPalette(c("blue","white", "red"))(100))

###----- Model on the effect of PCs on the presence of Laboulbeniales infection -----

pam1<-glmer(Labinf~pc.preds1+pc.preds2+(1|HOST), family="binomial", data=lb, 
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e9,npt=5)), na.action="na.fail")
summary(pam1)


################################################################################################

### Prevalence ###

###----- Principal component analyses (PCA) -----

# build matrix from relevant variables 
pre.prds<-c()
for(v in s.vars){
  pre.prds<-c(pre.prds, pr[,v])
}
preds<-matrix(pre.prds, ncol=nrow(asd[asd$p.adj<=0.05,]))

# PCA (on mean+SD-scaled variables to bring them to the same scale)
pc.preds.2<-prcomp(preds, scale=TRUE, center=TRUE)

# save principel components with Eigen-values higher than 1
pr$pc.prev1<-pc.preds.2$x[,1]
pr$pc.prev2<-pc.preds.2$x[,2]

###----- Model on the effect of PCs on the presence of Laboulbeniales infection -----

prem1<-glmer(prevalence~pc.prev1+pc.prev2+(1|host), family="binomial", data=pr, 
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e9,npt=5)), na.action="na.fail")
summary(prem1)
