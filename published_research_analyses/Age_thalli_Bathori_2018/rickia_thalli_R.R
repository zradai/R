
t<-read.csv("/media/zozoo/3191D18644760ED0/Tan/PhD/Others/Formi/age_and_thalli_csv.csv", sep=";")
r<-read.csv("/media/zozoo/3191D18644760ED0/Tan/PhD/Others/Formi/rickia_teratology_2016.csv", sep=";")

require("lme4")
require("multcomp")

#------------------------------------------------------------------------------------------------------

###----- Q1: how age affects thalli number? -----

q1 = glmer.nb(thalli~age+(age|site/colony_id), 
              nb.control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e9,npt=5)), data=t)
summary(q1)

# get factor-level comparison not shown in model summary
summary(glht(q1, linfct=mcp(age=c("age5 - age3 = 0")) ))

###----- Q2: how age affects proportion of teratological thalli? -----

q2 = glm(abnorm.prop ~ age, family="binomial", data=r)
summary(q2)
summary(glht(q2, linfct=mcp(age=c("age5 - age3 = 0")) ))

###----- Q3: how age affects proportion of melanised thalli? -----

q3 = glm(melan.prop ~ age, family="binomial", data=r)
summary(q3)
summary(glht(q3, linfct=mcp(age=c("age5 - age3 = 0")) ))


###----- Visualising age-dependent thalli distributions -----

t$Yp<-predict(q1, type="response")
par(mfrow=c(3,3))
for(s in unique(t$site)){
  for(c in unique(t$colony_id[t$site==s])){
    with(t[t$site==s & t$colony_id==c,], boxplot(thalli~age, cex=1.5, lwd=1.5, cex.lab=1.5, xlab="Age group", ylab="Thalli number", main=paste(s," - ",c,sep="")))
    with(t[t$site==s & t$colony_id==c,], points(c(1,2,3), c(mean(thalli[age=="age1"]), mean(thalli[age=="age3"]), mean(thalli[age=="age5"])),pch=15, cex=1.5))
    with(t[t$site==s & t$colony_id==c,], boxplot(Yp~age, add=TRUE, border="orange", cex=1.75, lwd=1.75, cex.lab=1.5, main=paste(s," - ",c,sep="")))
  }
}
par(mfrow=c(1,1))
