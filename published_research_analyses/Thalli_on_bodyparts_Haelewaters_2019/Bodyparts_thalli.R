
############################################################################################################################################
#### comparing thallus counts (on different positions) between the three species (and between distinct populations of the same species) ####
############################################################################################################################################

bb<-read.csv("/media/zozoo/3191D18644760ED0/Tan/PhD/Others/Formi/bodyparts/AA_all_bodyparts_2019_V6.csv", sep=";")
bb.rel<-read.csv("/media/zozoo/3191D18644760ED0/Tan/PhD/Others/Formi/bodyparts/AA_all_bodyparts_2019_V6_RELATIVE.csv", sep=";")

require("stats")
require("ggplot2")
require("lsmeans")
require("Rtsne")

############################################################################################

###----- Absolute thalli number -----

qt.1<-glm(full.thalli~species, data=bb, family="quasipoisson")
summary(qt.1)

# extracting unshown factor-level comparison
rg<-ref.grid(qt.1, data=expand.grid(full.thalli=bb$full.thalli, 
                                    species=unique(bb$species)) )
lsmeans(rg, list(pairwise~species))

# visualising thalli distribution across species
ggplot(data = bb, aes(x = species, y = full.thalli)) + 
  geom_boxplot(lwd=1) + theme_bw() + 
  labs(x="Myrmica species'", y="Total number of thalli") + 
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))


###----- Comparing absolute thalli number on each body part between species -----

# carry out Conover-Iman tests
CV.tests<-list()
for(c in 1:16){
  cn<-colnames(bb)[c+6]
  CV.tests[[cn]]<-conover.test(x = bb[,cn], g = bb$species)
}

# gather significance levels of group comparisons (with P-value adjustment)
sign.diffs<-data.frame(bodypart=0, comparison=0, t=0, P=0, P.adjusted=0)
cr<-1
for(i in 1:length(CV.tests)){
  cit<-CV.tests[[i]]
  i.nm<-names(CV.tests)[i]
  sign.diffs[cr:(cr+2),"bodypart"]<-i.nm
  sign.diffs[cr:(cr+2),"comparison"]<-cit$comparisons
  sign.diffs[cr:(cr+2),"t"]<-cit$T
  sign.diffs[cr:(cr+2),"P"]<-cit$P
  cr<-cr+3
}
sign.diffs$P.adjusted<-p.adjust(sign.diffs$P)

# applying t-SNE* and visualising general absolute thalli distribution over ant bodies
# *t-SNE: t-distributed stochastic neighbour embedding
tsne.vars<-bb[,7:22]
Labels<-bb$SpecLoc
colz<-rainbow(length(unique(bb$SpecLoc)))
names(colz)<-unique(bb$SpecLoc)
tsne <- Rtsne(as.matrix(tsne.vars)^(1/2), dims = 2, perplexity=100, verbose=TRUE, max_iter = 500, check_duplicates=FALSE)
df<-cbind(as.data.frame(tsne$Y), as.character(bb$SpecLoc))
colnames(df)<-c("X", "Y", "SpecLoc")
find_hull <- function(df) df[chull(df$X, df$Y), ]
hulls <- ddply(df, "SpecLoc", find_hull)
ggplot(data = df, aes(x = X, y = Y, colour=SpecLoc, fill = SpecLoc)) +
  geom_point(size=2) +
  labs(x = "X", y = "Y") + 
  theme_bw() + scale_fill_discrete(name="Species (country)", 
                                   breaks=c("rubra.AT", "sabuleti.NL", "scabrinodis.HU", "scabrinodis.NL"), 
                                   labels=c("M. rubra (AT)", "M. sabuleti (NL)", "M. scabrinodis (HU)", "M. scabrinodis (NL)")) +
  theme(legend.text=element_text(size=20), legend.title=element_text(size=20)) 


###----- Relative thalli number -----

###----- Comparing relative thalli number* on each body part between species -----
# *relative thalli number: for each body part the proportion of thalli counted on the given body part
#  is divided by the total number of thalli (counted on all body parts of the given individual), 

# carry out Conover-Iman tests
CV.tests.2<-list()
for(c in 1:16){
  cn<-colnames(bb.rel)[c+6]
  CV.tests.2[[cn]]<-conover.test(x = bb.rel[,cn], g = bb.rel$species)
}

# gather significance levels of group comparisons (with P-value adjustment)
sign.diffs.REL<-data.frame(bodypart=0, comparison=0, t=0, P=0, P.adjusted=0)
cr<-1
for(i in 1:length(CV.tests.2)){
  cit<-CV.tests.2[[i]]
  i.nm<-names(CV.tests.2)[i]
  sign.diffs.REL[cr:(cr+2),"bodypart"]<-i.nm
  sign.diffs.REL[cr:(cr+2),"comparison"]<-cit$comparisons
  sign.diffs.REL[cr:(cr+2),"t"]<-cit$T
  sign.diffs.REL[cr:(cr+2),"P"]<-cit$P
  cr<-cr+3
}
sign.diffs.REL$P.adjusted<-p.adjust(sign.diffs.REL$P)

# applying t-SNE and visualising general relative thalli distribution over ant bodies
tsne.vars.rel<-bb.rel[,7:22]
Labels<-bb.rel$SpecLoc
colz<-rainbow(length(unique(bb.rel$SpecLoc)))
names(colz)<-unique(bb.rel$SpecLoc)
tsne <- Rtsne(as.matrix(tsne.vars.rel)^(1/2), dims = 2, perplexity=100, verbose=TRUE, max_iter = 500, check_duplicates=FALSE)
df<-cbind(as.data.frame(tsne$Y), as.character(bb.rel$SpecLoc))
colnames(df)<-c("X", "Y", "SpecLoc")
find_hull <- function(df) df[chull(df$X, df$Y), ]
hulls <- ddply(df, "SpecLoc", find_hull)
ggplot(data = df, aes(x = X, y = Y, colour=SpecLoc, fill = SpecLoc)) +
  geom_point(size=3) +
  labs(x = "X", y = "Y") + 
  theme_bw()

###----- Infering on the potential origin of infection -----

bps = data.frame(species=0, bodypart=0, min.rel.thalli=0, max.rel.thalli=0)
cc = 1
for(s in unique(bb.rel$species)){
  for(n in colnames(bb.rel)[7:22]){
    bps[cc, "species"] = s
    bps[cc, "bodypart"] = n
    bps[cc, "min.rel.thalli"] = min(bb.rel[bb.rel$species==s, n])
    bps[cc, "max.rel.thalli"] = max(bb.rel[bb.rel$species==s, n])
    cc = cc + 1
  }
}

bps$bodypart[bps$species=="rubra" & bps$min.rel.thalli>0 & bps$max.rel.thalli[bps$species=="rubra"]==max(bps$max.rel.thalli[bps$species=="rubra"])]
bps$bodypart[bps$species=="sabuleti" & bps$min.rel.thalli>0 & bps$max.rel.thalli[bps$species=="sabuleti"]==max(bps$max.rel.thalli[bps$species=="sabuleti"])]

# no bodyparts in M. scabrinodis on which minimal relative thalli value is larger than zero; exclude this criterion
bps$bodypart[bps$species=="scabrinodis" & bps$max.rel.thalli[bps$species=="scabrinodis"]==max(bps$max.rel.thalli[bps$species=="scabrinodis"])]


# fitting binomial regression models to test whether relative thalli number decreases 
# with the increase of total thalli number
bb.rel$full.thalli.centeredSqrt<-sqrt(bb.rel$full.thalli)-mean(sqrt(bb.rel$full.thalli))
bb.rel$HeadDors.sqrt<-sqrt(bb.rel$Head..dorsal)

# for M. rubra and M. sabuleti (dorsal head part)
strt.m1<-glm(HeadDors.sqrt ~ full.thalli.centeredSqrt * SpecLoc, data = bb.rel[bb.rel$infected==1,], family="quasibinomial")
summary(strt.m1)
sim_slopes(model = strt.m1, pred = full.thalli.centeredSqrt, modx = SpecLoc, johnson_neyman = FALSE, centered = "none", digits = 3)
rg<-ref.grid(strt.m1, data=expand.grid(HeadDors.sqrt=bb.rel$HeadDors.sqrt, 
                                       full.thalli.centeredSqrt=bb.rel$full.thalli.centeredSqrt,
                                       SpecLoc=unique(bb.rel$SpecLoc)) )
lsmeans(rg, list(pairwise~SpecLoc))


# for M. scabrinodis (gaster-tergites)
bb.rel$GasterTergites.sqrt<-sqrt(bb.rel$Gaster..tergites)
strt.m2<-glm(GasterTergites.sqrt ~ full.thalli.centeredSqrt * SpecLoc, data = bb.rel[bb.rel$infected==1,], family="quasibinomial")
summary(strt.m2)
sim_slopes(model = strt.m2, pred = full.thalli.centeredSqrt, modx = SpecLoc, johnson_neyman = FALSE, centered = "none", digits = 3)

# when excluding those observations that are zero (18 out of 354: 5%) 
# the pattern is significant and is as expected
strt.m2.2<-glm(GasterTergites.sqrt ~ full.thalli.centeredSqrt * SpecLoc, data = bb.rel[bb.rel$infected==1 & bb.rel$Gaster..tergites!=0,], family="quasibinomial")
summary(strt.m2.2)

sim_slopes(model = strt.m2.2, pred = full.thalli.centeredSqrt, modx = SpecLoc, johnson_neyman = FALSE, centered = "none", digits = 3)
rg<-ref.grid(strt.m2.2, data=expand.grid(GasterTergites.sqrt=bb.rel$GasterTergites.sqrt, 
                                       full.thalli.centeredSqrt=bb.rel$full.thalli.centeredSqrt,
                                       SpecLoc=unique(bb.rel$SpecLoc)) )
lsmeans(rg, list(pairwise~SpecLoc))


