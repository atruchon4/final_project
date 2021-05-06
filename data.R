virphy <- read.nexus("virus_tree.nex")
host <- read.csv("host_type.csv")
currenthosts <- host[ ,6:9]
row.names(currenthosts) <- host[ ,1]
ancestralhosts <- host[ ,2:5]
row.names(ancestralhosts) <- host[ ,1]
currenthosttrait <- NULL

names(currenthosttrait) <- host[ ,1]
currenthosttrait[currenthosts[,1]==1] = 1
currenthosttrait[currenthosts[,2]==1] = 2
currenthosttrait[currenthosts[,3]==1] = 3


cleaned <- CleanData(virphy, currenthosttrait)
cleaned2 <- CleanData(virphy, ancestralhosttrait)
phydat <- phangorn::phyDat(cleaned$data, type = "USER", levels = c(1, 2, 3))
phydat2 <- phangorn::phyDat(cleaned2$data, type = "USER", levels = c("A", "PLM", "ALM", "AP", "AL", "AM"))
anc.p = phangorn::ancestral.pars(cleaned$phy, phydat)
anc.ml = ancestral.pml(pml(cleaned$phy, phydat), type = "ml")
par(mfcol = (c(1,2)))
plotAnc(cleaned$phy, anc.p, 1)
plotAnc(cleaned$phy, anc.ml, 1)

fitDiscrete(cleaned$phy, cleaned2$data, model = "ER")
fitDiscrete(cleaned$phy, cleaned2$data, model = "SYM")
fitDiscrete(cleaned$phy, cleaned2$data, model = "ARD")
fitDiscrete(cleaned$phy, cleaned2$data, model = "meristic")

models <- c("SYM", "ARD", "ER", "meristic")
results <- lapply(models, RunEachDiscreteFit, cleaned$phy, cleaned$data)

aicc1 <- results[[1]]$opt$aicc
aicc2 <- results[[2]]$opt$aicc
aicc3 <- results[[3]]$opt$aicc
aicc4 <- results[[4]]$opt$aicc
names(aicc.values)<-models
aicc.values<-aicc.values-min(aicc.values)
aicc.values
# Best model for each system is the meristic model, suggesting that transition from trait to trait is not universally distributing and requires more of an island hopping method.

plotSimmap(make.simmap(virphy, currenthosttrait), pts=FALSE, fsize=0.8)
plotSimmap(make.simmap(virphy, ancestralhosttrait), pts=FALSE, fsize=0.8)

ancestralhosttrait <- NULL
ancestralhosttrait <- (c("A", "A", "A", "AL", "AL", "AL", "PLM", "PLM", "PLM", "PLM", "PLM", "ALM", "AP", "AP", "ALM", "AM", "A", "PLM", "A", "A", "A"))
ancestralhosttrait[ancestralhosts[,3]==1] = 0
ancestralhosttrait[ancestralhosts[,3]==0] = 0

ancestralhosttrait <- (c(1, 1, 1, 1, 1, 1, 6, 6, 6, 6, 5, 5, 4, 4, 5, 5, 1, 6, 1, 1, 1))


# 1 Algae
# 2 Protozoa
# 3 Animal
# 4 Algae + Protozoa
# 5 Algae + Animal
# 6 Animal + Protozoa
# 7 Algae + Animal + Protozoa
# 2 Algae + Protozoa
# 3 Algae + Animal
# 4 Algae + Plant
# 5 Algae + Protozoa + Plant
# 6 Algae + Protozoa + Animal
# 7 Algae + Plant + Animal
# 8 Algae + Protozoa + Animal + Plant
# 9 2 Protozoa
# 10 6 Protozoa + Animal
# 11 5 Protozoa + Plant
# 12 8 Protozoa + Plant + Animal
# 13 7 Plant + Animal
# 14 3 Plant
# 15 4 Animal

anestry_animal <- (c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 1, 2, 1, 4, 3))
hostdf <- data.frame(host[ ,1], currenthosttrait, ancestralhosttrait)

rate.mat.9<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=10, model="ARD")
print(rate.mat.9)
pp.ard<-corHMM(virphy,hostdf[,c(1,3)],rate.cat=2,rate.mat=rate.mat.ard.4state,node.states="marginal")

rate.mat.gtr.4state<-rate.mat.ard.4state
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(1,4))
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(2,6))
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(3,8))
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(4,6))
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(5,7))
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(6,7))
print(rate.mat.gtr.4state)

print(rayDISC(virphy, hostdf , ntraits=1, rate.mat= rate.mat.gtr.4state, node.states="marginal", model="ARD"))

rate.mat.ard.24state<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=24, model="ARD")

rate.mat.ard.6state<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=6, model="ARD")

pp.6.state<-corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.fixed,node.states="marginal")

rate.mat.ard.6state.fixed <- rate.par.drop(rate.mat.ard.6state, c(2,5,8,9,10,11,14,17,19, 20, 22, 24, 25, 26, 27, 28))

rate.mat.ard.PLM <- rate.par.drop(rate.mat.ard.6state, c(11,14,26, 27))

pp.6.state.fixed<-corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.fixed,node.states="marginal")

pp.6.state.PLM<-corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.PLM,node.states="marginal")

rate.mat.ard.6state.loose <- rate.par.drop(rate.mat.ard.6state, c(10,11,14,20,25,26, 27))

pp.6.state.loose<-corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal")

rate.mat.ard.6state.noplant <- rate.par.drop(rate.mat.ard.6state, c(6,8,9,10,11,13,14,20,25,26, 27,29,30))


pp.6.state.noplant<-corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.noplant,node.states="marginal")

matrixnames <- c("fixed", "loose", "PLM", "noplant")
matrices <- c(rate.mat.ard.6state.fixed, rate.mat.ard.6state.loose)
matrixAICc <- c(pp.6.state.fixed$AICc, pp.6.state.loose$AICc, pp.6.state.PLM$AICc, NA)
names(matrixAICc) <- matrixnames

pp.fixed.A <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(1,0,0,0,0,0))
pp.fixed.AL <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(0,1,0,0,0,0))
pp.fixed.ALM <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(0,0,1,0,0,0))
pp.fixed.AM <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(0,0,0,1,0,0))
pp.fixed.AP <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(0,0,0,0,1,0))
pp.fixed.PLM <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(0,0,0,0,0,1))

anstatenames <- c("A", "AL", "ALM", "AM", "AP", "PLM")
rootAICc <- c(pp.fixed.A$AICc, pp.fixed.AL$AICc, pp.fixed.ALM$AICc, pp.fixed.AM$AICc, pp.fixed.AP$AICc,pp.fixed.PLM$AICc)
names(rootAICc) <- anstatenames

rate.mat.ard.4state <- rate.mat.maker(rate.cat=1, hrm = FALSE, ntraits=1, nstates=4, model ="ARD")

rate.mat.ard.3state <- rate.mat.maker(rate.cat=1, hrm = FALSE, ntraits=1, nstates=3, model ="ARD")
rate.mat.ard.3state
rate.mat.3state.drop <- rate.par.drop(rate.mat.ard.3state, c(2,5))
rate.mat.3state.sym <- rate.par.eq(rate.mat.3state.drop, c(1,2))
rate.mat.3state.sym <- rate.par.eq(rate.mat.3state.sym, c(2,3))

pp.cur.ard<-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.ard.3state,node.states="marginal")

pp.cur.drop1 <-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.3state.drop,node.states="marginal")

pp.cur.drop2 <-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.3state.sym,node.states="marginal")

curmatrixnames <- c("ARD", "MER", "SYM")
curstatenames <- c("Animal", "Protozoa", "Algae")

curmatricAICc <- c(pp.cur.ard$AICc, pp.cur.drop1$AICc, pp.cur.drop2$AICc)
names(curmatricAICc) <- curmatrixnames
curmatricAICc

pp.cur.root.an <-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.3state.sym,node.states="marginal", root.p = c(1,0,0))
pp.cur.root.prot <-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.3state.sym,node.states="marginal", root.p = c(0,1,0))
pp.cur.root.al <-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.3state.sym,node.states="marginal", root.p = c(0,0,1))

currootAICc <- c(pp.cur.root.an$AICc, pp.cur.root.prot$AICc, pp.cur.root.al$AICc)
names(currootAICc) <- curstatenames
currootAICc

names(ancestralhosttrait) <- host[ ,1]
