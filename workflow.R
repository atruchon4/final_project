# First we read the virus tree and the host data containing ancestral and current host traits.

virphy <- read.nexus("virus_tree.nex")
host <- read.csv("host_type.csv")

# Visualize the virus tree to make sure there aren't any polytomies. The tree looks approximately normal, many of the ranaviruses are extremely closely related, though the distance between them might just appear so slight because of the wide range of diversity of the other viruses.
VisualizeData(virphy, host)

# Here I use the host_history function to estimate the ancestral hosts based on current infection phenotypes. This is done in both a maximum likelihood method and a parsimony method to estimate the ancestry. The distinction between the two ancestral plots is shown here, but very little can be determined easily. For instance, the maximum likelihood method can determine that the Pandoraviruses ancestor likely was a protist and the Ostreococcus/Micromonas virus ancestor was likely algal. However, little can be established about viruses that aren't closely related.
lik <- host_hist(virphy, host)
currentphy <- lik$phy
currentdat <- lik$dat
current.anc.ml <- lik$anc.m
current.anc.p <- lik$anc.p

# Code for ancestry plots:
# 0 = Animal
# 1 = Protist
# 2 = Algae
plotAnc(currentphy, current.anc.p, 1)
plotAnc(currentphy, current.anc.ml, 1)

# I couldn't figure out how to get a function to work with my input data to do what I wanted, so this is hard-coded for now.
ancestralhosts <- host[ ,2:5]
row.names(ancestralhosts) <- host[ ,1]
ancestralhosttrait <- (c("A", "A", "A", "AL", "AL", "AL", "PLM", "PLM", "PLM", "PLM", "PLM", "ALM", "AP", "AP", "ALM", "AM", "A", "PLM", "A", "A", "A"))
names(ancestralhosttrait) <- host[ ,1]

# The same as above can be stated about inferring ancestral states based on DNA fossils in the host genomes. While the ancestor of many of the mimivirus associated viruses appears to be algal, confidence in this trait disappates over time. Many lineages appear in algal genomes, even the animal infection iridovirus, but there is about equal likelihood in many cases that deep branching ancestors were from the PLM category.

cleaned2 <- CleanData(virphy, ancestralhosttrait)
phydat2 <- phangorn::phyDat(cleaned2$data, type = "USER", levels = c("A", "PLM", "ALM", "AP", "AL", "AM"))
anc.anc.p = phangorn::ancestral.pars(cleaned$phy, phydat2)
anc.anc.ml = ancestral.pml(pml(cleaned$phy, phydat2), type = "ml")

# Key for Ancestry Plots:
# A = Algae
# PLM = Protist + Plant + Animal
# ALM = Alage + Plant + Animal
# AP = Algae + Protist
# AL = Algae + Plant
# AM = Algae + Animal
plotAnc(currentphy, anc.anc.p, 1)
plotAnc(currentphy, anc.anc.ml, 1)

# Simmaps were also plotted to see if any differences appeared from the ancestry plots.
plotSimmap(make.simmap(virphy, currenthosttrait), pts=FALSE, fsize=0.8)
plotSimmap(make.simmap(virphy, ancestralhosttrait), pts=FALSE, fsize=0.8)

# In order to assist in assigning an ancestral host based on this dataset, I decided to define rate matrices for both current and ancestral traits. Predicting the rates of transition could improve on the way these ancestral traits were estimated and might give us a better answer to the question on the origin of giant viruses. To do this first, I tried fitting discrete rate models to our traits.

CurrentModelFit <- FittingModelsDiscrete(currentphy, currentdat)
AncModelFit <- FittingModelsDiscrete(currentphy, cleaned2$data)

# The best model for each system is the meristic model, suggesting that transition from trait to trait is not universally distributing and requires more of an island hopping method. I then had to construct these matrices, and designed a few different approaches as described in protocol.txt.

rate.mat.ard.6state<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=6, model="ARD")

rate.mat.ard.6state.fixed <- rate.par.drop(rate.mat.ard.6state, c(2,5,8,9,10,11,14,17,19, 20, 22, 24, 25, 26, 27, 28))
rate.mat.ard.PLM <- rate.par.drop(rate.mat.ard.6state, c(11,14,26, 27))
pp.6.state.fixed<-corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.fixed,node.states="marginal")
pp.6.state.PLM<-corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.PLM,node.states="marginal")

rate.mat.ard.6state.loose <- rate.par.drop(rate.mat.ard.6state, c(10,11,14,20,25,26, 27))
pp.6.state.loose<-corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal")

rate.mat.ard.6state.noplant <- rate.par.drop(rate.mat.ard.6state, c(6,8,9,10,11,13,14,20,25,26, 27,29,30))
pp.6.state.noplant<-corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.noplant,node.states="marginal")

matrixnames <- c("fixed", "loose", "PLM", "noplant")
matrixAICc <- c(pp.6.state.fixed$AICc, pp.6.state.loose$AICc, pp.6.state.PLM$AICc, NA)
names(matrixAICc) <- matrixnames

matrixAICc

# The best matrix I constructed for this dataset based on AIC was the loose fitted, less restrictive matrix. This makes sense considering we are missing many taxa here and therefore transitions from individual to individual may include several hidden states along the tree. Likewise the no plant model was unsuccessful, though I plan on returning to this with more taxa and perhaps a more robust data set.

pp.6.state.loose

# Based o the predicted rates, it seems as though the common transitions are from ALM to AP, which apparently occurred for the Pandoraviruses, from AM to PLM, which appeared in the transition of iridovirus to ranavirusm and from AL to ALM. It is also possible for several traits to somewhat easily revert to only being found in algae. 

pp.fixed.A <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(1,0,0,0,0,0))
pp.fixed.AL <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(0,1,0,0,0,0))
pp.fixed.ALM <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(0,0,1,0,0,0))
pp.fixed.AM <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(0,0,0,1,0,0))
pp.fixed.AP <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(0,0,0,0,1,0))
pp.fixed.PLM <- corHMM(virphy,hostdf[,c(1,3)],rate.cat=1,rate.mat=rate.mat.ard.6state.loose,node.states="marginal", root.p = c(0,0,0,0,0,1))

# Here we test which root state is the best fit for this model.

anstatenames <- c("A", "AL", "ALM", "AM", "AP", "PLM")
rootAICc <- c(pp.fixed.A$AICc, pp.fixed.AL$AICc, pp.fixed.ALM$AICc, pp.fixed.AM$AICc, pp.fixed.AP$AICc,pp.fixed.PLM$AICc)
names(rootAICc) <- anstatenames

rootAICc
# Here we see that the most likely roots are either AP or AM. Historically it appears very likely that algae has always been one of the most common ancestors of giant viruses, though this method also points to animals as being a likely origin. As claims have been made that giant viruses evolved from smaller, animal infecting viruses into larger viruses that typically infect single celled organisms, this could provide evidence in this direction. 

rate.mat.ard.3state <- rate.mat.maker(rate.cat=1, hrm = FALSE, ntraits=1, nstates=3, model ="ARD")
rate.mat.ard.3state
rate.mat.3state.drop <- rate.par.drop(rate.mat.ard.3state, c(2,5))
rate.mat.3state.sym <- rate.par.eq(rate.mat.3state.drop, c(1,2))
rate.mat.3state.sym <- rate.par.eq(rate.mat.3state.sym, c(2,3))

# The same process is performed with the current states of the host infection in order to observe the ancestry from a different perspective. 

pp.cur.ard<-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.ard.3state,node.states="marginal")
pp.cur.drop1 <-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.3state.drop,node.states="marginal")
pp.cur.drop2 <-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.3state.sym,node.states="marginal")

curmatrixnames <- c("ARD", "MER", "SYM")
curstatenames <- c("Animal", "Protozoa", "Algae")

curmatricAICc <- c(pp.cur.ard$AICc, pp.cur.drop1$AICc, pp.cur.drop2$AICc)
names(curmatricAICc) <- curmatrixnames
curmatricAICc

# This method shows that a symmetrical merristic matrix is the best fit for the current host data. Unlike the ancestral host approach, having symmetrical values for host range may make sense in this case, as its less based on the DNA fossil being eliminated from a host and more based on the inividual virus and what it is able to infect.

pp.cur.drop2

# It appears as though the transition from algal host to protist and vice versa occurs more readily than the transition from animal host to protist. This stands to reason as algae and protists are more similar in structure and phylogeny. Likewise, the lack of transition between algae and animal makes sense as the two taxa are very divergent.

pp.cur.root.an <-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.3state.sym,node.states="marginal", root.p = c(1,0,0))
pp.cur.root.prot <-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.3state.sym,node.states="marginal", root.p = c(0,1,0))
pp.cur.root.al <-corHMM(virphy,hostdf[,c(1,2)],rate.cat=1,rate.mat=rate.mat.3state.sym,node.states="marginal", root.p = c(0,0,1))

currootAICc <- c(pp.cur.root.an$AICc, pp.cur.root.prot$AICc, pp.cur.root.al$AICc)
names(currootAICc) <- curstatenames
currootAICc

# Unlike the ancestral root calculation, this method predicts that the best ancestral host would be a protist. While AP is a possible root from the DNA fossil dataset, it isn't as likely as AM, making it odd that in this case protozoa appear to be better hosts than either animals or algae. In fact, animals are the least likely ancestral host in this case. 

# There are many ways this method could be improved upon given more time, and I think one of the best ways to do this would be to add more taxa. It is perhaps and error to only utilize actual giant viruses when trying to determine their origin, especially if they arose from other viruses. However, I believe it would also be a mistake to include all viruses here, as that may confound the results. A possible solution would be to employ the use of Poxviruses. Like giant viruses, poxviruses are large DNA viruses, and although they have been shown to be phylogenetically close to giant viruses, they are not necessarily within the same clade. 
# Additionally, adding more viruses to this phylogeny will improve the robustness of the data and may allow for the testing of other unique transition rate matrices, like the no plant matrix that was attempted here.
# Finally, employing different methods of data collection may also improve on the robustness. For instance, Eukaryotic DNA is typically present to some extent in giant viruses, though it is typically pretty variable in terms of the Eukaryotic taxa present. For this reason, it becomes hard to treat this Eukaryotic DNA as a discrete variable. However, it is also possible to treat the amount of each Eukaryotic taxa DNA as a continuous variable in relation to the percent within the virus genome. Though I didn't have time to do this for the purposes of this class, it would be interesting to use several continuous traits as a means of fitting predicting where a virus has been. This could also be used as a Bayesian predictor for the age of the virus. Building a more robust dataset overall would improve on these results and may help with these predictions.