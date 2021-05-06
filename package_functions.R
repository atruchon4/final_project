CleanData <- function(phy, data) {geiger::treedata(phy,data)
  #treedata() in Geiger is probably my favorite function in R.
}

VisualizeData <- function(phy, data) {
  plot(phy, cex = 0.5)
}

host_hist <- function(tree, data) {
  currenthosts <- data[ ,6:9]
  row.names(currenthosts) <- data[ ,1]
  ancestralhosts <- data[ ,2:5]
  row.names(ancestralhosts) <- data[ ,1]
  currenthosttrait <- NULL
  currenthosttrait[currenthosts[,1]==1] = 0
  currenthosttrait[currenthosts[,2]==1] = 1
  currenthosttrait[currenthosts[,3]==1] = 2
  names(currenthosttrait) <- data[ ,1]
  cleaned <- CleanData(tree, currenthosttrait)
  phy <- cleaned$phy
  dat <- cleaned$data
  phydat <- phangorn::phyDat(cleaned$data, type = "USER", levels = c(0, 1, 2))
  anc.p = phangorn::ancestral.pars(cleaned$phy, phydat)
  anc.m = ancestral.pml(pml(cleaned$phy, phydat), type = "ml")
  lik <- list("phy" = phy, "dat" = dat, "anc.p" = anc.p, "anc.m" = anc.m)
  return(lik)
}

RunEachDiscreteFit <- function(model, phy, data)
{
  discretefit <- fitDiscrete(phy, data, model)
  return(discretefit)
}

FittingModelsDiscrete <- function(phy, data){
  models <- c("SYM", "ARD", "ER", "meristic")
  results <- lapply(models, RunEachDiscreteFit, phy, data)
  aicc1 <- results[[1]]$opt$aicc
  aicc2 <- results[[2]]$opt$aicc
  aicc3 <- results[[3]]$opt$aicc
  aicc4 <- results[[4]]$opt$aicc
  aicc.values <- c(aicc1, aicc2, aicc3, aicc4)
  names(aicc.values)<-models
  aicc.values<-aicc.values-min(aicc.values)
  return(aicc.values)
}

corHMMsepratemodels <- function(matrix, phy, data){
  corHMM(phy,data,rate.cat=1,rate.mat=matrix,node.states="marginal")
}