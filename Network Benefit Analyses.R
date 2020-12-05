##### Benefit Analyses

# This code estimates the dispersal benefit ratio for seed
# and seedling survival experiments, then assesses relationships
# between network metrics and the dispersal benefit ratio, and
# finally presents figures showing these relationships as well as
# the relationship between species strength and the pulp:seed ratio.

#load libraries & source files
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("MCMCglmm","mcmcplots","coda","ggplot2","plyr","boot")

ipak(packages)


topwd <- getwd()
setwd(paste(topwd,"Input data/Plant benefit data",sep="/"))


##### Field Seedlings

field_seedlings <- read.csv("field_seedlings.csv",row.names=1,header=T)

###set priors
priors.sdlng <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V=1, nu=1, alpha.mu=0, alpha.V= 10^3),
                                G2 = list(V=1, nu=1, alpha.mu=0, alpha.V= 10^3)))

### Run Models
set.seed(11)
sdlng.mod <- MCMCglmm(cbind(numalive, numseedplant-numalive) ~ dist*species + centavgopen*species, 
                      random = ~island+site, 
                      family = "multinomial2", 
                      data = field_seedlings, 
                      prior = priors.sdlng, 
                      verbose = F, 
                      pr=TRUE, #store random effects
                      pl = TRUE, #store latent variable
                      burnin = 10000,
                      nitt = 300000,
                      thin = 50) 
# summary(sdlng.mod) #
# raftery.diag(sdlng.mod) #3746 samples
# autocorr(sdlng.mod$Sol) #
# autocorr(sdlng.mod$VCV) #
# allChains <- as.mcmc(cbind(sdlng.mod$Sol,sdlng.mod$VCV))  
# plot(allChains) #
# 
# #check for overdispersion
# mcmc.units <- colMeans(sdlng.mod$Liab) - predict(sdlng.mod, marginal = NULL, type = "terms") #fitted minus observed 
# par(mfrow=c(1,2))
# qqnorm(mcmc.units) ; qqline(mcmc.units)
# hist(mcmc.units, col = grey(.6), main = "Over-dispersion")
# 
# #Check distributional assumptions of random effects
# mcmc.int <- colMeans(sdlng.mod$Sol[,13:29]) # one value for each random effect parameter 
# qqnorm(mcmc.int) ; qqline(mcmc.int)
# hist(mcmc.int, col = grey(.6), main = "Intercept")
# 
# #Check to see if everything is significant.
# cbind(B = colMeans(sdlng.mod$Sol[,1:(length(colnames(sdlng.mod$Sol))-17)]), CI = HPDinterval(sdlng.mod$Sol[,1:(length(colnames(sdlng.mod$Sol))-17)]))
# 
### Graphing

#functions

MyLinesStuff <- function(x){
  OUT <- matrix(nrow = nrow(x), ncol=4) 
  for(i in 1:nrow(x)){
    xi <- x[i,]  
    OUT[i,3:4] <- quantile(xi, probs = c(0.025, 0.975))
    OUT[i,1] <- mean(xi)
    OUT[i,2] <- sd(xi)
  }
  colnames(OUT) <- c("mean", "se", "lowcri", "uppcri")
  OUT
}


# 
MyRatios<- function(x){
  OUT <- numeric(0)
  for(i in 1:length(levels(x$species))) {
    sppnear<-x[x$species==levels(x$species)[i] & x$dist=="near","Pi"]
    sppfar<-x[x$species==levels(x$species)[i] & x$dist=="far","Pi"]
    z<- log(sppfar / sppnear)
    OUT<-rbind(OUT,z)
  }
  OUT
}

#Predict the mean, ci's 5800 times (1 for each iteration)
#1) Get the  betas
beta<-sdlng.mod$Sol[,1:(length(colnames(sdlng.mod$Sol))-17)] #5800 values

#2) Create a grid of covariate values 
range(field_seedlings$centavgopen)
preddata <- with(field_seedlings,expand.grid(dist=levels(dist), 
                                             species=levels(species), 
                                             centavgopen=c(quantile(field_seedlings$centavgopen,0.1),
                                                           quantile(field_seedlings$centavgopen,0.5),
                                                           quantile(field_seedlings$centavgopen,0.9))))

#3) For this we need the X on a grid
X <- model.matrix(~ dist*species+centavgopen*species, data=preddata)

#4) #Calculate X * beta for each MCMC iteration, 4500 times, to get 4500 predictions of propsurvival for each distance and species combo!
eta.mcmc<-(X %*% t(beta))
mu.mcmc<-inv.logit(eta.mcmc) #logit link
dim(mu.mcmc)

#5) Calculate the mean, se, and 2.5% and 97.5% credible intervals for each of the species:dist combinations. x is a matrix of predicted values for the grid we created earlier (each species, near and far). the mean, se, and quantiles are calculated on each row. 

L <- MyLinesStuff(eta.mcmc)
L #Posterior mean, se and 95% CI for each fixed effect parameter

#6) combine with preddata, so can see which dist/spp each row is
predmcmc<-cbind(preddata,L) #new dataframe with posterior mean, se, and CI's

#7) get mean and se on original data scale (not link)
predmcmc$predpropsurv<-inv.logit(predmcmc$mean)
predmcmc$uci<-inv.logit(predmcmc$mean + 1.96 * predmcmc$se)
predmcmc$lci<-inv.logit(predmcmc$mean - 1.96 * predmcmc$se)

# #8) Graph
# limits<- aes(ymin=lci, ymax=uci) #using se's
# limits<-aes(ymin=inv.logit(lowcri), ymax=inv.logit(uppcri)) #using credible intervals
# #get very similar results
# 
# ggplot(predmcmc[predmcmc$centavgopen==0,], aes(x=dist, y=predpropsurv))+
#   geom_point()+
#   geom_errorbar(limits)+
#   facet_grid(.~species)

#9) Calculate differences between near and far, and create confidence intervals for these differences
preddata$eta.mcmc<-(X %*% t(beta))
preddata$Pi <- exp(preddata$eta.mcmc) / (1 + exp(preddata$eta.mcmc)) #converted to probabilities (i.e. on the scale of the predictor 0 to 1) 

diffs<-MyRatios(preddata[preddata$centavgopen==quantile(field_seedlings$centavgopen,0.5),])

diffs.sum<-MyLinesStuff(diffs)
# I think species is just the species list
species <- levels(field_seedlings$species)

species2<-as.data.frame(species)
sppdistdiffs<-cbind(species2, diffs.sum)

sppdistdiffs$species<-factor(sppdistdiffs$species, levels=c("aglaia","psychotria", "papaya", "premna"), ordered=F)
limits<- aes(ymin=lowcri, ymax=uppcri)
# ggplot(sppdistdiffs, aes(species,mean))+
#   geom_point()+
#   geom_errorbar(limits, width=0.35)+
#   theme_bw()+
#   geom_abline(intercept=0, slope=0)+ 
#   ylab("Dispersal Benefit Ratio Difference in survival near - far at mean light levels")+
#   xlab("")+
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),    text = element_text(size=20))+
#   scale_x_discrete(labels=c("Aglaia", "Psychotria","Carica", "Premna"))




##### Seedlings Benefit of light (from same model) #####


MyRatiosLight<- function(x){
  OUT <- numeric(0)
  for(i in 1:length(levels(x$species))) {
    sppfarll<-x[x$species==levels(x$species)[i] & x$dist=="far" & x$centavgopen == quantile(field_seedlings$centavgopen,0.1),"Pi"]
    sppfarhl<-x[x$species==levels(x$species)[i] & x$dist=="far" & x$centavgopen == quantile(field_seedlings$centavgopen,0.9),"Pi"]
    z<- log(sppfarhl / sppfarll)
    OUT<-rbind(OUT,z)
  }
  OUT
}

light.diffs<-MyRatiosLight(preddata)
light.diffs.sum<-MyLinesStuff(light.diffs)

spplightdiffs<-cbind(species2, light.diffs.sum)

spplightdiffs$species<-factor(spplightdiffs$species, levels=c("aglaia","papaya","premna","psychotria"), ordered=F)
limits<- aes(ymin=lowcri, ymax=uppcri)
# ggplot(spplightdiffs, aes(species, mean))+
#   geom_point()+
#   geom_errorbar(limits, width=0.35)+
#   theme_bw()+
#   geom_abline(intercept=0, slope=0)+ 
#   ylab("Dispersal Benefit Ratio (High vs low light)")+
#   xlab("")+
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),    text = element_text(size=20))+
#   scale_x_discrete(labels=c("Aglaia","Papaya", "Premna","Psychotria"))+ geom_hline(yintercept=1)






##### Field Seedlings

field_seeds <- read.csv("field_seeds.csv",row.names=1,header=T)

priorsseeds <- list(R = list(V = 1, nu = 0.002),
                    G = list(G1 = list(V=1, nu=1, alpha.mu=0, alpha.V= 10^3)))


### Run Models
set.seed(11)
seed.mod <- MCMCglmm(cbind(seedlings, seed_input-seedlings) ~ near_far*species + species*dens, 
                     random = ~sitepoint, 
                     family = "multinomial2", 
                     data = field_seeds, 
                     prior = priorsseeds, 
                     verbose = F, 
                     pr=TRUE, #store random effects
                     pl = TRUE, #store latent variable
                     burnin = 10000,
                     nitt = 300000,
                     thin = 50) 
summary(seed.mod) #
# raftery.diag(seed.mod) #3746 samples
# autocorr(seed.mod$Sol) #
# autocorr(seed.mod$VCV) #
# allChains <- as.mcmc(cbind(seed.mod$Sol,seed.mod$VCV))  
# plot(allChains) #

# #check for overdispersion
# mcmc.units <- colMeans(seed.mod$Liab) - predict(seed.mod, marginal = NULL, type = "terms") #fitted minus observed 
# par(mfrow=c(1,2))
# qqnorm(mcmc.units) ; qqline(mcmc.units)
# hist(mcmc.units, col = grey(.6), main = "Over-dispersion")
# 
# #Check distributional assumptions of random effects
# mcmc.int <- colMeans(seed.mod$Sol[,(7:59)]) # one value for each random effect parameter 
# qqnorm(mcmc.int) ; qqline(mcmc.int)
# hist(mcmc.int, col = grey(.6), main = "Intercept")

#Check to see if everything is significant.
cbind(B = colMeans(seed.mod$Sol[,1:(length(colnames(seed.mod$Sol))-3)]), CI = HPDinterval(seed.mod$Sol[,1:(length(colnames(seed.mod$Sol))-3)]))

# Need slightly different version of this function
MyRatios<- function(x){
  OUT <- numeric(0)
  for(i in 1:length(levels(x$species))) {
    sppnear<-x[x$species==levels(x$species)[i] & x$near_far=="near","Pi"]
    sppfar<-x[x$species==levels(x$species)[i] & x$near_far=="far","Pi"]
    z<- log(sppfar / sppnear)
    OUT<-rbind(OUT,z)
  }
  OUT
}



### Effect of distance on survival

#make predictions on original data (so only get predictions for actual light levels measured)
seed.modpred <- predict(seed.mod, interval="confidence", type="response", level=0.95, marginal=NULL) 
colnames(seed.modpred) <- paste("seed.mod",c("fit","lwr","upr"),sep=".")
seed.modmod<-cbind(field_seeds, seed.modpred)
ddply(seed.modmod, .(species, near_far), summarize, mean=mean(seed.mod.fit))


### create predictions "by hand"

#Predict the mean, ci's 5800 times (1 for each iteration)
#1) Get the  betas
beta.ss<-seed.mod$Sol[,1:(length(colnames(seed.mod$Sol))-54)] #5800 values

#2) Create a grid of covariate values 
range(field_seeds$dens) #86.75 to 97.25
mean(field_seeds$dens) #92.29
preddata.ss <- with(field_seeds,expand.grid(near_far=levels(near_far), 
                                            species=levels(species), 
                                            dens=c(quantile(field_seeds$dens,0.1),
                                                   quantile(field_seeds$dens,0.5),
                                                   quantile(field_seeds$dens,0.9))))

#3) For this we need the X.ss on a grid
X.ss <- model.matrix(~ near_far*species+species*dens, data=preddata.ss)

#4) #Calculate X.ss * beta.ss for each MCMC iteration, 4500 times, to get 4500 predictions of propsurvival for each distance and species combo!
eta.mcmc.ss<-(X.ss %*% t(beta.ss))
mu.mcmc.ss<-inv.logit(eta.mcmc.ss) #logit link
dim(mu.mcmc.ss)

#5) Calculate the mean, se, and 2.5% and 97.5% credible intervals for each of the species:dist combinations. x is a matrix of predicted values for the grid we created earlier (each species, near and far). the mean, se, and quantiles are calculated on each row. 

L.ss <- MyLinesStuff(eta.mcmc.ss)
L.ss #Posterior mean, se and 95% CI for each fixed effect parameter

#6) combine with preddata.ss, so can see which dist/spp each row is
predmcmc.ss<-cbind(preddata.ss,L.ss) #new dataframe with posterior mean, se, and CI's

#7) get mean and se on original data scale (not link)
predmcmc.ss$predpropsurv<-inv.logit(predmcmc.ss$mean)
predmcmc.ss$uci<-inv.logit(predmcmc.ss$mean + 1.96 * predmcmc.ss$se)
predmcmc.ss$lci<-inv.logit(predmcmc.ss$mean - 1.96 * predmcmc.ss$se)

#8) Graph
limits<- aes(ymin=lci, ymax=uci) #using se's
limits<-aes(ymin=inv.logit(lowcri), ymax=inv.logit(uppcri)) #using credible intervals
#get very similar results

# ggplot(predmcmc.ss[predmcmc.ss$dens==quantile(field_seeds$dens,0.5),], aes(x=near_far, y=predpropsurv))+
#   geom_point()+
#   geom_errorbar(limits)+
#   facet_grid(.~species)

#9) Calculate differences between near and far, and create confidence intervals for these differences
preddata.ss$eta.mcmc.ss<-(X.ss %*% t(beta.ss))
preddata.ss$Pi <- exp(preddata.ss$eta.mcmc.ss) / (1 + exp(preddata.ss$eta.mcmc.ss)) #converted to probabilities (i.e. on the scale of the predictor 0 to 1) 

diffs.ss<-MyRatios(preddata.ss[preddata.ss$dens==quantile(field_seeds$dens,0.5),])

diffs.ss.sum<-MyLinesStuff(diffs.ss)
# I think species is just the species list
species <- levels(field_seeds$species)

species2<-as.data.frame(species)
sppdistdiffs.ss<-cbind(species2, diffs.ss.sum)

sppdistdiffs.ss$species<-factor(sppdistdiffs.ss$species, levels=c("aglaia","premna"), ordered=F)
limits<- aes(ymin=lowcri, ymax=uppcri)
# ggplot(sppdistdiffs.ss, aes(species, mean))+
#   geom_point()+
#   geom_errorbar(limits, width=0.35)+
#   theme_bw()+
#   geom_abline(intercept=0, slope=0)+ 
#   ylab("Dispersal Benefit Ratio Difference in survival near - far at mean light levels")+
#   xlab("")+
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),    text = element_text(size=20))+
#   scale_x_discrete(labels=c("Aglaia", "Premna"))+ geom_hline(yintercept=1)




##### Seeds Benefit of light (from same model) #####


#Need slightly different version of this function
MyRatiosLight<- function(x){
  OUT <- numeric(0)
  for(i in 1:length(levels(x$species))) {
    sppfarll<-x[x$species==levels(x$species)[i] & x$near_far=="far" & x$dens == quantile(field_seeds$dens,0.9),"Pi"]
    sppfarhl<-x[x$species==levels(x$species)[i] & x$near_far=="far" & x$dens == quantile(field_seeds$dens,0.1),"Pi"]
    z<- log(sppfarhl / sppfarll)
    OUT<-rbind(OUT,z)
  }
  OUT
}


light.diffs.ss<-MyRatiosLight(preddata.ss)

light.diffs.ss.sum<-MyLinesStuff(light.diffs.ss)

spplightdiffs.ss<-cbind(species2, light.diffs.ss.sum)

spplightdiffs.ss$species<-factor(spplightdiffs.ss$species, levels=c("aglaia","premna"), ordered=F)
limits<- aes(ymin=lowcri, ymax=uppcri)
# ggplot(spplightdiffs.ss, aes(species, mean))+
#   geom_point()+
#   geom_errorbar(limits, width=0.35)+
#   theme_bw()+
#   geom_abline(intercept=0, slope=0)+ 
#   ylab("Dispersal Benefit Ratio (High vs low light)")+
#   xlab("")+
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),    text = element_text(size=20))+
#   scale_x_discrete(labels=c("Aglaia","Premna"))



##### Gap survival

gap_seeds <- read.csv("gap_seeds.csv",row.names=1,header=T)

#set priors
priorsgap <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V=1, nu=1, alpha.mu=0, alpha.V= 10^3)))



gap_seeds$gap_nongap <- as.factor(gap_seeds$gap_nongap)
gap_seeds$species <- as.factor(gap_seeds$species)
gap_seeds$fence_site <- as.factor(paste(gap_seeds$fence,gap_seeds$site,sep="_"))

### Run Models

set.seed(11)
gap.mod <- MCMCglmm(cbind(germ, seed_input-germ) ~ gap_nongap*species, 
                    random = ~fence_site, 
                    family = "multinomial2", 
                    data = gap_seeds, 
                    prior = priorsgap, 
                    verbose = F, 
                    pr=TRUE, #store random effects
                    pl = TRUE, #store latent variable
                    burnin = 10000,
                    nitt = 300000,
                    thin = 50) 
summary(gap.mod) #
# raftery.diag(gap.mod) #3746 samples
# autocorr(gap.mod$Sol) #
# autocorr(gap.mod$VCV) #
# allChains <- as.mcmc(cbind(gap.mod$Sol,gap.mod$VCV))  
# plot(allChains) #

# #check for overdispersion
# mcmc.units <- colMeans(gap.mod$Liab) - predict(gap.mod, marginal = NULL, type = "terms") #fitted minus observed 
# par(mfrow=c(1,2))
# qqnorm(mcmc.units) ; qqline(mcmc.units)
# hist(mcmc.units, col = grey(.6), main = "Over-dispersion")
# 
# #Check distributional assumptions of random effects
# mcmc.int <- colMeans(gap.mod$Sol[,c(7:18)]) # one value for each random effect parameter 
# qqnorm(mcmc.int) ; qqline(mcmc.int)
# hist(mcmc.int, col = grey(.6), main = "Intercept")
# 
# #Check to see if everything is significant.
# cbind(B = colMeans(gap.mod$Sol[,1:(length(colnames(gap.mod$Sol))-3)]), CI = HPDinterval(gap.mod$Sol[,1:(length(colnames(gap.mod$Sol))-3)]))


MyLinesStuff <- function(x){
  OUT <- matrix(nrow = nrow(x), ncol=4) 
  for(i in 1:nrow(x)){
    xi <- x[i,]  
    OUT[i,3:4] <- quantile(xi, probs = c(0.025, 0.975))
    OUT[i,1] <- mean(xi)
    OUT[i,2] <- sd(xi)
  }
  colnames(OUT) <- c("mean", "se", "lowcri", "uppcri")
  OUT
}

# Need slightly different version of this function
MyRatios<- function(x){
  OUT <- numeric(0)
  for(i in 1:length(levels(x$species))) {
    sppnongap<-x[x$species==levels(x$species)[i] & x$gap_nongap=="nongap","Pi"]
    sppgap<-x[x$species==levels(x$species)[i] & x$gap_nongap=="gap","Pi"]
    z<- log(sppgap / sppnongap)
    OUT<-rbind(OUT,z)
  }
  OUT
}

#make predictions on original data (so only get predictions for actual light levels measured)
gap.modpred <- predict(gap.mod, interval="confidence", type="response", level=0.95, marginal=NULL) 
colnames(gap.modpred) <- paste("gap.mod",c("fit","lwr","upr"),sep=".")
gap.modmod<-cbind(gap_seeds, gap.modpred)
ddply(gap.modmod, .(species, gap_nongap), summarize, mean=mean(gap.mod.fit))


### create predictions "by hand"

#Predict the mean, ci's 5800 times (1 for each iteration)
#1) Get the  betas
beta.gng<-gap.mod$Sol[,1:(length(colnames(gap.mod$Sol))-12)] #5800 values

#2) Create a grid of covariate values. Don't need this first part!
#range(ss.dat$dens) #86 to 97.25
#mean(ss.dat$dens) #92.27
preddata.gng <- with(gap_seeds,expand.grid(gap_nongap=levels(gap_nongap), species=levels(species)))

#3) For this we need the X.gng on a grid
X.gng <- model.matrix(~ gap_nongap*species, data=preddata.gng)

#4) #Calculate X.gng * beta.gng for each MCMC iteration, 4500 times, to get 4500 predictions of propsurvival for each distance and species combo!
eta.mcmc.gng<-(X.gng %*% t(beta.gng))
mu.mcmc.gng<-inv.logit(eta.mcmc.gng) #logit link
dim(mu.mcmc.gng)

#5) Calculate the mean, se, and 2.5% and 97.5% credible intervals for each of the species:dist combinations. x is a matrix of predicted values for the grid we created earlier (each species, near and far). the mean, se, and quantiles are calculated on each row. 

L.gng <- MyLinesStuff(eta.mcmc.gng)
L.gng #Posterior mean, se and 95% CI for each fixed effect parameter

#6) combine with preddata.gng, so can see which dist/spp each row is
predmcmc.gng<-cbind(preddata.gng,L.gng) #new dataframe with posterior mean, se, and CI's

#7) get mean and se on original data scale (not link)
predmcmc.gng$predpropsurv<-inv.logit(predmcmc.gng$mean)
predmcmc.gng$uci<-inv.logit(predmcmc.gng$mean + 1.96 * predmcmc.gng$se)
predmcmc.gng$lci<-inv.logit(predmcmc.gng$mean - 1.96 * predmcmc.gng$se)

#8) Graph
limits<- aes(ymin=lci, ymax=uci) #using se's
limits<-aes(ymin=inv.logit(lowcri), ymax=inv.logit(uppcri)) #using credible intervals
#get very similar results

# ggplot(predmcmc.gng, aes(x=gap_nongap, y=predpropsurv))+
#   geom_point()+
#   geom_errorbar(limits)+
#   facet_grid(.~species)

#9) Calculate differences between near and far, and create confidence intervals for these differences
preddata.gng$eta.mcmc.gng<-(X.gng %*% t(beta.gng))
preddata.gng$Pi <- exp(preddata.gng$eta.mcmc.gng) / (1 + exp(preddata.gng$eta.mcmc.gng)) #converted to probabilities (i.e. on the scale of the predictor 0 to 1) 

diffs.gng<-MyRatios(preddata.gng)

diffs.gng.sum<-MyLinesStuff(diffs.gng)
# I think species is just the species list
species <- levels(gap_seeds$species)

species2<-as.data.frame(species)
sppdistdiffs.gng<-cbind(species2, diffs.gng.sum)

sppdistdiffs.gng$species<-factor(sppdistdiffs.gng$species, levels=c("aglaia","papaya","premna"), ordered=F)
limits<- aes(ymin=lowcri, ymax=uppcri)
# ggplot(sppdistdiffs.gng, aes(species, mean))+
#   geom_point()+
#   geom_errorbar(limits, width=0.35)+
#   theme_bw()+
#   geom_abline(intercept=0, slope=0)+ 
#   ylab("Dispersal Benefit Ratio Difference in survival near - far at mean light levels")+
#   xlab("")+
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),    text = element_text(size=20))+
#   scale_x_discrete(labels=c("Aglaia","Papaya", "Premna"))+ geom_hline(yintercept=1)






##### Frugivore Handling

handling_seeds <- read.csv("handling_seeds.csv",row.names=1,header=T)


priorshandling <- list(R = list(V = 1, fix = 1),
                       G = list(G1 = list(V=1, nu=1, alpha.mu=0, alpha.V= 10^3),
                                G2 = list(V=1, nu=1, alpha.mu=0, alpha.V= 10^3)))


### Run Models

handling.mod <- MCMCglmm(germ ~ trt*species, 
                         random = ~island+site, 
                         family = "categorical", 
                         data = handling_seeds, 
                         prior = priorshandling, 
                         verbose = F, 
                         pr=TRUE, #store random effects
                         pl = TRUE, #store latent variable
                         burnin = 10000,
                         nitt = 300000,
                         thin = 50) 
summary(handling.mod) #
# raftery.diag(handling.mod) #3746 samples
# autocorr(handling.mod$Sol) #
# autocorr(handling.mod$VCV) #
# allChains <- as.mcmc(cbind(handling.mod$Sol,handling.mod$VCV))  
# plot(allChains) #

# #check for overdispersion
# mcmc.units <- colMeans(handling.mod$Liab) - predict(handling.mod, marginal = NULL, type = "terms") #fitted minus observed 
# par(mfrow=c(1,2))
# qqnorm(mcmc.units) ; qqline(mcmc.units)
# hist(mcmc.units, col = grey(.6), main = "Over-dispersion")
# 
# #Check distributional assumptions of random effects
# mcmc.int <- colMeans(handling.mod$Sol[,c(5:19)]) # one value for each random effect parameter 
# qqnorm(mcmc.int) ; qqline(mcmc.int)
# hist(mcmc.int, col = grey(.6), main = "Intercept")

#Check to see if everything is significant.
cbind(B = colMeans(handling.mod$Sol[,1:(length(colnames(handling.mod$Sol))-3)]), CI = HPDinterval(handling.mod$Sol[,1:(length(colnames(handling.mod$Sol))-3)]))




### Effect of handling on survival

#make predictions on original data (so only get predictions for actual light levels measured)
handling.modpred <- predict(handling.mod, interval="confidence", type="response", level=0.95, marginal=NULL) 
colnames(handling.modpred) <- paste("handling.mod",c("fit","lwr","upr"),sep=".")
handling.modmod<-cbind(handling_seeds, handling.modpred)
ddply(handling.modmod, .(species, trt), summarize, mean=mean(handling.mod.fit))


### create predictions "by hand"


#Predict the mean, ci's 5800 times (1 for each iteration)
#1) Get the  betas
beta.gut<-handling.mod$Sol[,1:4] #5800 values

#2) Create a grid of covariate values 
preddata.gut <- with(handling_seeds,expand.grid(trt=levels(trt), 
                                                species=levels(species)))  #alternative: centavgopen=seq(-4.62, 21.88, length = 30)

#3) For this we need the X.gut on a grid
X.gut <- model.matrix(~ trt*species, data=preddata.gut)

#4) #Calculate X.gut * beta.gut for each MCMC iteration, 4500 times, to get 4500 predictions of propsurvival for each distance and species combo!
eta.mcmc.gut<-(X.gut %*% t(beta.gut))
mu.mcmc.gut<-inv.logit(eta.mcmc.gut) #logit link
dim(mu.mcmc.gut)

#5) Calculate the mean, se, and 2.5% and 97.5% credible intervals for each of the species:dist combinations. x is a matrix of predicted values for the grid we created earlier (each species, near and far). the mean, se, and quantiles are calculated on each row. 

L.gut <- MyLinesStuff(eta.mcmc.gut)
L.gut #Posterior mean, se and 95% CI for each fixed effect parameter

#6) combine with preddata.gut, so can see which dist/species each row is
predmcmc.gut<-cbind(preddata.gut,L.gut) #new dataframe with posterior mean, se, and CI's

#7) get mean and se on original data scale (not link)
predmcmc.gut$predpropsurv<-inv.logit(predmcmc.gut$mean)
predmcmc.gut$uci<-inv.logit(predmcmc.gut$mean + 1.96 * predmcmc.gut$se)
predmcmc.gut$lci<-inv.logit(predmcmc.gut$mean - 1.96 * predmcmc.gut$se)

#8) Graph
limits<- aes(ymin=lci, ymax=uci) #using se's
limits<-aes(ymin=inv.logit(lowcri), ymax=inv.logit(uppcri)) #using credible intervals
#get very similar results

# ggplot(predmcmc.gut, aes(x=trt, y=predpropsurv))+
#   geom_point()+
#   geom_errorbar(limits)+
#   facet_grid(.~species)

#9) Calculate differences between near and far, and create confidence intervals for these differences
preddata.gut$eta.mcmc.gut<-(X.gut %*% t(beta.gut))
preddata.gut$Pi <- exp(preddata.gut$eta.mcmc.gut) / (1 + exp(preddata.gut$eta.mcmc.gut)) #converted to probabilities (i.e. on the scale of the predictor 0 to 1) 

# Need slightly different version of this function
MyRatios<- function(x){
  OUT <- numeric(0)
  for(i in 1:length(levels(x$species))) {
    speciesunhand<-x[x$species==levels(x$species)[i] & x$trt=="UH","Pi"]
    specieshand<-x[x$species==levels(x$species)[i] & x$trt=="HA","Pi"]
    z<- log(specieshand / speciesunhand)
    OUT<-rbind(OUT,z)
  }
  OUT
}

diffs.gut<-MyRatios(preddata.gut)

str(diffs.gut)

diffs.gut.sum<-MyLinesStuff(diffs.gut)

# I think species is just the species list
species <- levels(handling_seeds$species)

species2<-as.data.frame(c("pr","ps"))
colnames(species2)[1] <- "species"
sppdistdiffs.gut<-cbind(species2, diffs.gut.sum)

sppdistdiffs.gut$species<-factor(sppdistdiffs.gut$species, levels=c("pr","ps"), ordered=F)
limits<- aes(ymin=lowcri, ymax=uppcri)
# ggplot(sppdistdiffs.gut, aes(species, mean))+
#   geom_point()+
#   geom_errorbar(limits, width=0.35)+
#   theme_bw()+
#   geom_abline(intercept=0, slope=0)+ 
#   ylab("Dispersal Benefit Ratio Difference in survival near - far at mean light levels")+
#   xlab("")+
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),    text = element_text(size=20))+
#   scale_x_discrete(labels=c("Premna", "Psychotria"))+ geom_hline(yintercept=0)





##### Distance dependence, shadehouse

dd_shadehouse <- read.csv("dd_shadehouse_seeds.csv",row.names=1,header=T)


#set priors
priorsshade <- list(R = list(V = 1, fix = 1),
                    G = list(G1 = list(V=1, nu=1, alpha.mu=0, alpha.V= 10^3)))


### Run Models

mcan.shade.ss <- MCMCglmm(seedlings ~ treat*species, 
                          random = ~treeid, 
                          family = "categorical", 
                          data = dd_shadehouse, 
                          prior = priorsshade,
                          verbose = F, 
                          pr=TRUE, #store random effects
                          pl = TRUE, #store latent variable
                          burnin = 10000,
                          nitt = 300000,
                          thin = 50) 
summary(mcan.shade.ss) #
# raftery.diag(mcan.shade.ss) #3746 samples
# autocorr(mcan.shade.ss$Sol) #
# autocorr(mcan.shade.ss$VCV) #
# allChains <- as.mcmc(cbind(mcan.shade.ss$Sol,mcan.shade.ss$VCV))  
# plot(allChains) #

# #check for overdispersion
# mcmc.units <- colMeans(mcan.shade.ss$Liab) - predict(mcan.shade.ss, marginal = NULL, type = "terms") #fitted minus observed 
# par(mfrow=c(1,2))
# qqnorm(mcmc.units) ; qqline(mcmc.units)
# hist(mcmc.units, col = grey(.6), main = "Over-dispersion")
# 
# #Check distributional assumptions of random effects
# mcmc.int <- colMeans(mcan.shade.ss$Sol[,(length(colnames(mcan.shade.ss$Sol))-2):(length(colnames(mcan.shade.ss$Sol)))]) # one value for each random effect parameter 
# qqnorm(mcmc.int) ; qqline(mcmc.int)
# hist(mcmc.int, col = grey(.6), main = "Intercept") # Not beautiful, but difficult when only 3 random effects

#Check to see if everything is significant.
cbind(B = colMeans(mcan.shade.ss$Sol[,1:(length(colnames(mcan.shade.ss$Sol))-3)]), CI = HPDinterval(mcan.shade.ss$Sol[,1:(length(colnames(mcan.shade.ss$Sol))-3)]))

# Need slightly different version of this function
MyRatios<- function(x){
  OUT <- numeric(0)
  for(i in 1:length(levels(x$species))) {
    sppnear<-x[x$species==levels(x$species)[i] & x$treat=="near","Pi"]
    sppfar<-x[x$species==levels(x$species)[i] & x$treat=="far","Pi"]
    z<- log(sppfar / sppnear)
    OUT<-rbind(OUT,z)
  }
  OUT
}



### Effect of distance on survival, Shadehouse experiments ###

#make predictions on original data (so only get predictions for actual light levels measured)
mcan.shade.sspred <- predict(mcan.shade.ss, interval="confidence", type="response", level=0.95, marginal=NULL) 
colnames(mcan.shade.sspred) <- paste("mcan.shade.ss",c("fit","lwr","upr"),sep=".")
mcan.shade.ssmod<-cbind(dd_shadehouse, mcan.shade.sspred)
ddply(mcan.shade.ssmod, .(species, treat), summarize, mean=mean(mcan.shade.ss.fit))


### create predictions "by hand"


#Predict the mean, ci's 5800 times (1 for each iteration)
#1) Get the  betas
beta.ss<-mcan.shade.ss$Sol[,1:4] #5800 values

#2) Create a grid of covariate values 
preddata.ss <- with(dd_shadehouse,expand.grid(treat=levels(treat), 
                                              species=levels(species)))

#3) For this we need the X.ss on a grid
X.ss <- model.matrix(~ treat*species, data=preddata.ss)

#4) #Calculate X.ss * beta.ss for each MCMC iteration, 4500 times, to get 4500 predictions of propsurvival for each distance and species combo!
eta.mcmc.ss<-(X.ss %*% t(beta.ss))
mu.mcmc.ss<-inv.logit(eta.mcmc.ss) #logit link
dim(mu.mcmc.ss)

#5) Calculate the mean, se, and 2.5% and 97.5% credible intervals for each of the species:dist combinations. x is a matrix of predicted values for the grid we created earlier (each species, near and far). the mean, se, and quantiles are calculated on each row. 

L.ss <- MyLinesStuff(eta.mcmc.ss)
L.ss #Posterior mean, se and 95% CI for each fixed effect parameter

#6) combine with preddata.ss, so can see which dist/spp each row is
predmcmc.ss<-cbind(preddata.ss,L.ss) #new dataframe with posterior mean, se, and CI's

#7) get mean and se on original data scale (not link)
predmcmc.ss$predpropsurv<-inv.logit(predmcmc.ss$mean)
predmcmc.ss$uci<-inv.logit(predmcmc.ss$mean + 1.96 * predmcmc.ss$se)
predmcmc.ss$lci<-inv.logit(predmcmc.ss$mean - 1.96 * predmcmc.ss$se)

#8) Graph
limits<- aes(ymin=lci, ymax=uci) #using se's
limits<-aes(ymin=inv.logit(lowcri), ymax=inv.logit(uppcri)) #using credible intervals
#get very similar results

# ggplot(predmcmc.ss, aes(x=treat, y=predpropsurv))+
#   geom_point()+
#   geom_errorbar(limits)+
#   facet_grid(.~species)

#9) Calculate differences between near and far, and create confidence intervals for these differences
preddata.ss$eta.mcmc.ss<-(X.ss %*% t(beta.ss))
preddata.ss$Pi <- exp(preddata.ss$eta.mcmc.ss) / (1 + exp(preddata.ss$eta.mcmc.ss)) #converted to probabilities (i.e. on the scale of the predictor 0 to 1) 

diffs.shade<-MyRatios(preddata.ss)

diffs.shade.sum<-MyLinesStuff(diffs.shade)
# I think species is just the species list
species <- levels(dd_shadehouse$species)

species2<-as.data.frame(species)
sppdistdiffs.shade<-cbind(species2, diffs.shade.sum)

sppdistdiffs.shade$species<-factor(sppdistdiffs.shade$species, levels=c("psycho","premna"), ordered=F)
limits<- aes(ymin=lowcri, ymax=uppcri)
# ggplot(sppdistdiffs.shade, aes(species, mean))+
#   geom_point()+
#   geom_errorbar(limits, width=0.35)+
#   theme_bw()+
#   geom_abline(intercept=0, slope=0)+ 
#   ylab("Dispersal Benefit Ratio Difference in survival near - far at mean light levels")+
#   xlab("")+
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),    text = element_text(size=20))+
#   scale_x_discrete(labels=c("psycho", "Premna"))+ geom_hline(yintercept=1)






##### Compile benefit data and integrate with network data #####
library(bipartite)

setwd(paste(topwd,"Input data/Plant benefit data",sep="/"))
SD_Saipan <- read.csv("SD_Saipan.csv",row.names=1,header=T)
fruit.traits <- read.csv("pulp_to_seed.csv",row.names=1,header=T)



sppdistdiffs.ss$life_stage <- "seed_sdlng"
sppdistdiffs$life_stage <- "sdlng"

sppdistdiffs.ss$benefit <- "ddmort"
sppdistdiffs$benefit <- "ddmort"

spplightdiffs.ss$life_stage <- "seed_sdlng"
spplightdiffs$life_stage <- "sdlng"

spplightdiffs.ss$benefit <- "highlight"
spplightdiffs$benefit <- "highlight"

sppdistdiffs.gng$life_stage <- "seed_sdlng"
sppdistdiffs.gng$benefit <- "gap"

sppdistdiffs.gut$life_stage <- "seed_sdlng_shadehouse"
sppdistdiffs.gut$benefit <- "gut"

sppdistdiffs.shade$life_stage <- "seed_sdlng_shadehouse"
sppdistdiffs.shade$benefit <- "ddmort"

# Compile posterior estimates

diffsmat <- rbind(diffs.ss,diffs,light.diffs.ss,light.diffs,diffs.gng,diffs.gut,diffs.shade)
str(diffsmat)

benefit.ratio.summary <- rbind(sppdistdiffs.ss,sppdistdiffs,spplightdiffs.ss,spplightdiffs,sppdistdiffs.gng,sppdistdiffs.gut,sppdistdiffs.shade)
benefit.ratio.summary$benefit <- as.factor(benefit.ratio.summary$benefit)
benefit.ratio.summary$benefit_stage <- as.factor(paste(benefit.ratio.summary$benefit,benefit.ratio.summary$life_stage,sep="_"))
benefit.ratio.summary$species <- revalue(benefit.ratio.summary$species,c("premna"="pr","aglaia"="ag","papaya"="pa","psychotria"="ps","psycho"="ps"))

benefit.ratio.summary$total.interactions <- rowSums(SD_Saipan)[as.character(benefit.ratio.summary$species)]
benefit.ratio.summary$species.strength <- specieslevel(SD_Saipan,level="lower",index="species strength")[as.character(benefit.ratio.summary$species),]
benefit.ratio.summary$species.degree <- specieslevel(SD_Saipan,level="lower",index="degree")[as.character(benefit.ratio.summary$species),]
benefit.ratio.summary$shannon.diversity <- specieslevel(SD_Saipan,level="lower",index="partner diversity")[as.character(benefit.ratio.summary$species),]
benefit.ratio.summary$pulp.ratios <- fruit.traits[as.character(benefit.ratio.summary$species),]

library("lme4")
benefit.ratio.mod.null <- lmer(mean ~ 1 + (1+species.strength|benefit_stage),data=benefit.ratio.summary)
benefit.ratio.mod <- lmer(mean ~ species.strength + (1+species.strength|benefit_stage),data=benefit.ratio.summary)

anova(benefit.ratio.mod.null,benefit.ratio.mod)


benefit.ratio.mod.sd.null <- lmer(mean ~ 1 + (1+species.degree|benefit_stage),data=benefit.ratio.summary)
benefit.ratio.mod.sd <- lmer(mean ~ species.degree + (1+species.degree|benefit_stage),data=benefit.ratio.summary)

anova(benefit.ratio.mod.sd.null,benefit.ratio.mod.sd)

benefit.ratio.summary$total.interactions.cent <- benefit.ratio.summary$total.interactions-mean(benefit.ratio.summary$total.interactions)
benefit.ratio.mod.ti.null <- lmer(mean ~ 1 + (total.interactions.cent|benefit_stage),data=benefit.ratio.summary)
benefit.ratio.mod.ti <- lmer(mean ~ total.interactions + (total.interactions.cent|benefit_stage),data=benefit.ratio.summary)

anova(benefit.ratio.mod.ti.null,benefit.ratio.mod.ti)


benefit.ratio.mod.shan.null <- lmer(mean ~ 1 + (1+shannon.diversity|benefit_stage),data=benefit.ratio.summary)
benefit.ratio.mod.shan <- lmer(mean ~ shannon.diversity + (1+shannon.diversity|benefit_stage),data=benefit.ratio.summary)

anova(benefit.ratio.mod.shan.null,benefit.ratio.mod.shan)



estimates.benefit.ratio <- coef(benefit.ratio.mod)$benefit_stage

##### Figures

make.pdf <- F # Specify whether you want pdfs of these individually

setwd(paste(topwd,"Figures",sep="/"))

colT <- rgb(0,0,0,max=255,alpha=2)
col1 <- rgb(217,95,2,max=255,alpha=255)
col2 <- rgb(27,158,119,max=255,alpha=255)
col3 <- rgb(117,112,179,max=255,alpha=255)


if(make.pdf == T) {
  par(mfcol=c(1,1),pin=c(2,1.66))
} else {
  par(mfcol=c(4,2),pin=c(2,1.66))
}


bs <- levels(factor(benefit.ratio.summary$benefit_stage))
for(i in 1:length(bs)) {
  if(make.pdf == T) {
    pdf(paste("benefit",i,".pdf",sep=""),width=3,height=3,useDingbats=F)
  }
  
  ben.set <- subset(benefit.ratio.summary,benefit.ratio.summary$benefit_stage==bs[i])
  head(ben.set)
  min.bs <- floor(min(ben.set$mean))
  max.bs <- ceiling(max(ben.set$mean))
  
  min.bs <- floor(min(ben.set$lowcri))
  max.bs <- ceiling(max(ben.set$uppcri))
  
  min.bs <- ifelse(min.bs>-1,-1,min.bs)
  max.bs <- ifelse(max.bs<0,0,max.bs)
  if(make.pdf == T) {
  } else {
    if(i==5) plot.new()
  }
  
  plot(mean~species.strength,
       data=ben.set,
       pch=(14+c(1,2,2,2,2,1,2,1))[i],
       col=1,#"white",
       cex=1.2,
       frame.plot=F,
       ylim=c(min.bs,max.bs),
       xlim=c(0,3),#2.6
       ylab="",#"Dispersal Benefit Ratio",
       xlab="",
       xaxt="n",
       yaxt="n")
  if(make.pdf == T) {
    axis(1,at=c(0,1,2,3),labels=c(0,1,2,3))
  } else {
    if(i==4) axis(1,at=c(0,1,2,3),labels=c(0,1,2,3))
    if(i==7) axis(1,at=c(0,1,2,3),labels=c(0,1,2,3))
  }
  
  axis(2,at=c(min.bs,0,max.bs),labels=parse(text = sprintf("e^%d", c(min.bs,0,max.bs))),
       las=1)
  abline(h=0,lty=2)
  
  bs.inds <- which(benefit.ratio.summary$benefit_stage==bs[i])
  
  diffvec <- c()
  for(j in bs.inds){
    diffvec <- c(diffvec,diffsmat[j,])
  }
  
  points(x=jitter(rep(ben.set$species.strength,each=dim(diffsmat)[2]),amount=0.05),y=diffvec,pch=".",col=colT)
  
  legend("topleft",c("Distance, seedling",
                     "Distance, seeds",
                     "Distance, shadehouse seeds",
                     "Gap, seeds",
                     "Handling, shadehouse seeds",
                     "High light, seedlngs",
                     "High light, seeds")[i],
         bty="n",
         cex=.7)
  
  curve(estimates.benefit.ratio[i,1]+estimates.benefit.ratio[i,2]*x,
        add=T)#,col=c(col1,col1,col2,col3,col3)[i])
  if(make.pdf == T) dev.off()
  
}






###
### Figures
###

# Network figure
par(pin=c(2.5,2.5),mfrow=c(1,1))
if(make.pdf==T) pdf("SD_Saipan_visweb.pdf",width=3,height=3,useDingbats=F)
visweb(SD_Saipan)
if(make.pdf==T) dev.off()

# Investment and Dependence figures


sp.str.max <- max(specieslevel(SD_Saipan,level="lower",index="species strength"))

par(pin=c(1,1.5),mfrow=c(1,1))
if(make.pdf==T) pdf("Investment vs species strength and mean dependence.pdf",width=5,height=3,useDingbats=F)
par(mfrow=c(1,2))

par(xpd=T)
plot(unlist(specieslevel(SD_Saipan,index="species strength",level="lower")) ~ fruit.traits[rownames(SD_Saipan),],
     ylim=c(0,sp.str.max),
     frame.plot=F,
     pch=21,
     col="black",
     xlab="Investment (pulp:seed ratio)",
     ylab="Species strength",
     xlim=c(0,2),
     las=1,
     xaxt="n",
     yaxt="n")
axis(2,at=c(0,1.25,2.5),labels=c("0","1.25","2.5"),las=1)
axis(1,at=c(0,1,2),labels=c(0,1,2),las=1)

pulp.seed.mod <- lm(unlist(specieslevel(SD_Saipan,index="species strength",level="lower")) ~ fruit.traits[rownames(SD_Saipan),])
par(xpd=F)
curve(coef(pulp.seed.mod)[1] + coef(pulp.seed.mod)[2]*x,add=T)
summary(pulp.seed.mod)




#######



summary(benefit.ratio.mod)

newdat <- expand.grid(
  species.strength=seq(0,sp.str.max,length.out=100)
  , mean = 0
)
mm <- model.matrix(terms(benefit.ratio.mod),newdat)
newdat$mean <- predict(benefit.ratio.mod,newdat,re.form=NA)
## or newdat$mean <- mm %*% fixef(benefit.ratio.mod)
pvar1 <- diag(mm %*% tcrossprod(vcov(benefit.ratio.mod),mm))
cmult <- 1.96
newdat <- data.frame(
  newdat
  , plo = newdat$mean-cmult*sqrt(pvar1)
  , phi = newdat$mean+cmult*sqrt(pvar1)
)
#plot confidence
#par(pin=c(3,3))
plot(NULL,xlim=c(0,sp.str.max),ylim=c(log(0.5),log(25)), 
     frame.plot=F, 
     las=1,
     ylab="Dispersal benefit ratio",
     xlab="Species strength",
     xaxt="n",
     yaxt="n")
axis(1,at=c(0,1.25,2.5),labels=c("0","1.25","2.5"),las=1)
y.ax.lab <- c(0.5,1,5,25)
axis(2, at=log(y.ax.lab),labels=y.ax.lab, las=2)
polygon(x=c(rev(newdat$species.strength),newdat$species.strength),
        y=c(rev(newdat$plo),newdat$phi),col="grey",
        border=NA)
curve(fixef(benefit.ratio.mod)[1] + fixef(benefit.ratio.mod)[2] * x, add=T,
      lwd=2,lend="butt",xlim=c(0,max(benefit.ratio.summary$species.strength)))
curve(x*0,xlim=c(0,sp.str.max),add=T,lty=2)

if(make.pdf==T) dev.off()
