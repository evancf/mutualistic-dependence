# This code analyzes the diet of frugivores in seed dispersal networks, and
# performs a number of network simulations.
# The simulation code is modified the coextNumber function
# from Vieira and Almeida-Neto (2015) in Ecology Letters. These
# simulations therefore were built with for-loops, and are very slow.
# The simulations themselves are commented out and simulation output is
# read from file instead. See Network Simulations and Diet Figures
# code to produce figures.

topwd <- getwd()

##### Read in files, clean data

#load libraries & source files
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("bipartite","lme4","RColorBrewer")

ipak(packages)

source("Network Simulation Functions.R")

##### Load bird diet info

bird.traits <- read.csv("Input data/BirdFuncDat.csv",header=T)
colnames(bird.traits)[4] <- "Family"
mam.traits <- read.csv("Input data/MamFuncDat.csv",header=T)
colnames(mam.traits)[3] <- "Family"

# Merge birds and mammals into the same dataset (only need this if there are mammals)
# Mammals included here if future additional networks contain mammals
mam.merge <- mam.traits[,c(2,3:15,24)]
mam.merge$taxon <- "mammal"
bird.merge <- bird.traits[,c(8,4,10:19,21:22,36)]
bird.merge$taxon <- "bird"
traits <- rbind(mam.merge,bird.merge)
traits$genus <- gsub(' [A-z ]*', '', traits$Scientific)

##### Make database with network metrics

filevec.sd <- list.files(path="Input data",pattern="M_SD")

frug.sp <- c()
plant.sp <- c()
frug.species.strength <- c()
frug.tot.interactions <- c()
frug.shannon.diversity <- c()
frug.degree <- c()
frug.species.strength.rel1 <- c()
frug.degree.rel1 <- c()
frug.species.strength.lr1 <- c()
frug.degree.lr1 <- c()
net.id <- c()
binary <- c()
for(i in 1:length(filevec.sd)){
  
  imatrix.df <- sortweb(read.csv(paste("Input data",filevec.sd[i],sep="/"),header=T,row.names=1))
  frug.sp <- c(frug.sp,colnames(imatrix.df))
  plant.sp <- c(plant.sp,rownames(imatrix.df))
  frug.species.strength <- c(frug.species.strength,as.numeric(unlist(specieslevel(imatrix.df,index="species strength",level="higher"))))
  frug.tot.interactions <- c(frug.tot.interactions,as.numeric(colSums(imatrix.df)))
  frug.shannon.diversity <- c(frug.shannon.diversity,as.numeric(unlist(specieslevel(imatrix.df,index="partner diversity",level="higher"))))
  frug.degree <- c(frug.degree,as.numeric(unlist(specieslevel(imatrix.df,index="degree",level="higher"))))
  frug.species.strength.rel1 <- c(frug.species.strength.rel1,rel1(as.numeric(unlist(specieslevel(imatrix.df,index="species strength",level="higher")))))
  frug.degree.rel1 <- c(frug.degree.rel1,rel1(as.numeric(unlist(specieslevel(imatrix.df,index="degree",level="higher")))))
  frug.species.strength.lr1 <- c(frug.species.strength.lr1,logrel1(as.numeric(unlist(specieslevel(imatrix.df,index="species strength",level="higher")))))
  frug.degree.lr1 <- c(frug.degree.lr1,logrel1(as.numeric(unlist(specieslevel(imatrix.df,index="degree",level="higher")))))
  net.id <- c(net.id,rep(filevec.sd[i],dim(imatrix.df)[2]))
  binary <- c(binary,rep(ifelse(max(imatrix.df)==1,"binary","quant"),dim(imatrix.df)[2]))
}
frug.sp <- gsub("[[:punct:]]"," ",frug.sp)
frug.sp <- gsub("[[:space:]]*$","",frug.sp) #get rid of trailing spaces
plant.sp <- gsub("[[:space:]]*$","",plant.sp) #get rid of trailing spaces

#
net.traits <- as.data.frame(cbind(frug.sp,frug.species.strength,frug.tot.interactions,frug.shannon.diversity,
                                  frug.degree,frug.species.strength.rel1,
                                  frug.degree.rel1,frug.species.strength.lr1,
                                  frug.degree.lr1,net.id,binary))
colnames(net.traits) <- c("sp","species.strength","total.interactions","shannon.diversity","degree",
                          "species.strength.rel1","degree.rel1",
                          "species.strength.lr1","degree.lr1","net.id","binary")
net.traits.genus <- gsub(' [A-z ]*', '', net.traits$sp)
# Fix up some naming issues

naming.issues <- read.csv("Input data/naming.issues.csv",header=T)

for(i in 1:dim(naming.issues)[1]){
  bad.name <- naming.issues[i,2]
  good.name <- naming.issues[i,3]
  net.traits$sp <- gsub(bad.name,good.name,net.traits$sp)
}

bad.indices <- which(net.traits$sp %!in% traits$Scientific)
bad.sp <- net.traits$sp[bad.indices]# This tells us which from the networks aren't in the trait database
unknowns <- bad.sp[order(bad.sp)]



traits[which(substr(traits$Scientific,1,6)=="Columb"),"Scientific"]

# Ignore unidenfified ones or those identified to genus
#net.traits$sp[which(net.traits$sp %in% unknowns)] <- NA


# Network info

# Add in data to a dataframe including network description
net.descr <- read.csv("Output data/net.descr.csv",row.names=1,header=T)
net.descr$number_species <- rowSums(net.descr[,2:3]) # get total number of spp
colnames(net.descr)[1] <- "net_id_csv"
net.descr$net_id <- gsub(".csv","",net.descr$net_id_csv)
references <- read.csv("Input data/references.csv",row.names=1,header=T)
net.descr[,(dim(net.descr)[2]+1):(dim(net.descr)[2]+dim(references)[2]+1)] <- references[net.descr$net_id,]
rownames(net.descr) <- net.descr$net_id_csv

##### Add diet data

net.traits[,colnames(traits)] <- NA

for(i in 1:length(net.traits$sp)){
  if(i %in% bad.indices){
    # Could use this if you wanted to give species identified to genus the median values among congeners
    # a<-NA
    # b<-traits[which(net.traits.genus[i] == traits$genus)[1],2]
    # c<-apply(traits[which(net.traits.genus[i] == traits$genus),3:12],2,median)
    # d<-rep(NA,5)
    # 
    # net.traits[i,12:(11+length(colnames(traits)))]<-c(a,b,c,d)
    net.traits[i,12:(11+length(colnames(traits)))]<-rep(NA,dim(traits)[2])
    
  }else{
    net.traits[i,12:(11+length(colnames(traits)))]<-traits[which(net.traits$sp[i] == traits$Scientific),]
  }
  
}

colnames(net.traits)
str(net.traits)
# This is odd that it isn't giving numeric # Make this more efficient
net.traits$species.strength <- as.numeric(as.character(net.traits$species.strength))
net.traits$total.interactions <- as.numeric(as.character(net.traits$total.interactions))
net.traits$shannon.diversity <- as.numeric(as.character(net.traits$shannon.diversity))
net.traits$degree <- as.numeric(as.character(net.traits$degree))
net.traits$species.strength.rel1 <- as.numeric(as.character(net.traits$species.strength.rel1))
net.traits$degree.rel1 <- as.numeric(as.character(net.traits$degree.rel1))
net.traits$species.strength.lr1 <- as.numeric(as.character(net.traits$species.strength.lr1))
net.traits$degree.lr1 <- as.numeric(as.character(net.traits$degree.lr1))
net.traits$sp <- as.factor(net.traits$sp)
net.traits$binary <- as.factor(net.traits$binary)



##### Analysis of diet information

# Examine relationship between species strength and portion of fruit in diet

# Get clean data set, get columns to be used for logistic regression
net.traitsnoNA <- net.traits[complete.cases(net.traits),]
net.traitsnoNA$log.total.interactions <- log(net.traitsnoNA$total.interactions)
net.traitsnoNA$log.degree <- log(net.traitsnoNA$degree)
net.traitsnoNA$log.shannon.diversity <- log(net.traitsnoNA$shannon.diversity)
net.traitsnoNA$log.species.strength <- log(net.traitsnoNA$species.strength)
net.traitsnoNA$fruit.yes <- net.traitsnoNA$Diet.Fruit/10
net.traitsnoNA$fruit.no <- 10-net.traitsnoNA$Diet.Fruit/10

### Quantitative Network Analyses
### Run Models Species Strength
diet.mod <- glmer(cbind(fruit.yes,fruit.no) ~ log.species.strength + (log.species.strength|net.id), 
                  family = "binomial",
                  data = subset(net.traitsnoNA,binary=="quant"))
diet.mod.null <- glmer(cbind(fruit.yes,fruit.no) ~ 1 + (log.species.strength|net.id), 
                       family = "binomial",
                       data = subset(net.traitsnoNA,binary=="quant"))

# Likelihood ratio test - assess whether there is a positive relationship
anova(diet.mod,diet.mod.null)



### Run Models Shannon Diversity
diet.mod.shan <- glmer(cbind(fruit.yes,fruit.no) ~ shannon.diversity + (shannon.diversity|net.id), 
                       family = "binomial",
                       data = subset(net.traitsnoNA,binary=="quant"))
diet.mod.shan.null <- glmer(cbind(fruit.yes,fruit.no) ~ 1 + (shannon.diversity|net.id), 
                            family = "binomial",
                            data = subset(net.traitsnoNA,binary=="quant"))

# Likelihood ratio test - assess whether there is a positive relationship
anova(diet.mod.shan,diet.mod.shan.null)


### Run Models Species Degree
diet.mod.degree <- glmer(cbind(fruit.yes,fruit.no) ~ log.degree + (log.degree|net.id), 
                         family = "binomial",
                         data = subset(net.traitsnoNA,binary=="quant"))
diet.mod.degree.null <- glmer(cbind(fruit.yes,fruit.no) ~ 1 + (log.degree|net.id), 
                              family = "binomial",
                              data = subset(net.traitsnoNA,binary=="quant"))

# Likelihood ratio test - assess whether there is a positive relationship
anova(diet.mod.degree,diet.mod.degree.null)


### Run Models Total Interactions
diet.mod.totint <- glmer(cbind(fruit.yes,fruit.no) ~ log.total.interactions + (log.total.interactions|net.id), 
                         family = "binomial",
                         data = subset(net.traitsnoNA,binary=="quant"))
diet.mod.totint.null <- glmer(cbind(fruit.yes,fruit.no) ~ 1 + (log.total.interactions|net.id), 
                              family = "binomial",
                              data = subset(net.traitsnoNA,binary=="quant"))

# Likelihood ratio test - assess whether there is a positive relationship
anova(diet.mod.totint,diet.mod.totint.null)

summary(diet.mod) #


# Get estimates for later use 
diet.estimates.mean <- summary(diet.mod)$coefficients
diet.estimates <- coef(diet.mod)$net.id




### Binary Network Analyses



### Run Models Species Degree
diet.mod.binary.degree <- glmer(cbind(fruit.yes,fruit.no) ~ log.degree + (log.degree|net.id), 
                         family = "binomial",
                         data = subset(net.traitsnoNA,binary=="binary"))
diet.mod.binary.degree.null <- glmer(cbind(fruit.yes,fruit.no) ~ 1 + (log.degree|net.id), 
                              family = "binomial",
                              data = subset(net.traitsnoNA,binary=="binary"))

# Likelihood ratio test - assess whether there is a positive relationship
anova(diet.mod.binary.degree,diet.mod.binary.degree.null)

summary(diet.mod.binary.degree)


diet.estimates.binary.mean <- summary(diet.mod.binary.degree)$coefficients
diet.estimates.binary <- coef(diet.mod.binary.degree)$net.id

coef(diet.mod.binary.degree)
coef(diet.mod.degree)



















# #### Simulations
# 
# # Simulations using empirical Seed dispersal networks
# # and diet info to predict extinctions in the presence and absence
# # of empirical structure and of the observed species strength - dependence relationship
# 
# # "slope" is used to represent the observed relationship between species strength and dependence
# # "null" indicate scenarios where real networks are randomized
# # "obligate" refers to the assumption that all species are obligate mutualists,
# # regardless of species strength
# 
# 
# # Adapted from coextNumber function
# 
# ### Specify parameters for the for-loop
# nsims <- 10000 # This is the number of simulations to test for
# # each network for each scenario
# 
# # Specify how many to sample for detailed info on the species strength of
# # all species. Otherwise simulations get super slow and models take
# # forever to run.
# thinning <- seq(0,nsims,5)
# 
# ### Get set up for for-loop simulation
# 
# filevec.sd <- c("M_SD_001.csv",
#                 "M_SD_002.csv",
#                 "M_SD_003.csv",
#                 "M_SD_004.csv",
#                 "M_SD_005.csv",
#                 "M_SD_006.csv",
#                 "M_SD_009.csv",
#                 "M_SD_010.csv",
#                 "M_SD_012.csv",
#                 "M_SD_020.csv",
#                 "M_SD_023.csv")
# 
# netvec.sd <- gsub(".csv","",filevec.sd)
# 
# 
# # Set up output
# net.diff <- as.data.frame(matrix(NA,nrow=(nsims*length(filevec.sd)),ncol=6))
# colnames(net.diff) <- c("net_id",
#                         "ext_counts_slope",
#                         "ext_counts_null_slope",
#                         "ext_counts_obligate",
#                         "ext_counts_null_obligate",
#                         "sim_id")
# 
# 
# ## Set up output for getting the species strengths of each coextinct species
# 
# status.sp.str.obligate <- data.frame("species.strength"=NA,"coextinct"=NA, "network"=NA)
# status.sp.str.observed <- data.frame("species.strength"=NA,"coextinct"=NA, "network"=NA)
# 
# # Function to extract these species strength values
# coext.sp.strs <- function(x) {
# 
#   animals <- unlist(ifelse(class(x[[2]])=="data.frame",
#                            subset(x[[2]],degree_of_extinction>1),NA))
#   plants <- unlist(ifelse(class(x[[3]])=="data.frame",
#                           subset(x[[3]],degree_of_extinction>1),NA))
#   coext.vec <- c(animals,nanimal+plants)
# 
#   mat <- as.data.frame(matrix(c(anim.spst,plant.spst,
#                                 rep(0,nmut),
#                                 rep(netvec.sd[j],nmut)),
#                               nrow=nmut,
#                               ncol=3))
# 
#   mat[,2] <- factor(mat[,2],levels=c(0,1))
# 
#   colnames(mat) <- c("species.strength","coextinct","network")
#   if(is.na(coext.vec[1])){
#     mat
#   } else {
#     mat[coext.vec,2] <- 1
#     mat
#   }
# 
# }
# 
# ### Run simulation
# 
# time1 <- Sys.time()
# set.seed(11)
# for(j in 1:length(filevec.sd)) {
#   imatrix.df <- sortweb(read.csv(paste("Input data",filevec.sd[j],sep="/"),header=T,row.names=1))
#   imatrix <- sortweb(as.matrix(imatrix.df,rownames.force=F))
# 
# 
#   # Get log species strength for each animal and each plant species
#   anim.spst <- log(unlist(specieslevel(imatrix.df,index="species strength",level="lower")))
#   plant.spst <- log(unlist(specieslevel(imatrix.df,index="species strength",level="higher")))
# 
#   # Assign Ri for each species based on their species strength and the
#   # the bird diet model fit (using fixed.estimates.logit function)
#   ranim_slope <- fixed.pred.logit(diet.mod,anim.spst)
#   rplants_slope <- fixed.pred.logit(diet.mod,plant.spst)
# 
#   # Make null matrix for comparison
#   # These provide randomized null models (using null model r2dtable)
#   # that have the same row and column sums
#   nmatrix <- nullmodel(imatrix.df,N=1,method=1)[[1]]
# 
#   # We will use the same assigned Ri in the structure "null" (randomized) scenarios
#   ranim_slope_null <- ranim_slope
#   rplants_slope_null <- rplants_slope
# 
#   nanimal <- nrow(imatrix)
#   nplant <- ncol(imatrix)
#   nmut <- nanimal + nplant
# 
#   for(sim in 1:nsims){
# 
#     sim.ind <- sim+(j-1)*nsims # This gives the index for saving data
# 
#     #Choose target guild and target species to start potential cascade
#     #Use the same target species for all scenarios for an individual iteration of "sim"
#     guild <- sample(c('animal','plant'),1,F,c(nanimal,nplant))
#     if(guild=="animal"){
#       target <- sample(1:nanimal,1)#rep(1, nanimal)
#       # If you want to chose a different one for each of the scenarios, use
#       # the commented out bit of code here to give each species an equal chance
# 
#     }else{
#       target <- sample(1:nplant,1)#rep(1, nplant)#
#     }
#     nc.ecs <- netcascade(imatrix,
#                          ranim = ranim_slope,
#                          rplants = rplants_slope,
#                          targetGuild = guild,
#                          target = target)
#     net.diff$ext_counts_slope[sim.ind] <- sum(nc.ecs[[1]]$n_extinctions)
#     if(sim %in% thinning){
#       status.sp.str.observed <- rbind(status.sp.str.observed,coext.sp.strs(nc.ecs),make.row.names=F)
#     }
# 
# 
#     # Now null models that allows comparisons to randomized architecture
#     nc.ecns <- netcascade(nmatrix,
#                           ranim = ranim_slope,
#                           rplants = rplants_slope,
#                           targetGuild = guild,
#                           target = target)
#     net.diff$ext_counts_null_slope[sim.ind] <- sum(nc.ecns[[1]]$n_extinctions)
# 
#     # Now for the "obligate" case - topological approach. Use nested and non-nested networks
#     nc.eco <- netcascade(imatrix,
#                          ranim = rep(1,nanimal),
#                          rplants = rep(1,nplant),
#                          targetGuild = guild,
#                          target = target)
#     net.diff$ext_counts_obligate[sim.ind] <- sum(nc.eco[[1]]$n_extinctions)
#     if(sim %in% thinning){
#       status.sp.str.obligate <- rbind(status.sp.str.obligate,coext.sp.strs(nc.eco),make.row.names=F)
#     }
# 
#     nc.ecno <- netcascade(nmatrix,
#                           ranim = rep(1,nanimal),
#                           rplants = rep(1,nplant),
#                           targetGuild = guild,
#                           target = target)
#     net.diff$ext_counts_null_obligate[sim.ind] <- sum(nc.ecno[[1]]$n_extinctions)
# 
#     net.diff$net_id[sim.ind] <- filevec.sd[j]
#     net.diff$sim_id[sim.ind] <- paste(filevec.sd[j],sim,guild,target,sep="_")
# 
#   }
# 
#   print(paste("j =",j))
# }
# 
# time2 <- Sys.time()
# time2-time1
# # Save this to "Output data" folder
# setwd(paste(topwd,"Output data",sep="/"))
# write.csv(net.diff,"net.diff.csv")
# write.csv(status.sp.str.obligate,"status.sp.str.obligate.csv")
# write.csv(status.sp.str.observed,"status.sp.str.observed.csv")
# setwd(topwd)
# 































# If not performing simulation, read from "Output data" folder
setwd(paste(topwd,"Output data",sep="/"))
net.diff <- read.csv("net.diff.csv",row.names=1,header=T)
status.sp.str.observed <- read.csv("status.sp.str.observed.csv",row.names=1,header=T)
status.sp.str.obligate <- read.csv("status.sp.str.obligate.csv",row.names=1,header=T)

setwd(topwd)

# Get this into a useable format
net.diff.long <- reshape(net.diff,direction="long",
                         varying=list(2:5),
                         v.names="extinctions",
                         times=c("ext_counts_slope",
                                 "ext_counts_null_slope",
                                 "ext_counts_obligate",
                                 "ext_counts_null_obligate"))

net.diff.long$time <- as.factor(net.diff.long$time)
net.diff.long$time <- factor(net.diff.long$time, levels=c("ext_counts_null_obligate",
                                                          "ext_counts_obligate",
                                                          "ext_counts_null_slope",
                                                          "ext_counts_slope"))
net.diff.long$number.species <- net.descr[as.character(net.diff.long$net_id),"number_species"]

# p.ext refers to probability of coextinctions
p.ext <- tapply((net.diff.long$extinctions-1)/net.diff.long$number.species,
                list(net.diff.long$time,net.diff.long$net_id),
                mean)

# n.ext refers to number of coextinctions
n.ext <- tapply(net.diff.long$extinctions-1,
                list(net.diff.long$time,net.diff.long$net_id),
                mean)


# Determine percent decrease

1-mean(p.ext[3,]/p.ext[1,]) # This is overall decrease in coextinciton - obligate scenario
1-mean(p.ext[4,]/p.ext[2,]) # This is overall decrease in coextinciton - observed scenario
1-mean(p.ext[2,]/p.ext[1,]) # Decrease in coextinction conferred by non-random structure - obligate scenario
1-mean(p.ext[4,]/p.ext[3,]) # Decrease in coextinction conferred by non-random structure - observed scenario


# Analyze the relationship between species strength and coextinction probability 
# in the "observed" and "obligate" scenarios
status.sp.str.observed$species.strength <- as.numeric(status.sp.str.observed$species.strength)
status.sp.str.obligate$species.strength <- as.numeric(status.sp.str.obligate$species.strength)
status.sp.str.observed$coextinct <- as.numeric(status.sp.str.observed$coextinct)
status.sp.str.obligate$coextinct <- as.numeric(status.sp.str.obligate$coextinct)

observed.mod <- glmer(coextinct ~ species.strength + (species.strength|network), family="binomial", data=status.sp.str.observed)
obligate.mod <- glmer(coextinct ~ species.strength + (species.strength|network), family="binomial", data=status.sp.str.obligate)





##### Coextinction Simulations using 11 empirical networks

# These simulations are used to examine how non-random structure and the slope 
# each influence coextinction predictions. Instead of 
# defining the strength of the slope using empirical data on bird diets
# as above, these simulations predict coextinctions at varying strengths
# of the slope (relationship between species strength and mutualistic dependence).

###
### Adapted from coextNumber function
###

### Specify parameters for the for-loop
seq.length <- 100 # This is the number of different slopes to test
nsims <- 200 # This is the number of simulations to test for each slope
R <- c(1)#,0.75,0.5,0.25) # These are the maximum Ri to test.

### Get set up for for-loop
# Get sequence of numbers between 0 and 1 that will serve as the
# y intercept of the intrinsic dependence - species strength relationship
seq10 <- seq(from=1, to=0, length=seq.length)
filevec <- list.files(path="Input data",pattern="M_")

### Set up output

# Want columns to be network id, tradeoff strength, slope extincts, horiz extincts, slope null extincts, horiz null extincts
# "horiz" is equivalent to the "slope null" scenarios
# "tradoeff_strength" is equivalent to "slope value"
# "null" indicate scenarios where real networks are randomized

net.out <- as.data.frame(matrix(NA,nrow=(seq.length*length(filevec)),ncol=8))
colnames(net.out) <- c("net_id","tradeoff_strength",
                       "slope_ext","horiz_ext",
                       "slope_null_ext","horiz_null_ext",
                       "obligate_ext","obligate_null_ext")

net.list <- rep(list(net.out),length(R))

# time1 <- Sys.time()
# set.seed(101)
# 
# for(k in 1:length(R)){
#   for(j in 1:length(filevec)) {
#     imatrix.df <- sortweb(read.csv(paste("Input data",filevec[j],sep="/"),header=T,row.names=1))
#     imatrix <- sortweb(as.matrix(imatrix.df,rownames.force=F))
#     
#     # Temporary data frames for output of netcascade
#     
#     ext_counts_slope <- as.data.frame(c())
#     ext_counts_horiz <- as.data.frame(c())
#     
#     ext_counts_null_slope <- as.data.frame(c())
#     ext_counts_null_horiz <- as.data.frame(c())
#     
#     ext_counts_obligate <- as.data.frame(c())
#     ext_counts_null_obligate <- as.data.frame(c())
#     
#     # Get relative species degree for each animal and each plant species
#     anim.spst <- rel1(unlist(specieslevel(imatrix.df,index="species strength",level="lower")))
#     plant.spst <- rel1(unlist(specieslevel(imatrix.df,index="species strength",level="higher")))
#     
#     for(i in 1:seq.length) {
#       for(sim in 1:nsims){
#         
#         # Assign Ri for each species based on their species strength and the
#         # slope that is used for this i in seq.length
#         ranim_slope <- (R[k]*seq10[i]) + anim.spst*(R[k]-(R[k]*seq10[i]))
#         rplants_slope <- (R[k]*seq10[i]) + plant.spst*(R[k]-(R[k]*seq10[i]))
#         
#         # Assign Ri based on the null (horizontal) relationship between Ri and 
#         # species strength that integrates to the same value for the current 
#         # value of slope
#         ranim_horiz <- rep(((R[k]*seq10[i])+R[k])/2,nrow(imatrix.df))
#         rplants_horiz <- rep(((R[k]*seq10[i])+R[k])/2,ncol(imatrix.df))
#         
#         # Make null matrix for comparison
#         # These provide randomized null models (using null model r2dtable)
#         # that have the same row and column sums
#         nmatrix <- nullmodel(imatrix.df,N=1,method=1)[[1]]
#         
#         # We will use the same assigned Ri in the structure "null" (randomized) scenarios
#         
#         ranim_slope_null <- ranim_slope
#         rplants_slope_null <- rplants_slope
#         ranim_horiz_null <- ranim_horiz
#         rplants_horiz_null <- rplants_horiz
#         
#         #Choose target guild and target species to start potential cascade
#         #Use the same target species for all scenarios for an individual iteration of "sim"
#         guild <- sample(c('animal','plant'),1,F,c(nrow(imatrix),ncol(imatrix)))
#         if(guild=="animal"){
#           target <- sample(1:nrow(imatrix),1)#rep(1, nrow(imatrix))
#           # If you want to chose a different one for each of the scenarios, use 
#           # the commented out bit of code here to give each species an equal chance
#         }else{
#           target <- sample(1:ncol(imatrix),1)#rep(1, ncol(imatrix))#
#         }    
#         ext_counts_slope[sim,i] <- sum(netcascade(imatrix,
#                                                   ranim = ranim_slope, 
#                                                   rplants = rplants_slope, 
#                                                   targetGuild = guild, 
#                                                   target = target)[[1]]$n_extinctions)
#         
#         # Now null model that allows comparison for equal integral between species strength and intrinsic dependence
#         ext_counts_horiz[sim,i] <- sum(netcascade(imatrix,
#                                                   ranim = ranim_horiz, 
#                                                   rplants = rplants_horiz, 
#                                                   targetGuild = guild, 
#                                                   target = target)[[1]]$n_extinctions)
#         
#         # Now null models that allows comparisons to randomized architecture
#         ext_counts_null_slope[sim,i] <- sum(netcascade(nmatrix,
#                                                        ranim = ranim_slope,
#                                                        rplants = rplants_slope,
#                                                        targetGuild = guild,
#                                                        target = target)[[1]]$n_extinctions)
#         
#         ext_counts_null_horiz[sim,i] <- sum(netcascade(nmatrix,
#                                                        ranim = ranim_horiz, 
#                                                        rplants = rplants_horiz, 
#                                                        targetGuild = guild, 
#                                                        target = target)[[1]]$n_extinctions)
#         
#         # Now for the "obligate" case - topological approach. Use empirical and randomized networks
#         ext_counts_obligate[sim,i] <- sum(netcascade(imatrix,
#                                                      ranim = rep(1,nrow(imatrix)), 
#                                                      rplants = rep(1,ncol(imatrix)), 
#                                                      targetGuild = guild, 
#                                                      target = target)[[1]]$n_extinctions)
# 
#         ext_counts_null_obligate[sim,i] <- sum(netcascade(nmatrix,
#                                                           ranim = rep(1,nrow(imatrix)), 
#                                                           rplants = rep(1,ncol(imatrix)), 
#                                                           targetGuild = guild, 
#                                                           target = target)[[1]]$n_extinctions)
#         
#         
#       }
#       
#       
#       print(paste("j =",j,"and i =",i))
#     }
#     
#     
#     # Output ("net_id","tradeoff_strength","slope_ext","horiz_ext","slope_null_ext","horiz_null_ext")
#     
#     indices <- ((j-1)*seq.length+1):((j-1)*seq.length+seq.length)
#     net.list[[k]]$net_id[indices] <- rep(filevec[j],seq.length)
#     net.list[[k]]$tradeoff_strength[indices] <- rev(seq10)
#     
#     net.list[[k]]$slope_ext[indices] <- colSums(ext_counts_slope)
#     net.list[[k]]$horiz_ext[indices] <- colSums(ext_counts_horiz)
#     net.list[[k]]$slope_null_ext[indices] <- colSums(ext_counts_null_slope)
#     net.list[[k]]$horiz_null_ext[indices] <- colSums(ext_counts_null_horiz)
#     net.list[[k]]$obligate_ext[indices] <- colSums(ext_counts_obligate)
#     net.list[[k]]$obligate_null_ext[indices] <- colSums(ext_counts_null_obligate)
#     
#     ## These allow us to get individual comparisons between 
#     #net.list[[k]]$slope_v_obligate[indices] <- colSums(ext_counts_slope/ext_counts_obligate)
#     #net.list[[k]]$slope_v_horiz[indices] <- colSums(ext_counts_slope/ext_counts_horiz)
# 
#     
#   }
#   
#   
# }
# 
# time2 <- Sys.time()
# time2-time1
# 
# ###
# ### Second for-loop for network level parameters
# ###
# # Don't need this if reading it in above!
# 
# # Want other data frame just listing some descriptors of 
# net.descr <- data.frame(character(),numeric(),numeric(),numeric(),numeric(),stringsAsFactors=F)
# colnames(net.descr) <- c("net_id","anim_sp","plant_sp","connectance","weighted_nestedness")
# 
# for(j in 1:length(filevec)) {
#   imatrix.df <- sortweb(read.csv(paste("Input data",filevec[j],sep="/"),header=T,row.names=1))
#   
#   net.descr[j,] <- c(filevec[j],nrow(imatrix.df),ncol(imatrix.df),
#                      networklevel(imatrix.df,"connectance"),
#                      networklevel(imatrix.df,"weighted nestedness"))
# }
# setwd(paste(topwd,"Output data",sep="/"))
# write.csv(net.descr,"net.descr.csv")
# setwd(topwd)
# 
# ##
# ## Write output files to CSV
# ##
# 
# setwd(paste(topwd,"Output data",sep="/"))
# for(i in 1:length(R)){
#   write.csv(net.list[[i]], paste("net.out",R[i],"csv",sep="."))
# }
# setwd(topwd)


### Read in files
setwd(paste(topwd,"Output data",sep="/"))
for(i in 1:length(R)){
  net.list[[i]] <- read.csv(paste("net.out",R[i],"csv",sep="."),row.names=1,header=T)
}
setwd(topwd)
net.out <- net.list[[1]] # Choose the max R value to look at

net.out$slope_ext <- net.out$slope_ext-(nsims) #This -(nsims) is because we want to measure coextinctions. There are nsims extinctions (to seed the simulation)
net.out$horiz_ext <- net.out$horiz_ext-(nsims)
net.out$slope_null_ext <- net.out$slope_null_ext-(nsims)
net.out$horiz_null_ext <- net.out$horiz_null_ext-(nsims)
net.out$obligate_ext <- net.out$obligate_ext-(nsims)
net.out$obligate_null_ext <- net.out$obligate_null_ext-(nsims)

net.out$number.species <- net.descr[as.character(net.out$net_id),"number_species"]

mod1 <- glmer(cbind(slope_ext,nsims*number.species-slope_ext)~tradeoff_strength + (tradeoff_strength|net_id),family=binomial,data=net.out)
mod2 <- glmer(cbind(horiz_ext,nsims*number.species-horiz_ext)~tradeoff_strength + (tradeoff_strength|net_id),family=binomial,data=net.out)
mod3 <- glmer(cbind(slope_null_ext,nsims*number.species-slope_null_ext)~tradeoff_strength + (tradeoff_strength|net_id),family=binomial,data=net.out)
mod4 <- glmer(cbind(horiz_null_ext,nsims*number.species-horiz_null_ext)~tradeoff_strength + (tradeoff_strength|net_id),family=binomial,data=net.out)

# Larger intecept values (less negative) are larger intercepts on figure
# Larger TS coefs (less negative) are shallower curves (more negative=steeper)
