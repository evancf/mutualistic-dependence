#(April 23 2014).

netcascade <- function(imatrix,ranim,rplants,deadPlants=NULL, deadAnimals=NULL, targetGuild,target,return.matrix=F){
  
  #---------ARGUMENT CHECKS-----------------------------------
  
  if(class(imatrix)!="matrix" || (nrow(imatrix)+ncol(imatrix))<3){stop("'imatrix' must be an object of class 'matrix', with animal species on rows, plant species on columns and at least three species overall")}
  if(class(ranim)!="numeric" || class(rplants) != "numeric" || max(c(max(ranim),max(rplants)))>1 || min(c(min(ranim),min(rplants)))<0){stop("'ranim' & 'rplants' must be numeric vectors with values ranging between 0 and 1")}
  if((targetGuild%in%c("animal","plant"))==F){stop('Invalid target guild for primary extinction. Valid targets guilds are "animal" and "plant"')}
  if(is.numeric(target)==F){stop('Invalid value for the "target" argument. You may specify a single species by entering its row or column number or you may use a vector of relative probabilites for all species in the target guild.')}
  if(is.null(deadAnimals)==F && class(deadAnimals)!= "integer"){stop("deadAnimals must be either NULL or an integer vector specifying the row numbers of animals considered to be extinct on the original matrix")}
  if(is.null(deadPlants)==F && class(deadPlants)!= "integer"){stop("deadPlants must be either NULL an integer vector specifying the column numbers of plants considered to be extinct on the original matrix")}
  if(length(ranim)!= nrow(imatrix)){stop("The length of vector'ranim' must be equal to number of rows (i.e. animal species) in 'imatrix'")}
  if(length(rplants)!= ncol(imatrix)){stop("The length of vector'rplants' must be equal to number of columns (i.e. plant species) in 'imatrix'")}
  
  #---------DEFINING SOME VARIABLES---------------------------
  nanim <- nrow(imatrix)
  npla <- ncol(imatrix)
  plants <- 1:npla;
  animals <- 1:nanim
  plantNA <- 1:npla
  animNA <- 1:nanim
  plantNA[deadPlants] <- NA
  animNA[deadAnimals] <- NA 
  
  degree_when_lost_plants <- c()
  degree_when_lost_animals <- c()
  
  #----------CALCULATING DEPENDENCE MATRICES-------------------
  M <- array(0,dim=c(nanim,npla,2))
  for(i in 1:npla){
    M[,i,1] <- imatrix[,i]/sum(imatrix[,i])
  } #matrix of plant dependence on each animal
  for(i in 1:nanim){
    M[i,,2] <- imatrix[i,]/sum(imatrix[i,])
  } #matrix of animal dependence on each plant
  
  #-----------CHOOSING TARGET SPECIES FOR PRIMARY EXTINCTION---
  coext_animals <- c()
  coext_plants <- c()
  if(length(target)==1){
    if(targetGuild=="animal"){
      if(target %in% deadAnimals){stop('Specified target species for the primary extinction is already extinct')}
      coext_animals <- target
      degree_when_lost_animals <- 1 #stores the degree of the extinction event of every animal species lost during the coextinction cascade. 
    }
    if(targetGuild=="plant"){
      if(target %in% deadPlants){stop('Specified target species for the primary extinction is already extinct')}
      coext_plants <- target
      degree_when_lost_plants <- 1
    }
  }else{
    nspecies <- switch(targetGuild,animal = nanim, plant = npla)
    if(length(target)==nspecies){
      if(targetGuild =="animal"){
        alive <- animals[is.na(animNA)==F]
        coext_animals <- sample(c(alive,0),1,prob = c(target[is.na(animNA)==F],0))
        degree_when_lost_animals <- 1
      }
      if(targetGuild =="plant"){
        alive <- plants[is.na(plantNA)==F]
        coext_plants <- sample(c(alive,0),1,prob = c(target[is.na(plantNA)==F],0))
        degree_when_lost_plants <- 1
      }
    }else{
      stop('Length of "target" must be 1 (specifying a single species within the target guild) or else be equal to the number of species in the target guild (specifying probabilities of primary extinction for each species in the target guild)')
    }
  }
  imatrix[coext_animals,] <- 0 
  imatrix[,coext_plants] <- 0
  
  lostanimals <- coext_animals #final list of animals which were "alive" in the original community but became extinct during this primary extinction + extinction cascade
  lostplants <- coext_plants 
  
  #-------------------CASCADE LOOP---------------------------
  equilibrium <- FALSE
  degree <- 1
  degree_table <- data.frame(degree,guild=factor(targetGuild,levels=c("animal","plant")),n_extinctions=1)
  while(equilibrium == FALSE){
    ext_animals <- coext_animals
    ext_plants <- coext_plants
    plantNA[ext_plants] <- NA
    animNA[ext_animals] <- NA
    
    aleft <- animals[is.na(animNA)==F]
    pleft <- plants[is.na(plantNA)==F]
    
    if(length(ext_animals)>0){
      for(i in 1:length(ext_animals)){
        unlucky <- rplants[pleft]*M[ext_animals[i],pleft,1] > runif(length(pleft))
        coext_plants = c(coext_plants,pleft[unlucky])
      }
      coext_animals <- c()
      coext_plants <- unique(coext_plants)
      plantNA[coext_plants] <- NA
      lostplants <- c(lostplants,coext_plants)
      imatrix[,coext_plants] <- 0
      for(i in 1:npla){
        if(sum(imatrix[,i])==0){
          M[,i,1] <- 0
        }else{
          M[,i,1] <- imatrix[,i]/sum(imatrix[,i])
        }
      }
      if(length(coext_plants)>0){
        degree <- degree + 1
        degree_when_lost_plants <- c(degree_when_lost_plants, rep(degree,length(coext_plants)))
        degree_table[degree,] <- data.frame(degree,"plant",length(coext_plants))
      }
    }else{
      for(i in 1:length(ext_plants)){
        unlucky <- ranim[aleft]*M[aleft,ext_plants[i],2] > runif(length(aleft))
        coext_animals <- c(coext_animals, aleft[unlucky])
      }
      coext_plants = c();
      coext_animals <- unique(coext_animals)
      lostanimals <- c(lostanimals,coext_animals)
      animNA[coext_animals] <- NA
      imatrix[coext_animals,] <- 0
      for(i in 1:nanim){
        if(sum(imatrix[i,])==0){
          M[i,,2] <- 0
        }else{
          M[i,,2] <- imatrix[i,]/sum(imatrix[i,])
        }
      }
      if(length(coext_animals)>0){
        degree <- degree + 1
        degree_when_lost_animals <- c(degree_when_lost_animals, rep(degree,length(coext_animals)))
        degree_table[degree,] <- data.frame(degree,"animal",length(coext_animals))  
      }
    }
    equilibrium <- equilibrium + (length(coext_plants)+length(coext_animals))==0
  }
  
  #-------------------OUTPUT---------------------------
  
  if(return.matrix==T){
    return(list(interaction_matrix = imatrix, lost_animals = lostanimals, lost_plants = lostplants))  
  }else{
    if(length(lostanimals)>0){
      spp_data_animals <- data.frame(lost_animal = lostanimals,degree_of_extinction=degree_when_lost_animals)
    }else{
      spp_data_animals <- "No animal species were lost"
    }
    if(length(lostplants)>0){
      spp_data_plants <- data.frame(lost_plant = lostplants,degree_of_extinction=degree_when_lost_plants)
    }else{
      spp_data_plants <- "No plant species were lost"
    }
    return(list(cascade_data=degree_table,animal_species_data=spp_data_animals,plant_species_data=spp_data_plants))
  }
}
