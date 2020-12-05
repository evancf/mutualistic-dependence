#####
##### Functions for network simulations, model output analysis
#####


##### Define misc functions, model output functions

'%!in%' <- function(x,y)!('%in%'(x,y))

rel1 <- function(x){
  x/max(x)
}

logrel1 <- function(x){
  (log(x)-min(log(x)))/max(log(x)-min(log(x)))
}

# For use in lme4 predictions

random.lines.logit <- function(x,col=1,lwd=1,factor=1){
  y <- as.data.frame(coef(x)[1])
  for(i in 1:dim(y)[1]){
    curve(plogis(y[i,1] + y[i,2]*x)*factor, add=T, col=col)
  } 
}

fixed.lines.logit <- function(x,col=1,lwd=1,factor=1,...){
  y <- fixef(x)
  curve(plogis(y[1] + y[2]*x)*factor, add=T, lwd=lwd, col=col) 
}

random.pred.logit <- function(mod,val){
  y <- as.data.frame(coef(mod)[1])
  vec <- NULL
  ran.names <- rownames(coef(mod)$net_id)
  for(i in 1:dim(y)[1]){
    vec <- c(vec,as.numeric(plogis(y[i,1] + y[i,2]*val)))
  } 
  out <- as.data.frame(cbind(vec,
                             rep(ran.names,each=length(val)),
                             rep(val,times=length(ran.names))))
  out[,1] <- as.numeric(as.character(out[,1]))
  out[,3] <- as.numeric(as.character(out[,3]))
  out
}

fixed.pred.logit <- function(mod,val){
  y <- fixef(mod)
  return(as.numeric(plogis(y[1] + y[2]*val)))
}



##### Define netcascade functions
# Vieira and Almeida-Neto Ecology Letters (2015) 18: 144–152

source("netcascade (April 2014).R")
