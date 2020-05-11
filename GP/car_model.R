##
# regular glm model for election prediction 
#
#
# D: 2018-07-08
##
graphics.off()
rm(list=ls())
library(rlist)
library(rstan)
rstan_options(auto_write = TRUE)
library(spdep)
library(maptools)
save.fig=T
source('./util/CARutil.R')

image_dir <- "../images/"
shp_file  <- "../data/shp/val_2014.shp"
nc.sids   <- readShapePoly(shp_file)
file_out <- c('../data/SD.csv') #store
y_name    <- c('sd_tal')
N_name    <- c('rost_gilti')
cov_name  <- c('med', 'skand', 'age')        
lans      <- unique(nc.sids@data$län) #c('Stockholms län','Blekinge län')   # nc.sids@data$län



nc.sids@data$med <- nc.sids@data$med/10^5
X <- nc.sids@data[,c(y_name, N_name, cov_name)]
X <- cbind(rep(1, dim(X)[1]) ,
           X)
names(X)[1] <- 'intercept'

fitted.models <- list()
for(i in 14:15){
  lan <- lans[i]
  cat('i = ',i,', lan = ', as.character( lans[i]), '\n')
  index <- nc.sids@data$län%in%lan==TRUE
  poly_in <- nc.sids[index,]
  X <- poly_in@data[,c(y_name, N_name, cov_name)]
  X <- cbind(rep(1, dim(X)[1]) ,
             X)
  names(X)[1] <- 'intercept'
  poly_in <- poly_in[complete.cases(X),]
  X       <- X[complete.cases(X),]
  if('med'%in%cov_name){
    poly_in <- poly_in[X$med>0,]
    X       <- X[X$med>0,]
  }
  
  fitted.models[[i]]  <-  model.stan(poly_in,
                                     X,
                                     y_name,
                                     N_name,
                                     path = './util/',
                                     iter   = 20000,
                                     warmup = 0.5)

  X <- X[, !( names( X) %in% c( y_name, N_name))]
  SumObj <- summary(fitted.models[[i]], probs=c(0.025,0.5,0.975))
  Betas <- SumObj$summary[1:dim(X)[2],c('2.5%','50%','97.5%')]
  out_res <- data.frame(gid = poly_in$gid)
  for(j in 1:dim(X)[2]){
    name_ <- paste(names(X)[j],'_025',sep='')
    out_res[,name_] <- Betas[j,'2.5%']
    name_ <- paste(names(X)[j],'_50',sep='')
    out_res[,name_] <- Betas[j,'50%']
    name_ <- paste(names(X)[j],'_975',sep='')
    out_res[,name_] <- Betas[j,'97.5%']
  }
  out_res[,'sigma_025'] <- SumObj$summary['sigma',c('2.5%')]
  out_res[,'sigma_50']  <- SumObj$summary['sigma',c('50%')]
  out_res[,'sigma_975'] <- SumObj$summary['sigma',c('97.5%')]
  out_res[,'tau_025']   <- SumObj$summary['tau',c('2.5%')]
  out_res[,'tau_50']    <- SumObj$summary['tau',c('50%')]
  out_res[,'tau_975']   <- SumObj$summary['tau',c('97.5%')]
  SumObj <- summary(fitted.models[[i]], pars = 'w',probs=c(0.025,0.5,0.975))
  objs <- c('w','p_coef','p_est')
  for(ii in 1:3){
    SumObj <- summary(fitted.models[[i]], pars = objs[ii],probs=c(0.025,0.5,0.975))
    out_res[, paste(objs[ii],'_025',sep='')] <- SumObj$summary[,'2.5%']
    out_res[, paste(objs[ii],'_50',sep='')] <- SumObj$summary[,'50%']
    out_res[, paste(objs[ii],'_975',sep='')] <- SumObj$summary[,'97.5%']
  }
  if(i==1){
  write.table(x=  out_res, 
              file =  file_out,
              append = F,
              quote  = F,
              sep = ';',
              row.names = F)
  }else{
    write.table(x=  out_res, 
                file =  file_out,
                append = T,
                quote  = F,
                sep = ';',
                row.names = F,
                col.names = F)
    
  }
  #list.save(fitted.models, '../../../temp/stored/save_SD.rds')
}
