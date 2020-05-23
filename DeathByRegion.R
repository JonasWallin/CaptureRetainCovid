###
# Simple function examining time series for each region
#
###

require(Matrix)
library(stringr)
source("stolen_function.R")

path.to.files <- file.path("data", "FHM")
fhm_files = list_fhm_files(folder = path.to.files)
death_dts <- c()
for(i in 1:length(fhm_files)) { 
  temp_frame = load_fhm_region(fhm_files[i])
  temp_frame$date = as.Date(fhm_files[i]%>% str_match_all("[0-9]{4}-[0-9]{2}-[0-9]{2}") %>%unlist)
  death_dts = rbind(death_dts,temp_frame)
  
}
death_VG <- death_dts[death_dts$Region=="Västra Götaland",c("TD","date")]
death_Skane <- death_dts[death_dts$Region=="Skåne",c("TD","date")]
death_Stockholm <- death_dts[death_dts$Region=="Stockholm",c("TD","date")]
death_dalarna <- death_dts[death_dts$Region=="Dalarna" ,c("TD","date")]
death_vast <- death_dts[death_dts$Region=="Västmanland" ,c("TD","date")]
 
x11()
par(mfrow=c(3,2))
plot(diff(death_VG$TD),xlab='day',ylab='rep', main = 'VG')
plot(diff(death_Skane$TD),xlab='day',ylab='rep', main='Skane')
plot(diff(death_Stockholm$TD),xlab='day',ylab='rep',main='Stockholm')
plot(diff(death_dalarna$TD),xlab='day',ylab='rep',main='Dalarna')
plot(diff(death_vast$TD),xlab='day',ylab='rep',main="Västmanland")