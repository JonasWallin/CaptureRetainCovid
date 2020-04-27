###
# load data
##

require(Matrix)
source("stolen_function.R")
download_latest_fhm()
path.to.files <- file.path("data", "FHM")
fhm_files = list_fhm_files(folder = path.to.files)
death_dts <- c()
for(i in 1:length(fhm_files)) 
  death_dts = rbind(death_dts,load_fhm(fhm_files[i]))

# res <- death_dts %>% dplyr::group_by(date, publication_date) %>% dplyr::spread(publication_date, value=N)
res <- death_dts %>% dplyr::group_by(date, publication_date)
death_dt <- data.table(death_dts)
setkey(death_dt, publication_date, date)
death_dt <-death_dt[!is.na(date) & publication_date > "2020-04-02" & date > "2020-04-02"]
res2 <- death_dt%>%tidyr::spread(publication_date, N)
detected <- as.matrix(res2[,2:dim(res2)[2]])
detected <- cbind(matrix(NA,dim(detected)[1],dim(detected)[1]-dim(detected)[2] ) , detected)
detected[lower.tri(detected)]=NA
colnames(detected) <- c(res2[,1]$date)
result = list(detected = detected, dates = res2[,1]$date, dates_report=unique(death_dt$publication_date))
save(result, file = "result.RData")

