###
#
#
# funktioner tagna eller inspirerade av:
# <https://github.com/adamaltmejd/covid>
##
require(readxl)
require(ggplot2)
require(hrbrthemes)
require(gdtools)

download_latest_fhm <- function(folder = file.path("data", "FHM")) {

  
  DL <- download.file("https://www.arcgis.com/sharing/rest/content/items/b5e7488e117749c19881cce45db13f7e/data",
                      destfile = file.path(folder, "FHM_latest.xlsx"), method = "curl", extra = c("-L"), quiet = TRUE)
  if (DL != 0) { stop("File download error.") }
  
  # Check archived files for latest record
  latest_record <- max(as.Date(gsub("^.*(2020-[0-9]{2}-[0-9]{2}).xlsx", "\\1", list.files(folder, pattern = "^Folkhalso"))))
  
  # Check if new download is newer than latest record, in that case, archive it.
  new_record <- get_record_date(file.path(folder, "FHM_latest.xlsx"))
  
  if (latest_record < new_record) {
    file.copy(file.path(folder, "FHM_latest.xlsx"),
              file.path("data", "FHM", paste0("Folkhalsomyndigheten_Covid19_", new_record, ".xlsx")))
  }
}

list_fhm_files <- function(folder = file.path("data", "FHM")) {
  list.files(folder, pattern = "^Folkhalso", full.names = TRUE)
}


load_fhm <- function(f) {
  require(data.table)
  require(readxl)
  
  DT <- data.table((
    read_excel(path = f, sheet = 2, col_types = c("text", "numeric"))
  ))
  
  setnames(DT, c("date", "N"))
  
  DT[date == "Uppgift saknas" | date == "uppgift saknas", date := NA]
  
  if (can_be_numeric(DT[, date])) {
    DT[, date := as.Date(as.numeric(date), origin = "1899-12-30")]
  } else {
    DT[, date := as.Date(date)]
  }
  
  DT[is.na(N), N := 0]
  
  DT[, publication_date := get_record_date(f)]
  
  return(as.data.frame(DT))
}


load_fhm_region <- function(f) {
  require(data.table)
  require(readxl)
  
  DT <- data.table((
    read_excel(path = f, sheet = 4, col_types = c("text", "numeric","numeric","numeric","numeric"))
  ))
  
  setnames(DT, c("Region", "TA","TAp","TI","TD"))
  

  return(as.data.frame(DT))
}

can_be_numeric <- function(x) {
  # Check if vector can be converted to numeric
  stopifnot(is.atomic(x) || is.list(x)) # check if x is a vector
  numNAs <- sum(is.na(x))
  numNAs_new <- suppressWarnings(sum(is.na(as.numeric(x))))
  return(numNAs_new == numNAs)
}



get_record_date <- function(f) {
  sheets <- excel_sheets(f)
  return(as.Date(sub("^FOHM ", "", sheets[length(sheets)]), format="%d %b %Y"))
}


set_default_theme <- function() {
  require(hrbrthemes)
  
  fam <- "sans"
  if (font_family_exists(font_family = "Arial")) fam <- "Arial"
  if (font_family_exists(font_family = "EB Garamond")) fam <- "EB Garamond"
  
  theme_ipsum(base_family = fam) %+replace%
    theme(
      plot.title = element_text(size = rel(2), face = "plain", hjust = 0, margin = margin(0,0,5,0)),
      plot.subtitle = element_text(size = rel(1), face = "plain", hjust = 0, margin = margin(0,0,5,0)),
      legend.background = element_rect(fill = "grey90", color = "grey80"),
      legend.margin = margin(5,5,5,5),
      legend.direction = "vertical",
      legend.position = "right",
      axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1.2),
      
      # Panels
      plot.background = element_rect(fill = "#f5f5f5", color = NA), # bg of the plot
      panel.border = element_blank(),
      panel.grid.major = element_line(linetype = "dotted", color = "grey60", size = 0.2),
      panel.grid.minor = element_line(linetype = "dotted", color = "grey80", size = 0.2)
    )
}

plot.predReport <- function(result, CI,  true.day=0, ymax = NULL){
  
  Reported <- result$detected
  N <- dim(Reported)[2]
  reported_so_far <- c()
  for(i in 1:N){
    reported_so_far <- c(reported_so_far,max(Reported[i,i:N]))
  }
  default_theme = set_default_theme()
  pred.data <- data.frame(date = c(result$dates),
                          deaths =reported_so_far,
                          CIl      = CI[1,],
                          CIu      = CI[2,])
  if(true.day>0){
    pred.data$type     = as.factor(c(rep("låst", true.day), 
                                      rep("ej låst", N-true.day)))
    ggfig <- ggplot(data = pred.data, aes(y = deaths, x = date, fill = type)) +
      geom_bar( stat="identity") + 
      default_theme+
      scale_x_date(date_breaks = "4 day", expand = c(0, 0)) +
      scale_fill_manual(values=c(  alpha("gray",0.8),"red")) + 
      geom_errorbar(aes(ymin = CIl, ymax = CIu), width = 1) +
      guides(fill=guide_legend(title="låsta dagar"))
  }else{
    ggfig <- ggplot(data = pred.data, aes(y = deaths, x = date)) +
      geom_bar( stat="identity") + 
      default_theme+
      scale_x_date(date_breaks = "4 day", expand = c(0, 0)) +
      scale_fill_manual(values=c(  alpha("gray",0.8))) + 
      geom_errorbar(aes(ymin = CIl, ymax = CIu), width = 1)
  }
  if(is.null(ymax)==F)
  {
    ggfig <- ggfig + coord_cartesian(ylim = c(0,ymax))
  }
  return(ggfig)
}
#'
#' Here we remove the fixed day and only display reported days
#' CI and median 
#' CI  - (3 x N) lower CI, median, upper CI,
#' 
plot.predReport2 <- function(result, CI,  ymax = NULL){
  
  Reported <- result$detected
  N <- dim(Reported)[2]
  reported_so_far <- c()
  for(i in 1:N){
    reported_so_far <- c(reported_so_far,max(Reported[i,i:N]))
  }
  default_theme = set_default_theme()
  pred.data <- data.frame(date = c(result$dates),
                          deaths =reported_so_far,
                          CIl      = CI[1,],
                          median   = CI[2,],
                          CIu      = CI[3,])

    ggfig <- ggplot(data = pred.data, aes(y = deaths, x = date)) +
      geom_bar( stat="identity") + 
      default_theme+
      scale_x_date(date_breaks = "4 day", expand = c(0, 0)) +
      scale_fill_manual(values=c(  alpha("gray",0.8))) + 
      geom_errorbar(aes(ymin = CIl, ymax = CIu), width = 1)
    ggfig <- ggfig  + geom_line(
                            mapping=aes(y = median, x = date),
                            lwd=1.,
                            color='black')+ scale_color_discrete(name = "Y series", labels = c("Y2"))
  if(is.null(ymax)==F)
    ggfig <- ggfig + coord_cartesian(ylim = c(0,ymax))
  
  return(ggfig)
}

#'
#' Here we remove the fixed day and only display reported days
#' CI and median 
#' CI  - (3 x N) lower CI, median, upper CI,
#' 
plot.predReport3 <- function(result, CI,  ymax = NULL){
  
  Reported <- result$detected
  N <- dim(Reported)[2]
  reported_so_far <- c()
  for(i in 1:N){
    reported_so_far <- c(reported_so_far,max(Reported[i,i:N]))
  }
  default_theme = set_default_theme()
  pred.data <- data.frame(date = c(result$dates),
                          reports =reported_so_far,
                          CIl      = CI[1,],
                          median   = CI[2,],
                          CIu      = CI[3,])
  
  ggfig <- ggplot(data = pred.data, aes(y = reports, x = date)) +
    geom_bar( stat="identity") + 
    default_theme+
    scale_x_date(date_breaks = "4 day", expand = c(0, 0)) +
    scale_fill_manual(values=c(  alpha("gray",0.8))) + 
    geom_errorbar(aes(ymin = CIl, ymax = CIu), width = 1)
  if(is.null(ymax)==F)
    ggfig <- ggfig + coord_cartesian(ylim = c(0,ymax))
  
  return(ggfig)
}