preprocessing <- function(remove_20perc_mets,output='') {
  
# remove_20perc_mets (TRUE/FALSE) 
#       TRUE: Metabolites found in less than 80% of samples are removed
#       FALSE: No metabolites are removed
# output: output == metnames -> only the metabolite names are used as output
#         output != metnames -> the full dataset is used as output
  
# Filter and process data set for further analysis 

# load library 
library(tidyverse)
library(readxl)

# load original data
  getwd()
data <- read_excel("C:/Users/mspg/National University of Ireland, Galway/Group_MSP - COVID-19 patient data/PROCESSED_COVID19_METABOLOMICS_Common_SAMPLES_RAW_NO_IMP.xlsx")

# filter data
data <- data %>%
  # Remove samples from Cork
  filter(Location != 'Cork') %>%
  # remove on post-covid and control samples
  slice(which(!grepl('Long|control', COVID_SEVERITY))) %>%
  filter(`Days Post Hospitalisation` <9)

# Get first metabolite column
carn <- which(colnames(data)=='carnitine')

# Add sample numbers
sample =  sub('.*(?=.$)', '', data$COVID_SEVERITY, perl=T)
sample[grep('[^0-9]', sample)] <- "1"
data$sample <- as.numeric(sample)
data <- data %>% relocate(sample, .after = colnames(data)[carn-1])

# Remove patient numbers
data$Patient_ID <- sub('_[0-9]', '', data$Patient_ID)



# Get all metabolites in the data
metabolites <- data[,(carn+1):ncol(data)]

# return metabolites
if (output=='metnames' & !remove_20perc_mets) { 
  return(colnames(metabolites))
} else if (output!='metnames' & !remove_20perc_mets) {
  return(data)
}

if (remove_20perc_mets) {
  # Remove all metabolites that were present in less then 20% of the total samples
  eightyperc <- function(x) {apply(is.na(x), 2, sum)/dim(x)[1] >= .2}
  metabolites1<-metabolites[,-which(eightyperc(metabolites))]
  
  if (output=='metnames') {
    return(colnames(metabolites1))  
  }

  # Add new metabolites to the data variable
  data <- cbind(data[,1:(carn)],metabolites1)
  return(data)
}

# OUTPUT data
return(data)
}