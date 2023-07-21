library(lme4)
library(purrr)
library(doParallel)
library(plm)

# load the 80% metabolites
source('preprocessing.R')
metnames <- preprocessing(TRUE,'metnames')

# process imputed data metabolites from imputed data
data <- read_excel("PROCESSED_COVID19_METABOLOMICS_Common_SAMPLES_IMPUTED.xlsx",.name_repair = 'minimal')

# filter data
data <- data %>%
  # Remove samples from Cork
  filter(Location != 'Cork') %>%
  # remove on post-covid and control samples
  slice(which(!grepl('Long|control', COVID_SEVERITY))) %>%
  filter(`Days Post Hospitalisation` <9) 

# Filter on selected metabolites
carn <- which(colnames(data)=='carnitine')
metabolites <- data[,carn:ncol(data)]
metabolites <- metabolites %>% select(metnames)
data <- cbind(data[,1:(carn-1)], metabolites)

# Add sample numbers
sample =  sub('.*(?=.$)', '', data$COVID_SEVERITY, perl=T)
sample[grep('[^0-9]', sample)] <- "1"
data$sample <- as.numeric(sample)
data <- data %>% relocate(sample, .after = colnames(data)[carn-1])

# Remove patient numbers
data$Patient_ID <- sub('_[0-9]', '', data$Patient_ID)

# Remove patients with only one measurement
data <- data[data$Patient_ID %in% data$Patient_ID[duplicated(data$Patient_ID)],]

# remove cov-
data$Patient_ID <- as.numeric(sub('COV-', '', data$Patient_ID))

colnames(data) <- make.names(colnames(data))


# The script should also adjust for later death (geneva/stgallen) or moderate-severe (ticino) status
# if geneva/stgallen
data$Death.due.to.COVID.19[data$Location == 'Ticino' & data$COVID_SEVERITY == 'COVID_severe'] <- "Yes"
data$Death.due.to.COVID.19[data$Location == 'Ticino' & data$COVID_SEVERITY == 'COVID_moderate'] <- "No"

# Manually add the missing age and sex info <- This is added on 10-july-2023 however the metabolomics data already had this included when it was run.
data$AGE[data$Patient_ID=='COV-165'] <- 57
data$GENDER[data$Patient_ID=='COV-165'] <- 'Male'

data$AGE[data$Patient_ID=='COV-140'] <- 68
data$GENDER[data$Patient_ID=='COV-140'] <- 'Male'

data$AGE[data$Patient_ID=='COV-53'] <- 65
data$GENDER[data$Patient_ID=='COV-53'] <- 'Female'

data$AGE[data$Patient_ID=='COV-54'] <- 69
data$GENDER[data$Patient_ID=='COV-54'] <- 'Female'

data$AGE[data$Patient_ID=='COV-56'] <- 70
data$GENDER[data$Patient_ID=='COV-56'] <- 'Male'


logscaleLocation <- function (loc) {
  mets <- data %>%
    filter(Location == loc) %>%
    select(carnitine:X.26119) %>%
    apply(.,2,function(x) {log(x)})
  
  metadat <- data %>% select(1:ends_with("sample")) %>% filter(Location == loc)
  data <- cbind(metadat,mets)
}

# logscale data for each location separately
scaled_list <- lapply(unique(data$Location), logscaleLocation)
data1 <- do.call(rbind.data.frame, scaled_list)

# Get column names
metnamesOri <- metnames
metnames <- colnames(select(data1, carnitine:X.26119))



lmeres <- function(data, metnames, i) { 
  # get metabolite name
  met <- metnames[i]
  
  # Formula
  formula <- glue::glue('{met} ~  Days.Post.Hospitalisation + Death.due.to.COVID.19 + BMI + AGE + GENDER')
  # Run regression
  fit <- plm(formula, data=data, model="random", index=c("Patient_ID"))
  # Calculate 95% CI
  ci <- confint(fit)
  
  # concatinate results
  result <- summary(fit)$coef['Days.Post.Hospitalisation',]
  result <- c(ci['Days.Post.Hospitalisation',],result[c(1,2,3)])
}
lmeresAll <- function(i,data2,metnames) {lmeres(data2,metnames,i)}

getAllLocLme <- function (x) {
  # filter on location  
  data2 <- data1[data1$Location==unique(data1$Location)[x],]
  
  # calculate lme's
  res <- do.call(rbind, lapply(1:length(metnames),lmeresAll,data2,metnames))
  
  # process results
  res <- res %>% 
    as_tibble(res) %>% 
    mutate(metabolite = metnames) %>% 
    relocate(metabolite, .before = `2.5 %`) %>% 
    mutate(location = unique(data$Location)[x])
  return(res)
}

c1<-makeCluster(3)
registerDoParallel(c1)
res <- foreach(i=1:3, .combine = rbind, .packages = c('plm','tidyverse','purrr')) %dopar% {
  result <- getAllLocLme(i)
}

res$metabolite <- rep(metnamesOri,times = 3)
write.csv(res, 'lme_results_profile_cis_status.csv')

################################################################################
# Meta analysis of regression outcomes 

library("robumeta")
library("metafor")
library("dplyr")

# Read regression results
d <- read.csv('lme_results_profile_cis_status.csv')
# Select variables of interest
d <- d[, 2:ncol(d)]
colnames(d)[2:5] <- c("cil", "cir", "coef", "stderr")
d$var <- d$stderr^2

# Get results per location
dGeneva <- d[d$location == "Geneva", ]
dStGallen <- d[d$location == "St Gallen", ]
dTicino <- d[d$location == "Ticino", ]

# Preallocate data
dat <- data.frame(metabolite=character(),
                  estimateGeneva=double(),
                  estimateStGallen=double(),
                  estimateTicino=double(),
                  CI_Geneva=double(),
                  CI_StGallen=double(),
                  CI_Ticino=double(),
                  SE=double(),   
                  estimateMeta=double(),
                  estimatepval=double(),
                  estimateFDR=double(),
                  ConfInt=double(),
                  tau2=double(),
                  QE=double(),
                  Qpval=double(),
                  QFDR=double(),
                  stringsAsFactors=FALSE)

mets <- unique(d$metabolite)

# Perform meta-analysis
for (i in 1:length(mets)){
  
  dat[nrow(dat)+1, "metabolite"] <- mets[i]
  dat[nrow(dat), "estimateGeneva"] <- d$coef[d$location=="Geneva" & d$metabolite == mets[i]]
  dat[nrow(dat), "estimateStGallen"] <- d$coef[d$location=="St Gallen" & d$metabolite == mets[i]]
  dat[nrow(dat), "estimateTicino"] <- d$coef[d$location=="Ticino" & d$metabolite == mets[i]]
  
  dat[nrow(dat), "CI_Geneva"] <- paste("(", d$cil[d$location=="Geneva" & d$metabolite == mets[i]], ",", d$cir[d$location=="Geneva" & d$metabolite == mets[i]], ")", sep="")
  dat[nrow(dat), "CI_StGallen"] <- paste("(", d$cil[d$location=="St Gallen" & d$metabolite == mets[i]], ",", d$cir[d$location=="St Gallen" & d$metabolite == mets[i]], ")", sep="")
  dat[nrow(dat), "CI_Ticino"] <- paste("(", d$cil[d$location=="Ticino" & d$metabolite == mets[i]], ",", d$cir[d$location=="Ticino" & d$metabolite == mets[i]], ")", sep="")
  
  dprim <- d[d$metabolite == mets[i], ]
  res <- rma(coef, var, data=dprim, method="FE")
  
  dat[nrow(dat), "estimateMeta"]  <- round(res$beta[1,1], digits=5)
  dat[nrow(dat), "SE"]  <- round(res$se, digits=5)
  dat[nrow(dat), "ConfInt"]  <- paste("(", round(res$ci.lb, digits=5), ",", round(res$ci.ub, digits=5), ")", sep="")
  dat[nrow(dat), "tau2"]  <- round(res$tau2, digits=5)
  dat[nrow(dat), "estimatepval"]  <- res$pval
  dat[nrow(dat), "QE"]  <- round(res$QE, digits=5)
  dat[nrow(dat), "Qpval"]  <- res$QEp
  
}

# Add p-value corrections
dat$estimateFDR <- p.adjust(dat$estimatepval, method = "fdr", n = length(dat$estimatepval))
dat$QFDR <- p.adjust(dat$Qpval, method = "fdr", n = length(dat$Qpval))

df = subset(dat, select = -c(13) )
dfnew <- df[which(df$estimateFDR < 0.05 & df$QFDR >0.05),]

# Write results
write.csv(df, "CovidMetaanalysisFE_Estimates.csv'.csv", row.names = FALSE)