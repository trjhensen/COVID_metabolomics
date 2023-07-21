# Metabolite depletion analysis and meta-analysis (Not included in the manuscript)
################################################################################
# Linear mixed-effect regressions

library(lme4)
library(purrr)
library(doParallel)
library(here)

# load the 80% metabolites
source('preprocessing.R')
data <- preprocessing(FALSE)
colnames(data) <- make.names(colnames(data))

lmeres <- function(data1, metnames, i) { 
  # get metabolite name
  met <- metnames[i]
  # Formula
  formula <- glue::glue('{met} ~  Days.Post.Hospitalisation + BMI + AGE + GENDER + (1 | Patient_ID)')
  
  # Calculate confidence interval for the time points
  cn <- c('2.5 %','97.5 %','Estimate','Std. Error','z value','Pr(>|z|)', 'singular_fit')
  err <- matrix(c(NaN,NaN,NaN,NaN,NaN,NaN,NaN), nrow = 1, ncol = 7, dimnames = list(c(""),cn))
  
  trycat = possibly(.f = glmer, otherwise = err)
  
  fit <- trycat(formula, data = data1, family = binomial)
  
  if(any(class(fit) == 'matrix')) {return(err)}
  else {
    
    # get results
    result <- summary(fit)$coef['Days.Post.Hospitalisation',]
    
    # get confidence interval
    ci <- confint(fit,method ="Wald")
    
    # concatinate results
    result <- c(ci['Days.Post.Hospitalisation',],result[c(1,2,3,4)])
    
    # check for singularity
    result['singular_fit'] <- isSingular(fit, tol = 1e-4)
    
    return(result)
  }
  
}
lmeresAll <- function(i,data2,metnames) {lmeres(data2,metnames,i)}
getAllLocLme <- function (x) {
  # filter on location  
  data2 <- filter(data,Location==unique(data$Location)[x])
  
  f <- function() apply(!is.na(metabolites), 2, sum)/dim(metabolites)[1]
  # filter on metabolites with between 20 and 80% missing values
  metabolites <- data2 %>%  select(carnitine:X.26119)
  #metabolites <- metabolites %>% select(which(f()>.2 & f()<.8))
  # convert metabolite values into true (measured) and false (NA)
  metabolites <- is.na(metabolites)
  data2 <- cbind(data2[,1:10],metabolites)
  
  metnames <- colnames(metabolites)
  
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
res <- do.call(rbind, lapply(1:3,getAllLocLme))

# save results
write.csv(res, 'log_glme_results_wald_cis.csv')

# Get statistics on the number of compounds for which a result could be obtained
singCounts <-res %>% 
  group_by(location) %>% 
  summarise(count = sum(singular_fit,na.rm = TRUE),
            length = length(singular_fit), 
            fraction_singular = count/length,
            errors = sum(is.na(singular_fit)))
write.csv( singCounts, 'log_glme_number_of_singular_fits_wald.csv')

################################################################################
# Meta-analysis

# Load data and restore metabolite names
depletions <- read.csv('log_glme_results_wald_cis.csv')

# Load missing values file
miss_vals <- read.csv('missing_values.csv') # Create in script getMissingMetabolitePercentages.R
# Get metabolite names and replicate three times and restore metnames
depletions$metabolite <- rep(miss_vals$metabolite, times = 3)

# filter on metabolites with no NAs in all locations

# Check if there is a singular fit
depletions$singular_fit[is.na(depletions$singular_fit)] <- 1

# for each metabolite, if z-value != NA, add to list

# Get all metabolite names
metabolites <- unique(depletions$metabolite)

# Create function to find which metabolites should be kept
addMetKeep <- function(depletions,metabolites,i) {
  met <- metabolites[i]
  keep <- depletions %>% 
    filter(metabolite == met) %>%
    mutate(keep = !any(is.na(z.value)))
  return(keep)
}

# Add TRUE or false to find which metabolites should be removed
depKeep <- foreach(i=1:length(metabolites), .combine = rbind, .packages = 'dplyr') %do% {
  result <- addMetKeep(depletions,metabolites,i)
}

# Remove all metabolites with FALSE
deplFilt <- depKeep[depKeep$keep,]

# Process dataframe
deplFilt <- select(deplFilt, -c('X','keep'))
colnames(deplFilt)[2:5] <- c("cil", "cir", "coef", "stderr")

# Add variance
deplFilt$var <- deplFilt$stderr^2

Check which metabolites have a singular value in at least one location
# for each metabolite, check if the value is singular
mets <- unique(deplFilt$metabolite)

checkSingRes <- function(i){
  result <- deplFilt %>% 
    filter(metabolite == mets[i]) %>%
    summarise(any(singular_fit==1))
  return(result)
}
singVal <- as.numeric(do.call('rbind',sapply(1:length(mets),checkSingRes)))


Create result table and perform meta analysis
# result table
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
                  singularVal_Present = double(),
                  stringsAsFactors=FALSE)


for (i in 1:length(mets)){
  
  dat[nrow(dat)+1, "metabolite"] <- mets[i]
  dat[nrow(dat), "estimateGeneva"] <- deplFilt$coef[deplFilt$location=="Geneva" & deplFilt$metabolite == mets[i]]
  dat[nrow(dat), "estimateStGallen"] <- deplFilt$coef[deplFilt$location=="St Gallen" & deplFilt$metabolite == mets[i]]
  dat[nrow(dat), "estimateTicino"] <- deplFilt$coef[deplFilt$location=="Ticino" & deplFilt$metabolite == mets[i]]
  
  dat[nrow(dat), "CI_Geneva"] <- 
    paste("(", deplFilt$cil[deplFilt$location=="Geneva" & deplFilt$metabolite == mets[i]], ",", 
          deplFilt$cir[deplFilt$location=="Geneva" & deplFilt$metabolite == mets[i]], ")", sep="")
  
  dat[nrow(dat), "CI_StGallen"] <- 
    paste("(", deplFilt$cil[deplFilt$location=="St Gallen" & deplFilt$metabolite == mets[i]], ",", 
          deplFilt$cir[deplFilt$location=="St Gallen" & deplFilt$metabolite == mets[i]], ")", sep="")
  
  dat[nrow(dat), "CI_Ticino"] <- 
    paste("(", deplFilt$cil[deplFilt$location=="Ticino" & deplFilt$metabolite == mets[i]], ",", 
          deplFilt$cir[deplFilt$location=="Ticino" & deplFilt$metabolite == mets[i]], ")", sep="")
  
  dprim <- deplFilt[deplFilt$metabolite == mets[i], ]
  res <- rma(coef, var, data=dprim, method="FE")
  
  dat[nrow(dat), "estimateMeta"]  <- round(res$beta[1,1], digits=5)
  dat[nrow(dat), "SE"]  <- round(res$se, digits=5)
  dat[nrow(dat), "ConfInt"]  <- paste("(", round(res$ci.lb, digits=5), ",", round(res$ci.ub, digits=5), ")", sep="")
  dat[nrow(dat), "tau2"]  <- round(res$tau2, digits=5)
  dat[nrow(dat), "estimatepval"]  <- res$pval
  dat[nrow(dat), "QE"]  <- round(res$QE, digits=5)
  dat[nrow(dat), "Qpval"]  <- res$QEp
  
}

# Add fdr values
dat$estimateFDR <- p.adjust(dat$estimatepval, method = "fdr", n = length(dat$estimatepval))
dat$QFDR <- p.adjust(dat$Qpval, method = "fdr", n = length(dat$Qpval))

# Add singular value
dat$singularVal_Present<- singVal

# remove tau2
df = subset(dat, select = -c(13) )
# Get all metabolites that are consistent
dfnew <- df[which(df$estimateFDR < 0.05 & df$QFDR >0.05),]

# Get missing values file
miss_vals <- read.csv('../results/missing_values.csv')
# select columns
miss_vals <- miss_vals[,c('metabolite','Total','metabolite_removed')]

# merge tables
df_full <- merge(df,miss_vals, by = 'metabolite',all.x = TRUE)

# Save results
write.csv(df_full, "depletion_metaAnalyses.csv", row.names = FALSE)
