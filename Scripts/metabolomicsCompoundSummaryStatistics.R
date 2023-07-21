source('preprocessing.R')
data <- preprocessing(FALSE)
metnames <- preprocessing(FALSE,'metnames')

# Get all metabolite names
metabolites <- data[,which(colnames(data) %in% metnames)]

metabolites <- cbind(data$Location, metabolites)
colnames(metabolites)[1]<-'Location'

f <- function(x) {sum(is.na(x))/length(x)}

missing_values <- metabolites %>% 
  group_by(Location) %>% 
  group_map(~ apply(., MARGIN=2, FUN=f)) %>% 
  do.call(cbind.data.frame, .)

missing_values$total <- apply(metabolites[,2:ncol(metabolites)], MARGIN=2, FUN=f)

# Make new column names
colnames(missing_values)<- c('Geneva', 'St Gallen', 'Ticino', 'Total')
missing_values <- as_tibble(missing_values,rownames='metabolite')

# Add a column denoting if the metabolite is removed
missing_values$metabolite_removed <- missing_values$Total>0.2

# save results
write_csv(missing_values,'missing_values.csv')
