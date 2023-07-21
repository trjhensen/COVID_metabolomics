################################################################################
# Figure 1A
# Metabolite distributions of significant compounds
# load libraries
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(ggrepel)
library(stringr)

# Get full data metabolite names and super-pathway info
metNames <- read_xlsx('APCM-02-21MD+ MERGED DATA TABLES_Correct.xlsx',
                      sheet = 'Chemical Annotation Common')

# rename super pathway in Biochemical family
metNames <- rename(metNames, `Biochemical family` = SUPER_PATHWAY)

# CONVERT `Biochemical family` names
metNames$`Biochemical family` <- gsub('Amino Acid','Amino Acids',metNames$`Biochemical family`)
metNames$`Biochemical family` <- gsub('Carbohydrate','Carbohydrates',metNames$`Biochemical family`)
metNames$`Biochemical family` <- gsub('Nucleotide','Nucleotides',metNames$`Biochemical family`)
metNames$`Biochemical family` <- gsub('Energy','Energy metabolism',metNames$`Biochemical family`)
metNames$`Biochemical family` <- gsub('Lipid','Lipids',metNames$`Biochemical family`)
metNames$`Biochemical family` <- gsub('Peptide','Peptides',metNames$`Biochemical family`)
metNames$`Biochemical family`[which(is.na(metNames$`Biochemical family`))] <- 'Unnamed Compounds'

# Read results
regRes <- read.csv('CovidMetaanalysisFE_Estimates.csv') %>% filter(QFDR>0.05 & estimateFDR<0.05)
regMets <- regRes$metabolite

sig <- metNames %>% filter(CHEMICAL_NAME  %in% regMets)
sig1<- sig %>% group_by(`Biochemical family`, SUB_PATHWAY) %>% summarise(sig_ab = n()) %>% arrange(desc(sig_ab))

source('./scripts/preprocessing.R')
eightypercmets <- preprocessing(TRUE, 'metnames')
analysed <- metNames %>% filter(CHEMICAL_NAME  %in% eightypercmets)
an1<- analysed %>% group_by(`Biochemical family`, SUB_PATHWAY) %>% summarise(an_ab = n()) %>% arrange(desc(an_ab))

total <- merge(sig1,an1) %>% 
  mutate(expected = an_ab * (sum(sig_ab)/sum(an_ab))) %>%
  mutate(enrichment_ratio = (sig_ab / expected)) 
total$SUB_PATHWAY[is.na(total$SUB_PATHWAY)]<-'Unnamed Compounds'

s = 6
total %>% 
  filter(`Biochemical family` != 'unnamed metabolites' & `Biochemical family` != 'Partially Characterized Molecules') %>%
  arrange(`Biochemical family`,sig_ab) %>%
  ggplot() + 
  geom_bar(aes(x = as_factor(SUB_PATHWAY), y=an_ab, fill = ' Analysed metabolites per pathway'),stat = 'identity')+
  geom_bar(aes(x = as_factor(SUB_PATHWAY), y=sig_ab, fill=`Biochemical family`),stat = 'identity') +
  labs(y = 'Metabolites', 
       x = '',
       title = 'Distribution of significant metabolites per biochemical family') +
  scale_fill_manual(values=c("#bcbcbc",colorRampPalette(brewer.pal(8, "Spectral"))(8))) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold',hjust = 0.5),
    axis.text.x = element_text(size = s*1.3,angle = 30, hjust = 1, colour = 'black'),
    axis.text.y = element_text(size = s*2, face = 'bold'),
    legend.text = element_text(size = s*2),
    legend.title = element_text(face = 'bold')
  )

png(file="BarPathwayDistr.png", width=12, height=12, units="in", res=300)
s = 6
total %>%
  #filter(`Biochemical family` != 'unnamed metabolites' & `Biochemical family` != 'Partially Characterized Molecules') %>%
  arrange(`Biochemical family`,desc(sig_ab)) %>%
  ggplot() +
  geom_col(aes(y = fct_rev(as_factor(SUB_PATHWAY)), x=an_ab, fill = 'All analysed metabolites'))+
  geom_col(aes(y = fct_rev(as_factor(SUB_PATHWAY)), x=sig_ab, fill=`Biochemical family`), col = 'black')+
  labs(x = 'Metabolites',
       y = '',
       title = 'Number of serum metabolites that changed over time') +
  scale_fill_manual(values=c("#bcbcbc",colorRampPalette(brewer.pal(8, "Spectral"))(10))) +
  theme_bw() +
  theme(
    plot.title = element_text(size = s*3, face = 'bold',hjust = 0.5),
    axis.text.y = element_text(size = s*1.5, face = 'bold'),
    axis.text.x = element_text(size = s*2, face = 'bold'),
    axis.title.x = element_text(size = s*2, face = 'bold'),
    legend.text = element_text(size = s*1.6),
    legend.title = element_text(size = s*2.5, face = 'bold'),
    legend.position = "bottom") +
  guides(fill=guide_legend(title=""),
         color=guide_legend(nrow=2, byrow=FALSE))
dev.off()
################################################################################
# Figure 1B,1C,1D,and 1E
# load libraries
library(tidyverse)
library(readxl)
library(RColorBrewer)

# Get full data metabolite names and super-pathway info
metNames <- read_xlsx('APCM-02-21MD+ MERGED DATA TABLES_Correct.xlsx',
                      sheet = 'Chemical Annotation Common')

# rename super pathway in biochemical class
metNames <- rename(metNames, `Biochemical class` = SUPER_PATHWAY)

# CONVERT `Biochemical class` names
metNames$`Biochemical class` <- gsub('Amino Acid','Amino Acids',metNames$`Biochemical class`)
metNames$`Biochemical class` <- gsub('Carbohydrate','Carbohydrates',metNames$`Biochemical class`)
metNames$`Biochemical class` <- gsub('Nucleotide','Nucleotides',metNames$`Biochemical class`)
metNames$`Biochemical class` <- gsub('Energy','Energy metabolism',metNames$`Biochemical class`)
metNames$`Biochemical class` <- gsub('Lipid','Lipids',metNames$`Biochemical class`)
metNames$`Biochemical class` <- gsub('Peptide','Peptides',metNames$`Biochemical class`)
metNames$`Biochemical class`[which(is.na(metNames$`Biochemical class`))] <- 'unnamed metabolites'

# Get all metabolites used in this analysis

# load the 80% metabolites
source('./scripts/preprocessing.R')
eightypercmets <- preprocessing(TRUE, 'metnames')

# Obtain all significant metabolites names
regRes <- read.csv('./results/CovidMetaanalysisFE_EstimatesTimUpdate0805.csv') %>% filter(QFDR>0.05 & estimateFDR<0.05)
regMets <- regRes$metabolite

# Get absolute and relative metabolite counts
metNumbers <-metNames %>%
  group_by(`Biochemical class`) %>%
  summarise(meas_ID = n(), 
            an_ID = sum(CHEMICAL_NAME  %in% eightypercmets),
            drop_ID = sum(!CHEMICAL_NAME  %in% eightypercmets),
            sig_ID = sum(CHEMICAL_NAME  %in% regMets)) %>%
  mutate(`Measured metabolites` = meas_ID/sum(meas_ID),
         `Analysed metabolites` = an_ID/sum(an_ID),
         `Dropped metabolites` = drop_ID/sum(drop_ID),
         `Significant metabolites`= sig_ID/sum(sig_ID))

# Get the absolute and relative metabolite counts
metAbs <- metNumbers[,1:5] %>% pivot_longer(!`Biochemical class`, names_to = "samples", values_to = "absCount")
metNorm <- metNumbers[,c(1,6:9)] %>% pivot_longer(!`Biochemical class`, names_to = "samples", values_to = "relCount")
metNumbers1 <- cbind(metAbs[,c(1,3)], metNorm[,2:3])

# Add column with the total metabolite numbers
metNumbers1 <- metNumbers1 %>% 
  add_column(total = rep(NA, times = dim(metNumbers1)[1])) %>%
  mutate(total = case_when(
    samples=='Measured metabolites' ~  paste0('n=',as.character(sum(metNumbers$meas_ID))),
    samples=='Analysed metabolites' ~ paste0('n=',as.character(sum(metNumbers$an_ID))),
    samples=='Dropped metabolites' ~ paste0('n=',as.character(sum(metNumbers$drop_ID))),
    samples=='Significant metabolites' ~ paste0('n=',as.character(sum(metNumbers$sig_ID)))))

metNumbers1$samples <- factor(metNumbers1$samples, levels=c('Measured metabolites','Dropped metabolites','Analysed metabolites','Significant metabolites'))

# Get colors

png(file="donutCharts.png", width=4, height=16, units="in", res=400)
s = 4
ggplot(metNumbers1, aes(x = 3, y = relCount, fill = `Biochemical class`)) +
  geom_col() +
  labs(title = 'Distributions of biochemical families') +
  geom_text(aes(label = absCount),
            position = position_stack(vjust = 0.5), colour = 'black', fontface = 'bold',size =s*1.2) +
  geom_text(aes(label = total),
            x=0,y=0, size =s*1.2) + 
  facet_wrap(~samples, ncol=1) +
  coord_polar(theta = "y") +
  xlim(c(0.2, 3 + 0.5)) +
  labs(x=NULL, y=NULL) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(8, "Spectral"))(10)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold'),
    strip.text = element_text(size = s*4, face = 'bold'),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = s*3),
    legend.position = "bottom",
    panel.grid = element_blank())

dev.off()   

################################################################################
# Figure 2
# Create a vulcano plot of all regression results
regRes <- read.csv('CovidMetaanalysisFE_Estimates.csv') 

# Get full data metabolite names and super-pathway info
metNames <- read_xlsx('APCM-02-21MD+ MERGED DATA TABLES_Correct.xlsx',
                      sheet = 'Chemical Annotation Common')

# rename super pathway in Biochemical family
metNames <- rename(metNames, `Biochemical family` = SUPER_PATHWAY)

# CONVERT `Biochemical family` names
metNames$`Biochemical family` <- gsub('Amino Acid','Amino Acids',metNames$`Biochemical family`)
metNames$`Biochemical family` <- gsub('Carbohydrate','Carbohydrates',metNames$`Biochemical family`)
metNames$`Biochemical family` <- gsub('Nucleotide','Nucleotides',metNames$`Biochemical family`)
metNames$`Biochemical family` <- gsub('Energy','Energy metabolism',metNames$`Biochemical family`)
metNames$`Biochemical family` <- gsub('Lipid','Lipids',metNames$`Biochemical family`)
metNames$`Biochemical family` <- gsub('Peptide','Peptides',metNames$`Biochemical family`)
metNames$`Biochemical family`[which(is.na(metNames$`Biochemical family`))] <- 'Unnamed Compounds'

metNames <- metNames %>% select(c(`Biochemical family`, SUB_PATHWAY, CHEMICAL_NAME))
colnames(metNames)[3] <- 'metabolite'

# Get regression results table
resTab <- merge(metNames,regRes,by.y = 'metabolite') %>%
  select(-c('CI_Geneva', 'CI_StGallen','CI_Ticino','SE','ConfInt')) %>%
  write.csv(file = 'regressionSummary.csv',row.names = FALSE)

regRes1 <- merge(metNames,regRes,by.y = 'metabolite') %>% 
  mutate(status = rep(NA,times = NROW(regRes)),
         status = if_else(estimateFDR<0.05 & QFDR>0.05 & estimateMeta>0, 'Consistent increase', status),
         status = if_else(estimateFDR<0.05 & QFDR>0.05 & estimateMeta<0, 'Consistent decrease', status),
         status = if_else(estimateFDR<0.05 & QFDR<0.05 & estimateMeta>0, 'Inconsistent increase', status),
         status = if_else(estimateFDR<0.05 & QFDR<0.05 & estimateMeta<0, 'Inconsistent decrease', status),
         status = if_else(estimateFDR>0.05, 'No change', status),
         status = factor(status, levels = c('Consistent increase',
                                            'Inconsistent increase',
                                            'Consistent decrease',
                                            'Inconsistent decrease',
                                            'No change')))
regRes1 <- regRes1 %>% mutate(
  delabel = if_else(estimateFDR<6.09e-18,regRes1$metabolite,NA),
  delabel = str_remove(delabel, "\\(.+"))

graphics.off()
png('meta_vulcanoplot.png', width = 9, height = 8, unit = 'in', res = 200)
s=20
cbPalette <- c('#E32528','#E3AEAF','#2977E6','#A5C8E8','#CDCECE')
theme_set(theme_classic(base_size = s))
regRes1 %>% 
  ggplot(aes(x=estimateMeta,y=-log10(estimateFDR),label = delabel)) +
  geom_point(size = 3,aes(colour = status)) +
  geom_text_repel(max.overlaps = Inf ,size = s*0.25) +
  scale_colour_manual(values=cbPalette) +
  
  geom_hline(aes(yintercept = -log10(0.05)), col = "black", linetype = 'dashed') +
  geom_text(aes(0.3,-log10(0.05),label = 'FDR<0.05', vjust = -0.5)) +
  labs(
    title = 'Time-dependent serum metabolite concentration changes',
    x='Estimate',
    y='-log10(p)') +
  theme(
    axis.title.y = element_text(face = "bold",  size =s*0.8, color = 'black'),
    axis.title.x = element_text(hjust = 0.5, face = "bold", size = s*0.8, color = 'black'),
    plot.title = element_text(hjust = 0, size = s, face = 'bold'),
    plot.subtitle=element_text(hjust=0.5),
    legend.position = "bottom",
    legend.title = element_blank()) +
  guides(color=guide_legend(nrow=2, byrow=TRUE))
dev.off()

################################################################################
# Visualise meta-analysis results by forest plots
# Figure 3,4,5,6,7,8,9B, and S3

# Load libraries 
library(tidyverse)
library(stringr)
library(gridExtra)
library(forestplot)

# Load meta analysis results
metaRes <- read.csv('CovidMetaanalysisFE_Estimates.csv')

#metaboliteList <- mets

# Uncomment the metabolite names to visualise them
metaboliteList <- c(
            # FIGURE 3
            # 'perfluorooctanoate (PFOA)',
            # 'propyl 4-hydroxybenzoate sulfate',
            # 'methyl-4-hydroxybenzoate sulfate',
            # '2,4-di-tert-butylphenol'
            
            # FIGURE 4
            # 'salicylate',
            # '2-methoxyacetaminophen glucuronide*',
            # '4-acetaminophen sulfate',
            # '3-(methylthio)acetaminophen sulfate'

            # FIGURE 5
            # 'S-methylcysteine sulfoxide',
            # 'carotene diol (1)',
            # 'carotene diol (3)',
            # 'erythritol',
            # 'maltol sulphate',
            # 'vanillic alcohol sulphate',
            # 'catechol sulphate',
            # '4-methylcatechol sulphate',
            # 'guaiacol sulphate',
            # '4-vinylguaiacol sulphate',
            # 'homostachydrine',
            # 'benzoate',
            # 'hippurate'
            # 'caffeine',
            # 'paraxanthine',
            # 'theobromine',
            # 'theophylline',
            # '1-methylxanthine'

            # FIGURE 6
            # 'glycodeoxycholate',
            # 'glycodeoxycholate 3-sulfate',
            # 'glycochenodeoxycholate',
            # 'glycochenodeoxycholate 3-sulfate',
            # 'glycochenodeoxycholate glucuronide (1)',
            # 'lithocholate sulfate (1)',
            # 'glycolithocholate',
            # 'glycolithocholate sulfate*',
            # 'taurolithocholate 3-sulfate',
            # 'tauroursodeoxycholate',
            # 'glycoursodeoxycholate',
            # 'glycoursodeoxycholic acid sulfate (1)',
            # 'cholate',
            # 'hyocholate',
            # 'glycohyocholate',
            # '3-formylindole',
            # 'indolelactate',
            # 'indoleacetate',
            # 'anthranilate',
            # '8-methoxykynurenate'
            
            # FIGURE 7
            # 'carboxyethyl-GABA',
            # 'fibrinopeptide A (2-15)**',
            # 'fibrinopeptide A (3-16)**',
            # 'fibrinopeptide A (4-15)**',
            # 'fibrinopeptide A (5-16)*',
            # 'fibrinopeptide A (8-16)**',
            # 'fibrinopeptide A, des-ala(1)*',
            # 'fibrinopeptide B (1-9)**',
            # 'fibrinopeptide B (1-11)**',
            # 'fibrinopeptide B (1-12)**',
            # 'fibrinopeptide B (1-13)**',
            # 'cholesterol sulfate',
            # 'bilirubin degradation product, C16H18N2O5 (1)**',
            # 'bilirubin degradation product, C16H18N2O5 (2)**',
            # 'bilirubin degradation product, C16H18N2O5 (3)**',
            # 'bilirubin degradation product, C16H18N2O5 (4)**',
            # 'bilirubin degradation product, C17H18N2O4 (1)**',
            # 'bilirubin degradation product, C17H18N2O4 (2)**',
            # 'bilirubin degradation product, C17H18N2O4 (3)**',
            # 'bilirubin degradation product, C17H20N2O5 (1)**',
            # 'bilirubin degradation product, C17H20N2O5 (2)**',
            # '3-methylglutaconate',
            # 'pyridoxal',
            # 'retinol (vitamin A)',
            # 'alpha-tocopherol',
            # 'gamma-tocopherol/beta-tocopherol',
            # '2-O-methylascorbic acid'
            
            # FIGURE 8
            # 'behenoyl sphingomyelin (d18:1/22:0)*',
            # 'palmitoyl sphingomyelin (d18:1/16:0)',
            # 'sphingomyelin (d17:2/16:0, d18:2/15:0)*',
            # 'sphingomyelin (d18:1/17:0, d17:1/18:0, d19:1/16:0)',
            # 'sphingomyelin (d18:1/19:0, d19:1/18:0)*',
            # 'sphingomyelin (d18:1/20:0, d16:1/22:0)*',
            # 'sphingomyelin (d18:1/21:0, d17:1/22:0, d16:1/23:0)*',
            # 'sphingomyelin (d18:1/22:2, d18:2/22:1, d16:1/24:2)*',
            # 'sphingomyelin (d18:1/25:0, d19:0/24:1, d20:1/23:0, d19:1/24:0)*',
            # 'sphingomyelin (d18:2/16:0, d18:1/16:1)*',
            # 'sphingomyelin (d18:2/18:1)*',
            # 'glycerophosphoserine*',
            # 'palmitoyl-sphingosine-phosphoethanolamine (d18:1/16:0)',
            # 'glycerophosphorylcholine (GPC)',
            # 'sphinganine-1-phosphate',
            # 'valine',
            # 'arginine',
            # 'lysine',
            # 'histidine',
            # 'glycine',
            # 'serine'
  
            # FIGURE 9B
            # 'aspartate',
            # 'fumarate',
            # 'arginine',
            # 'urea',
            # 'ornithine',
            # 'citrulline',
            # 'citrate',
            # 'isocitrate',
            # 'alpha-ketoglutarate',
            # 'succinate',
            # 'malate'
            
            # Figure S3
            # 'glycerophosphorylcholine (GPC)',
            # 'cysteine s-sulfate',
            # 'glycochenodeoxycholate glucuronide (1)',
            # 'pristanate',
            # 'sarcosine',
            # 'N-acetylglutamate',
            # '1-(1-enyl-stearoyl)-2-arachidonoyl-GPE (P-18:0/20:4)*',
            # 'N1-methylinosine',
            # '1-oleoyl-2-docosahexaenoyl-GPE (18:1/22:6)*',
            # 'N-acetylglucosamine/N-acetylgalactosamine',
            # '3-(methylthio)acetaminophen sulfate*',
            # 'bilirubin degradation product, C17H20N2O5 (1)**',
            # '2-methoxyhydroquinone sulfate (1)',
            # '11beta-hydroxyandrosterone glucuronide',
            # 'beta-sitosterol',
            # 'bilirubin degradation product, C17H18N2O4 (3)**',
            # 'N-acetylproline'
  )
            

# Get significant metabolites
metaboliteList <-severity$Metabolite[1:17]
            
metaboliteList <- gsub('\\*','',metaboliteList)
metaboliteList <- gsub('sulfate','sulphate',metaboliteList)

titles <- c('Meta-analysed regressions of metabolites associated with environmental exposures',
            'Meta-analysed regressions of metabolites associated with drug metabolism',
            'Meta-analysed regressions of diet related compounds',
            'Meta-analysed regressions of metabolites associated with host-gut crosstalk',
            'Meta-analysed regressions of metabolites related to physiological functioning',
            'Meta-analysed metabolite concentration changes over time in the urea and TCA cycles',
            'Meta-analysed regressions of metabolites with disease-dependent trajectories')
title <- titles[7]

# Create function that automatically creates and a forest plot 
getForestPlot <- function(metaRes, metaboliteList, title, w, h) {
# Input: metaRes, metabolite list, title
# Output: Forest plot (plt)
# Remove any stars
metaRes$metabolite <- gsub('\\*','',metaRes$metabolite)
metaRes$metabolite <- gsub('sulfate','sulphate',metaRes$metabolite)

# Get meta analysis results for metabolites in the list
metaRes1 <- metaRes %>% 
  select(metabolite, estimateMeta, estimateFDR, ConfInt) %>%
  filter(metaRes$metabolite %in% metaboliteList) %>%
  arrange(factor(metabolite, levels = metaboliteList))

# remove *
metaRes1$metabolite <- gsub('\\*','',metaRes1$metabolite)
#metaRes1$metabolite <- gsub('\\(1\\)','',metaRes1$metabolite)
#metaRes1$metabolite <- gsub('\\(3\\)','',metaRes1$metabolite)
metaRes1$metabolite <- gsub('\\(GPC\\)','',metaRes1$metabolite)
metaRes1$metabolite <- gsub('\\(P-18:0\\/20:4\\)','',metaRes1$metabolite)
metaRes1$metabolite <- gsub('\\((d18:1\\/16:0)\\)','',metaRes1$metabolite)
metaRes1$metabolite <- gsub('\\(18:1\\/22:6\\)','',metaRes1$metabolite)
metaRes1$metabolite <- gsub('sulfate','sulphate',metaRes1$metabolite)

# Get lower and upper confidence intervals and round values
res <- 3
metaRes1$CI_lower <- str_extract(metaRes1$ConfInt,'\\((.*)\\,')
metaRes1$CI_lower <- as.numeric(substr(metaRes1$CI_lower,2,nchar(metaRes1$CI_lower)-1))
metaRes1$CI_upper <- str_extract(metaRes1$ConfInt,'\\,(.*)\\)')
metaRes1$CI_upper <- as.numeric(substr(metaRes1$CI_upper,2,nchar(metaRes1$CI_upper)-1))

# Create new confint variable
roundCI <- (\(x) as.character(round(x,digits = 3)))
metaRes1$ConfInt <- paste0('[',roundCI(metaRes1$CI_lower),',',roundCI(metaRes1$CI_upper),']')

# Factorise metabolite names
colnames(metaRes1) <- c('metabolite','mean','FDR','CI','lower','upper')

# VISUALISE
plt <-
  metaRes1 |> 
  mutate(est = sprintf("%.3f", mean), .after = metabolite) |> 
  forestplot(fn.ci_norm = local({
    i = 0
    b_clrs = rep('black',NROW(metaRes1))
    b_clrs[metaRes1$FDR<0.05]<-'#d32f2f'
    function(..., clr.line, clr.marker){
      i <<- i + 1
      fpDrawCircleCI(...,  clr.marker = b_clrs[i])
    }
  }),
  labeltext = c(metabolite), 
             title = title,
             xlab = "Estimate", 
             hrzl_lines = TRUE,
             ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             xticks.digits = 4,
             lwd.zero = 1.5,
             grid = TRUE,
             boxsize = 0.5,
             xticks = c(-0.3, -0.2,-0.1, 0, 0.1, 0.2, 0.3),
             txt_gp = fpTxtGp(cex=1.5,ticks=gpar(cex=1.5), xlab = gpar(cex=1.5)),
             axes = gpar(cex = 1.5)) |>
  fp_add_header(metabolite = 'Metabolite') |> 
  fp_set_style(zero = "gray50") #|> 
  #fp_add_lines(h_2 = gpar(col = "black", columns = 1:2, lty = 1),
               #h_6 = gpar(col = "black", columns = 1:2, lty = 1),
               #h_11 = gpar(col = "black", columns = 1:2, lty = 1),
               #h_24 = gpar(col = "black", columns = 1:2, lty = 1))
               # h_44 = gpar(col = "black", columns = 1:2, lty = 1),
               # h_71 = gpar(col = "black", columns = 1:2, lty = 1))

# Save plot
png(file=str_c(title,'FP.png'), width=w, height=h, units="in", res=300)
print(plt)
dev.off()

return(plt)
}

# run function
plt <- getForestPlot(metaRes, metaboliteList, title, 8*2, 3+length(metaboliteList)*0.5)
print(plt)
################################################################################
# Create Figure 10

# Load libraries
library(tidyverse)
library(readxl)
library(ggrepel)

# Load results on moderate/severe covid cases
severity <- read_excel('Supplementary tables.xlsx', sheet = 'Table S7')

# Get significant metabolites
mets <-severity$Metabolite[1:17]

# Load real data
path <-r'(C:\Users\mspg\National University of Ireland, Galway\Group_MSP - COVID-19 patient data\PROCESSED_COVID19_METABOLOMICS_Common_SAMPLES_IMPUTED.xlsx)'
data <- read_excel(path) %>% filter(Location == 'Ticino')

# Remove patients without a second timepoint
data$Patient_ID <- sub('_[0-9]', '', data$Patient_ID)
data <- data[data$Patient_ID %in% data$Patient_ID[duplicated(data$Patient_ID)],]

# Get sample timepoints
carn <- which(colnames(data)=='carnitine')
sample =  sub('.*(?=.$)', '', data$COVID_SEVERITY, perl=T)
sample[grep('[^0-9]', sample)] <- "1"
data$Timepoint <- as.numeric(sample)
data <- data %>% relocate(Timepoint, .after = colnames(data)[carn-1])

# remove >sample 2
data <- data[data$Timepoint<=2,]
data$Timepoint <- as.factor(data$Timepoint)

# Remove _TX in COVID_SEVERITY
data$COVID_SEVERITY <- gsub('_T[0-9]','',data$COVID_SEVERITY)
data$COVID_SEVERITY <- gsub('COVID_','',data$COVID_SEVERITY)

# Log transform and scale value
data<- data %>% mutate(across(carnitine:`X-26119`, log)) %>% mutate(across(carnitine:`X-26119`, scale))

# Select metabolites
data <- data %>% select(c('Patient_ID', 'COVID_SEVERITY', 'Timepoint', names(data)[which(names(data) %in% mets)]))  %>%# Filter on 18 metabolites
  pivot_longer(cols = -c('Patient_ID','COVID_SEVERITY','Timepoint'), names_to = 'metabolite',values_to = 'value') # sort metabolites on significance

# Tweak metabolite names
# remove *
data$metabolite <- gsub('\\*','',data$metabolite)
data$metabolite <- gsub('\\(1\\)','',data$metabolite)
data$metabolite <- gsub('\\(3\\)','',data$metabolite)
data$metabolite <- gsub('\\(GPC\\)','',data$metabolite)
data$metabolite <- gsub('\\(P-18:0\\/20:4\\)','',data$metabolite)
data$metabolite <- gsub('\\((d18:1\\/16:0)\\)','',data$metabolite)
data$metabolite <- gsub('\\(18:1\\/22:6\\)','',data$metabolite)
data$metabolite <- gsub('sulfate','sulphate',data$metabolite)

# Create grouped boxplot
png('diseaseSevDependentMets.png',width=9, height = 12, unit = 'in', res = 200)
data %>% 
  ggplot(aes(x=COVID_SEVERITY, y = value, fill = Timepoint)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  facet_wrap(~metabolite,nrow =6) +
  labs(title = 'Serum metabolites with disease-dependent trajectories in Ticino', 
       x= '',
       y='Logscaled and z-transformed serum concentration') +
  theme(axis.text.x = element_text(face='bold'))
dev.off()

################################################################################
# This script creates Figure S3

# load libraries
library(ggsci)
library(readxl)
library(tidyverse)

# load metabolomics dataset
data <- read_excel("PROCESSED_COVID19_METABOLOMICS_Common_SAMPLES_RAW_NO_IMP.xlsx")

# filter data
data <- data %>%
  # Remove samples from Cork
  filter(Location != 'Cork') %>%
  # remove on post-covid and control samples
  slice(which(!grepl('Long|control', COVID_SEVERITY)))

# Get first metabolite column
carn <- which(colnames(data)=='carnitine')

# Add sample numbers
sample =  sub('.*(?=.$)', '', data$COVID_SEVERITY, perl=T)
sample[grep('[^0-9]', sample)] <- "1"
data$sample <- as_factor(sample)
data <- data %>% relocate(sample, .after = colnames(data)[carn-1])

# Remove patient numbers
data$Patient_ID <- sub('_[0-9]', '', data$Patient_ID)
data$COVID_SEVERITY <- sub('_T[0-9]', '', data$COVID_SEVERITY)

# Remove geneva patient without metadata
data<-data %>% filter(! Patient_ID %in% 'COV-140')

# Get the number of patients per location
data %>%
  group_by(Location) %>%
  summarise(total = length(unique(Patient_ID)))

# Get the number of samples per location
data %>%
  group_by(Location) %>%
  summarise(total = n())

# Get the number of samples restricted to 9 days for each location per severity
before <-
  data %>%
  group_by(Location, COVID_SEVERITY, `Death due to COVID-19`,sample) %>%
  summarise(total = n())

after <- data %>% 
  filter(`Days Post Hospitalisation`<9) %>%
  group_by(Patient_ID) %>%
  mutate(sample = as.numeric(sample)) %>%
  filter(max(sample)>1) %>%
  ungroup() %>%
  group_by(Location, COVID_SEVERITY, `Death due to COVID-19`,sample) %>%
  summarise(analysed = n())

# Get the number of samples per location
after %>%
  group_by(Location) %>%
  summarise(total = sum(analysed))

# Process data further to create a visualisations of the sample distributions
data1 <- 
  merge(before, after, by = c('Location', 'COVID_SEVERITY', 'Death due to COVID-19', 'sample'), all = TRUE) %>%
  pivot_longer(cols = total:analysed, names_to = 'type',values_to = 'value') %>%
  mutate(type = factor(type,levels = c('total','analysed')), Timepoint = as_factor(sample))

# Visualise data
library(ggplot2)
png('sampleDistributions.png', width = 12, height = 7, unit = 'in', res = 200)
s=14
ggplot(data1) + 
  geom_col(aes(y=value, x=type,fill=Timepoint)) + 
  facet_grid(.~Location) +
  labs(title = 'Sample distribution per location',  x = "", y = 'Number of samples', legend = 'Timepoint') +
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw() + 
  theme(
    plot.title = element_text(size=s*1.2, face = 'bold'),
    strip.text = element_text(size=s, face = 'bold'),
    axis.text.x = element_text(size=s),
    axis.title.y = element_text(size=s),
    axis.title.x = element_text(size=s),
    legend.text = element_text(size=s),
    legend.title = element_text(size=s))
dev.off()