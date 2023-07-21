# CovidMetabolomicsAnalysis

### Readme for scripts used in CovidMetabolomicsAnalysis ###
### Tim Hensen, 07/2023 ###

This folder contains the scripts needed to execute the data analyses performed in Hensen et al., “The effects of hospitalisation on the serum metabolome in COVID-19 patients” 

Before running the analyses, the raw data was prepared in “metabolomicsDataPreparation.mlx”.

Time-dependent linear mixed-effect modelling and subsequent meta-analyses are performed in “Regressions&MetaAnalysis.R”. All visualisations in the manuscript, with exceptions for Figure 11 and Figure S2, are obtained by running “Visualisations.R”.

Analyses on COVID-19 severity and mortality and analyses on metabolite-pair bivariate distributions are obtained by running the STATA scripts: “analyses_Gallen_Geneva_death.do”, “analyses_Ticino_moderate_severe.do”, and “combined_analyses_longitudinal.do”. Visualisations for Figure 11 are created by running “analyses_Ticino_moderate_severe.do” and “combined_analyses_longitudinal.do”.

“MetaboliteDepletionAnalyses.R” runs metabolite depletion regressions and subsequent meta-analyses which not included in the manuscript results, but are mentioned in the methods.

The datasets used in this study are available upon reasonable request. The data is not publicly available due to privacy.
