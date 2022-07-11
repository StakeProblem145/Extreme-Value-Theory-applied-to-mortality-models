# Bachelor Thesis
## Extreme Value Theory applied to mortality models

### Data
The HMD is against copies of their data in public repositories thus the user has to manually download the dataset themself.

The HMD data for Canada can be downloaded here:
https://mortality.org/Country/Country?cntr=CAN

The Deaths and Exposure-to-risk are necessary.
Afterwards the script found in the folder Scripts can be used to generate the dataset for further use. The output has to be placed in the Data folder.


The Swedish data may not be disclosed.


### Program Structure
graduatePoisson.R is the core of the thesis. Here the functions to fit via Maximum Likelihood Estimation are implemented.  
plottingFunctions.R  implements functions to plot fits, residuals and parameters.  
helperFunctions.R small functions which have been used for the analysis.

AnalysisAggregated.R Analysis of the aggreagated datasets.  
AnalysisCohort.R Analysis of the individual cohorts.  
AnalysisCohortParameters.R Analysis of the cohort parameters.  
