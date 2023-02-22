# Education, Epigenetic Age Acceleration, and Dementia 

## Title: Additive interaction and mediation-interaction decomposition: DNA methylation age acceleration, education, and cognitive impairment in the Health and Retirement Study

## Citation Information
Ware EB, Higgins C, Harris S, Fisher JD, Bakulski KM. Additive interaction and mediation-interaction decomposition: DNA methylation age acceleration, education, and cognitive impairment in the Health and Retirement Study. In review.

This Github repository contains the data management and analytic scripts to produce the following manuscript:[DNA methylation age acceleration, education, and cognitive impairment in the Health and Retirement Study](https://www.medrxiv.org/content/10.1101/2022.12.23.22283919v1.full-text)

## Abstract

### Importance: 
Dementia represents a significant and increasing public health burden. DNA methylation
age acceleration may be associated with dementia and dementia risk factors, such as education, but
investigating their impact on dementia is necessary.
Objective: To evaluate the association of educational attainment on dementia and cognitive
impairment through DNA methylation age acceleration, while accommodating exposure-mediator
interaction effects.
### Design: 
In the 2016 Health and Retirement Study wave, we evaluated six epigenetic clocks, including
GrimAge, with Langa-Weir classified dementia, cognitive impairment non-dementia, and normal
cognition. Age acceleration was the residual between chronological age of participants and DNA
methylation clock, dichotomized at zero. To understand the joint association of low education (≤12
years) and DNA methylation age acceleration in relation to cognitive impairment, we used weighted
logistic regression and calculated interaction on the additive scale adjusting for chronological age, sex,
race/ethnicity, and cell type composition. We performed four-way mediation and interaction
decomposition analysis to estimate the: 1) controlled direct effect of education on cognition, 2)
interaction reference, 3) interaction mediation, and 4) the pure indirect effect of DNA methylation age
acceleration on cognition.
### Setting: 
Analysis was conducted on a subsample of Health and Retirement Study participants in the
2016 Venous Blood Study (N=3,724).
### Results: 
Both GrimAge acceleration (OR=1.6 95%CI 1.3 – 2.1) and low educational attainment (OR=
2.4 95%CI 1.9 – 3.0) were associated with higher odds of cognitive impairment, non-dementia in a
mutually adjusted logistic model. We found additive interaction associations between low education
and GrimAge acceleration on dementia. We observed that 6-8% of the association of education on
dementia was mediated through GrimAge acceleration. While mediation effects were small, the portion
of the association of education due to additive interaction with GrimAge acceleration was between 23.6
and 29.2%.
### Conclusions and Relevance: 
Accelerated DNA methylation age was associated with increased odds of
cognitive impairment and we observed more than additive interaction effects between education and
age acceleration on dementia. These results support the interplay of social disadvantage and biological
aging processes on impaired cognition.

## Funding
This work was supported by grants from the National Institute on Aging (R01 AG055406, R01 AG067592, P30 AG072931, R01 AG055654, R01 AG067592-01S1, R01AG06759201, R25AG053227); and the National Institute of Minority Health and Disparities (R01 MD013299); and the National Institute of Environmental Health Sciences (P42ES017198, P30ES017885); and the National Center for Advancing Translational Sciences (UL1TR002240).


## Content

This repo contains the STATA script used for data management, and data analysis of the manuscript. 
We used STATA for our main models and analysis of measures of interaction in the additive scale, and mediation-interaction decomposition analysis. 
We built our own code in STATA to calculate measures of interaction in the additive scale and the mediation-interaction decomposition analysis. 
We decided not rely on existing packages (i.e., STATA med4way: https://pubmed.ncbi.nlm.nih.gov/30452641/) 
because this package does not accommodate for survey weights. 

Health and Retirement Study survey data are publicly available (https://hrs.isr.umich.edu/data-products), 
and genetic data are available through dbGaP (https://dbgap.ncbi.nlm.nih.gov; phs000428.v2.p2)
Original Article

A preprint version of this article can be found at medRxiv (https://www.medrxiv.org/content/10.1101/2022.12.23.22283919v1)

## Script files
Export File Tables: including the tables related to this project

Stata Script: including the STATA .do file which contains code on creating epigenetic clock age acceleration 
variables and tables



