##Clean up R

rm(list=ls(all=TRUE))

####################################
##Install  packages
install.packages("devtools")
library(devtools)
install_github("jrs95/nlmr")
library(nlmr)
install_github("WSpiller/RadialMR")
library(devtools)
install_github("qingyuanzhao/mr.raps")

install.packages(c("ggplot2", "ggrepel", "igraph"))

library(devtools)

devtools::install_github('MRCIEU/TwoSampleMR')

library(TwoSampleMR)
library(RadialMR)
library(mr.raps)
library(ggplot2)

##

getwd()

setwd("...data/")
###################################

#Variant level analysis

###################################



##First, read in with the non-standard headings

bmi_exp_dat<-read_exposure_data(
  filename="....csv",
  sep = ",",
  snp_col = "snp",
  beta_col = "betaexposure",
  se_col = "seexposure",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "pvalexposure",
  units_col = "unitsexposure",
  #gene_col = "Gene",
  samplesize_col = "samplesizeexposure"
)



##Clump this file - 

bmi_exp_dat <- clump_data(bmi_exp_dat)



cost_out_dat<-read_outcome_data(
  filename="...betas_bmi_snp_pheno_genetic.csv",
  snps = bmi_exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta_outcome",
  se_col = "se_outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p",
  units_col = "unit_output",
  #gene_col = "Gene",
  samplesize_col = "samplesize_out"
)


##Harmonise data 

dat <- harmonise_data(
  exposure_dat = bmi_exp_dat, 
  outcome_dat = cost_out_dat
)





##IDentify methods

mr_method_list()


res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))
resegger <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
resmedian <- mr(dat, method_list=c("mr_ivw","mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median"))	
resmode <- mr(dat, method_list=c("mr_ivw", "mr_simple_mode", "mr_weighted_mode"))			
resraps <- mr(dat, method_list = c("mr_ivw","mr_raps"), parameters = list(over.dispersion = FALSE, loss.function = "l2"))
resivw<-mr(dat, method_list=c("mr_ivw"))


######################## Main results ########################

resegger
resmedian
resmode
resraps


#Standardise using sd=4.6 kg/m2 following Budu-Aggrey et al, https://www.biorxiv.org/content/10.1101/265629v1.abstract

standard.resegger<- data.frame("outcome"=resegger$outcome, "exposure"=resegger$exposure, "method"=resegger$method, "nsnp"=resegger$nsnp, "b"=resegger$b/4.6, "se"=resegger$se/4.6, "pval"=resegger$pval)      
standard.resmedian<- data.frame("outcome"=resmedian$outcome, "exposure"=resmedian$exposure, "method"=resmedian$method, "nsnp"=resmedian$nsnp, "b"=resmedian$b/4.6, "se"=resmedian$se/4.6, "pval"=resmedian$pval)      
standard.resmode<- data.frame("outcome"=resmode$outcome, "exposure"=resmode$exposure, "method"=resmode$method, "nsnp"=resmode$nsnp, "b"=resmode$b/4.6, "se"=resmode$se/4.6, "pval"=resmode$pval)      
standard.resraps<- data.frame("outcome"=resraps$outcome, "exposure"=resraps$exposure, "method"=resraps$method, "nsnp"=resraps$nsnp, "b"=resraps$b/4.6, "se"=resraps$se/4.6, "pval"=resraps$pval)      
standard.res.summary<- data.frame("outcome"=res.summary$outcome, "exposure"=res.summary$exposure, "method"=res.summary$method, "nsnp"=res.summary$nsnp, "b"=res.summary$b/4.6, "se"=res.summary$se/4.6, "pval"=res.summary$pval)      
                                      
standard.res.summary                
                  
                  


##Steiger test for directionality - assesses causal direction ie is BMI causing cost or is cost causing BMI 

out <- directionality_test(dat)


# Heterogeneity and Egger intercept
mr_het <- mr_heterogeneity(dat)
mr_het
mr_egger_int <- mr_pleiotropy_test(dat)
mr_egger_int



# single SNP analyses
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw", "mr_egger_regression", "mr_penalised_weighted_median", "mr_weighted_mode"))
res_single

# leave one out analyses
res_loo <- mr_leaveoneout(dat)
res_loo

#######GRAPHICS 

##Plain scatter 


#Scatter plot

main.scatter<- mr_scatter_plot(res.summary, dat)
main.scatter
ggsave(main.scatter[[1]], file="main_scatter.png", width=7, height=7)

ivw.scatter.only<-mr_scatter_plot(resivw,dat)
ivw.scatter.only
ggsave(ivw.scatter.only[[1]], file="ivw_scatter.png", width=7, height=7)

#FOREST PLOT 


res_single <- mr_singlesnp(dat)
forest.plot <- mr_forest_plot(res_single)
forest.plot[[1]]
ggsave(forest.plot[[1]], file="forest.png", width=7, height=7)


##Funnel plot

res_single <- mr_singlesnp(dat)
funnel.plot <- mr_funnel_plot(res_single)
funnel.plot[[1]]
ggsave(loo.plot[[1]], file="funnel.png", width=7, height=7)

################## OTHER SENSITIVTY ANALYSES  ####################################################

##Multivariable MR 

id_exposure <- c(2, 999)
exposure_dat <- mv_extract_exposures(id_exposure)

cost_out_mv_dat<-read_outcome_data(
  filename="...betas_bmi_whr_bodyfat_snp_pheno_genetic.csv",
  snps = exposure_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta_outcome",
  se_col = "se_outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p",
  units_col = "unit_output",
  #gene_col = "Gene",
  samplesize_col = "samplesize_out"
)



mvdat <- mv_harmonise_data(exposure_dat, cost_out_mv_dat)
resmv <- mv_multiple(mvdat)
data.resmv<-data.frame(resmv)
standard.resmv<- data.frame("exposure"=data.resmv$result.exposure, "nsnp"=data.resmv$result.nsnp, "b"=data.resmv$result.b/4.6, "se"=data.resmv$result.se/4.6,  "pval"=data.resmv$result.pval)      
standard.resmv

##Non-LINEAR MOdels


library(nlmr)

setwd("...data/")

nldata<- read.csv("grs_bmi_ipd_obs_pca_complete.csv",header=T)



attach(nldata)



#Covariates on sex and first ten principal components

covar=data.frame(c1=sex,c2=pca1, c3=pca2, c4=pca3, c5=pca4, c6=pca7, c7=pca6, c8=pca7, c9=pca8, c10=pca9, c11=pca10)

genes=bmi_grs

### Analyses


fp100 = fracpoly_mr(cost_person_year,bmi_combined, genes, covar, family="gaussian", q=100, d=1, ci="model_se", fig=T)
summary(fp100)





################################################################################################################

#ANALYSING SUB COMPONENTS OF COST
#ELECTIVE/NON-ELECTIVE/OTHER 

################################################################################################################

##ELECTIVE

library(TwoSampleMR)
library(RadialMR)
library(mr.raps)
library(ggplot2)

##

getwd()

setwd("...data/")
###################################

#Variant level analysis

###################################


bmi_exp_dat<-read_exposure_data(
  filename="...betas_bmi_snp_pheno_genetic.csv",
  sep = ",",
  snp_col = "snp",
  beta_col = "betaexposure",
  se_col = "seexposure",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "pvalexposure",
  units_col = "unitsexposure",
  #gene_col = "Gene",
  samplesize_col = "samplesizeexposure"
)



##Clump this file  

bmi_exp_dat <- clump_data(bmi_exp_dat)

cost_out_dat<-read_outcome_data(
  filename="...betas_bmi_snp_pheno_genetic.csv",
  snps = bmi_exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta_outcome_elective",
  se_col = "se_outcome_elective",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p_elective",
  units_col = "unit_output",
  #gene_col = "Gene",
  samplesize_col = "samplesize_out"
)


##Harmonise data 

dat <- harmonise_data(
  exposure_dat = bmi_exp_dat, 
  outcome_dat = cost_out_dat
)

##IDentify methods

mr_method_list()


res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))
resegger <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
resmedian <- mr(dat, method_list=c("mr_ivw","mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median"))	
resmode <- mr(dat, method_list=c("mr_ivw", "mr_simple_mode", "mr_weighted_mode"))			
resraps <- mr(dat, method_list = c("mr_ivw","mr_raps"), parameters = list(over.dispersion = FALSE, loss.function = "l2"))
resivw<-mr(dat, method_list=c("mr_ivw"))

resegger
resmedian
resmode
resraps


#Use sd=4.6 kg/m2 as above

standard.resegger<- data.frame("outcome"=resegger$outcome, "exposure"=resegger$exposure, "method"=resegger$method, "nsnp"=resegger$nsnp, "b"=resegger$b/4.6, "se"=resegger$se/4.6, "pval"=resegger$pval)      
standard.resmedian<- data.frame("outcome"=resmedian$outcome, "exposure"=resmedian$exposure, "method"=resmedian$method, "nsnp"=resmedian$nsnp, "b"=resmedian$b/4.6, "se"=resmedian$se/4.6, "pval"=resmedian$pval)      
standard.resmode<- data.frame("outcome"=resmode$outcome, "exposure"=resmode$exposure, "method"=resmode$method, "nsnp"=resmode$nsnp, "b"=resmode$b/4.6, "se"=resmode$se/4.6, "pval"=resmode$pval)      
standard.resraps<- data.frame("outcome"=resraps$outcome, "exposure"=resraps$exposure, "method"=resraps$method, "nsnp"=resraps$nsnp, "b"=resraps$b/4.6, "se"=resraps$se/4.6, "pval"=resraps$pval)      
standard.res.summary<- data.frame("outcome"=res.summary$outcome, "exposure"=res.summary$exposure, "method"=res.summary$method, "nsnp"=res.summary$nsnp, "b"=res.summary$b/4.6, "se"=res.summary$se/4.6, "pval"=res.summary$pval)      

standard.res.summary                




# Heterogeneity and Egger intercept
mr_het <- mr_heterogeneity(dat)
mr_het
mr_egger_int <- mr_pleiotropy_test(dat)
mr_egger_int



###############
#NON-ELECTIVE


bmi_exp_dat<-read_exposure_data(
  filename="...betas_bmi_snp_pheno_genetic.csv",
  sep = ",",
  snp_col = "snp",
  beta_col = "betaexposure",
  se_col = "seexposure",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "pvalexposure",
  units_col = "unitsexposure",
  #gene_col = "Gene",
  samplesize_col = "samplesizeexposure"
)



##Clump this file - 

bmi_exp_dat <- clump_data(bmi_exp_dat)

cost_out_dat<-read_outcome_data(
  filename="...betas_bmi_snp_pheno_genetic.csv",
  snps = bmi_exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta_outcome_nonelective",
  se_col = "se_outcome_nonelective",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p_nonelective",
  units_col = "unit_output",
  #gene_col = "Gene",
  samplesize_col = "samplesize_out"
)


##Harmonise data 

dat <- harmonise_data(
  exposure_dat = bmi_exp_dat, 
  outcome_dat = cost_out_dat
)





##IDentify methods

mr_method_list()

##cOMMENT 30 oCTOBER 2018 - THERE IS NO NEED TO SEPARATELY SPECIFY MR_IVw_RE AS THIS IS DEFAULT MR_IVW

res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))
resegger <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
resmedian <- mr(dat, method_list=c("mr_ivw","mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median"))	
resmode <- mr(dat, method_list=c("mr_ivw", "mr_simple_mode", "mr_weighted_mode"))			
resraps <- mr(dat, method_list = c("mr_ivw","mr_raps"), parameters = list(over.dispersion = FALSE, loss.function = "l2"))
resivw<-mr(dat, method_list=c("mr_ivw"))

resegger
resmedian
resmode
resraps

########### Express two-sample MR results in terms of units of BMI rather than SD ########################


#Use sd=4.6 kg/m2 as above

standard.resegger<- data.frame("outcome"=resegger$outcome, "exposure"=resegger$exposure, "method"=resegger$method, "nsnp"=resegger$nsnp, "b"=resegger$b/4.6, "se"=resegger$se/4.6, "pval"=resegger$pval)      
standard.resmedian<- data.frame("outcome"=resmedian$outcome, "exposure"=resmedian$exposure, "method"=resmedian$method, "nsnp"=resmedian$nsnp, "b"=resmedian$b/4.6, "se"=resmedian$se/4.6, "pval"=resmedian$pval)      
standard.resmode<- data.frame("outcome"=resmode$outcome, "exposure"=resmode$exposure, "method"=resmode$method, "nsnp"=resmode$nsnp, "b"=resmode$b/4.6, "se"=resmode$se/4.6, "pval"=resmode$pval)      
standard.resraps<- data.frame("outcome"=resraps$outcome, "exposure"=resraps$exposure, "method"=resraps$method, "nsnp"=resraps$nsnp, "b"=resraps$b/4.6, "se"=resraps$se/4.6, "pval"=resraps$pval)      
standard.res.summary<- data.frame("outcome"=res.summary$outcome, "exposure"=res.summary$exposure, "method"=res.summary$method, "nsnp"=res.summary$nsnp, "b"=res.summary$b/4.6, "se"=res.summary$se/4.6, "pval"=res.summary$pval)      

standard.res.summary                




# Heterogeneity and Egger intercept
mr_het <- mr_heterogeneity(dat)
mr_het
mr_egger_int <- mr_pleiotropy_test(dat)
mr_egger_int



###############
#OTHER


bmi_exp_dat<-read_exposure_data(
  filename="...betas_bmi_snp_pheno_genetic.csv",
  sep = ",",
  snp_col = "snp",
  beta_col = "betaexposure",
  se_col = "seexposure",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "pvalexposure",
  units_col = "unitsexposure",
  #gene_col = "Gene",
  samplesize_col = "samplesizeexposure"
)



##Clump this file - 

bmi_exp_dat <- clump_data(bmi_exp_dat)



cost_out_dat<-read_outcome_data(
  filename="...betas_bmi_snp_pheno_genetic.csv",
  snps = bmi_exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta_outcome_other",
  se_col = "se_outcome_other",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p_other",
  units_col = "unit_output",
  #gene_col = "Gene",
  samplesize_col = "samplesize_out"
)


##Harmonise data .

dat <- harmonise_data(
  exposure_dat = bmi_exp_dat, 
  outcome_dat = cost_out_dat
)





##IDentify methods

mr_method_list()


res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))
resegger <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
resmedian <- mr(dat, method_list=c("mr_ivw","mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median"))	
resmode <- mr(dat, method_list=c("mr_ivw", "mr_simple_mode", "mr_weighted_mode"))			
resraps <- mr(dat, method_list = c("mr_ivw","mr_raps"), parameters = list(over.dispersion = FALSE, loss.function = "l2"))
resivw<-mr(dat, method_list=c("mr_ivw"))

resegger
resmedian
resmode
resraps



#Use sd=4.6 kg/m2 as above

standard.resegger<- data.frame("outcome"=resegger$outcome, "exposure"=resegger$exposure, "method"=resegger$method, "nsnp"=resegger$nsnp, "b"=resegger$b/4.6, "se"=resegger$se/4.6, "pval"=resegger$pval)      
standard.resmedian<- data.frame("outcome"=resmedian$outcome, "exposure"=resmedian$exposure, "method"=resmedian$method, "nsnp"=resmedian$nsnp, "b"=resmedian$b/4.6, "se"=resmedian$se/4.6, "pval"=resmedian$pval)      
standard.resmode<- data.frame("outcome"=resmode$outcome, "exposure"=resmode$exposure, "method"=resmode$method, "nsnp"=resmode$nsnp, "b"=resmode$b/4.6, "se"=resmode$se/4.6, "pval"=resmode$pval)      
standard.resraps<- data.frame("outcome"=resraps$outcome, "exposure"=resraps$exposure, "method"=resraps$method, "nsnp"=resraps$nsnp, "b"=resraps$b/4.6, "se"=resraps$se/4.6, "pval"=resraps$pval)      
standard.res.summary<- data.frame("outcome"=res.summary$outcome, "exposure"=res.summary$exposure, "method"=res.summary$method, "nsnp"=res.summary$nsnp, "b"=res.summary$b/4.6, "se"=res.summary$se/4.6, "pval"=res.summary$pval)      

standard.res.summary                




# Heterogeneity and Egger intercept
mr_het <- mr_heterogeneity(dat)
mr_het
mr_egger_int <- mr_pleiotropy_test(dat)
mr_egger_int
