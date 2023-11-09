library(TwoSampleMR)
library(data.table)
library("readxl")

# Instruments of Iceland
Ins <- read_excel("Instruments_Iceland.xlsx",1)

# exposure
Ins<-format_data(Ins, type = "exposure", header = TRUE,
                 phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                 se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                 other_allele_col = "other_allele", pval_col = "pval")
# outcome
## Format the outcome data
perform_format <- function(root_data){
  outcome <- as.data.frame(fread(root_data))
  mycols <- c("chromosome", "base_pair_location", "other_allele", "effect_allele", 
              "effect_allele_frequency", "beta", "standard_error", "p_value")
  gwas <- outcome[, mycols]
  colnames(gwas) <- c("CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P")
  myfile <- format_sumstats(gwas, nThread = 12, ref_genome = "GRCh37", dbSNP=144)
  gwas <- as.data.frame(fread(myfile)) # please note A2 is effect allele
  return (gwas)
  }
## pancreatic iron content
PIC_gwas <- perform_format("GCST90016676_buildGRCh37.tsv.gz")
PIC_gwas <- format_data(PIC_gwas, "outcome", snp_col = "SNP", beta_col = "BETA", se_col = "SE",
                    eaf_col = "FRQ", effect_allele_col = "A2", other_allele_col = "A1", pval_col = "P")
## pancreatic cancer
# FinnGen database
pan_gwas <- as.data.frame(data.table::fread("D:/2023/Tasks/pQTL/0_pancreatic_cancer_data/finngen_R8_C3_PANCREAS_EXALLC.gz"))
pan_gwas <- format_data(pan_gwas, "outcome", snp_col = "rsids", beta_col = "beta", se_col = "sebeta",
                    eaf_col = "af_alt", effect_allele_col = "alt", other_allele_col = "alt", pval_col = "pval")
# Harmonise
dat <- NULL
dat <- harmonise_data(
  exposure_dat = Ins, 
  outcome_dat = PIC_gwas
)

##run the MR and sensitivity analyses 
mr_results <- NULL
mr_hetero <- NULL
mr_pleio <- NULL
mr_single <- NULL
try(mr_results <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw")))  # main MR analysis
mr_hetero <- mr_heterogeneity(dat) # heterogeneity test across instruments
mr_pleio <- mr_pleiotropy_test(dat) # MR-Egger intercept test  
try(mr_single <- mr_singlesnp(dat)) #single SNP MR using Wald ratio

fwrite(dat, paste0("harmonise.csv"))
fwrite(mr_results, paste0("mr.csv"))
fwrite(mr_hetero, paste0("mr_hetero.csv"))
fwrite(mr_pleio, paste0("mr_pleio.csv"))
fwrite(mr_single, paste0("mr_single.csv"))