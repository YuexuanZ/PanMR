library(TwoSampleMR);library(data.table)
library(ieugwasr)
library(plinkbinr)
library("readxl")
mainMR_risk_factors <- function(exp_dat, out_dat, expName, outName, samplesize.outcome) {
  if(!is.null(out_dat)){
    mydata <- harmonise_data(exp_dat, out_dat)
    fwrite(mydata, paste0(expName, "-", outName, "_IV.csv"))
    
    if (is.na(samplesize.outcome)) {
      mydata$samplesize.outcome <- 1000
    }else{mydata$samplesize.outcome <- samplesize.outcome}
    if(dim(mydata[which(mydata$mr_keep==TRUE),])[1]>=1) {
      mydata$pval.outcome[which(mydata$pval.outcome == 0)] <- 1e-200
      res <- mr(mydata, method_list = c("mr_wald_ratio", "mr_ivw_mre", "mr_ivw_fe", 
                                        "mr_weighted_median", "mr_egger_regression", 
                                        "mr_weighted_mode"
      )
      )
      steiger <- directionality_test(mydata)
      res$exposure_r2 <- steiger[1,5]
      res$outcome_r2 <- steiger[1,6]
      res$direction <- steiger[1,7]
      res$PRESSO_beta <- NA
      res$PRESSO_se <- NA
      res$PRESSO_pval <- NA
      res$IVW_Qpval <- NA
      res$egger_pleiotropy_pval <- NA
      if(dim(mydata)[1] > 3) {
        presso <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",SdExposure = 'se.exposure',
                                      SdOutcome = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                                      data = mydata, NbDistribution = 1000,  SignifThreshold = 0.05)
        res$PRESSO_beta <- presso$`Main MR results`$`Causal Estimate`[2]
        res$PRESSO_se <- presso$`Main MR results`$`Sd`[2]
        res$PRESSO_pval <- presso$`Main MR results`$`P-value`[2]
        het <- mr_heterogeneity(mydata)
        pleio <- mr_pleiotropy_test(mydata)
        res$IVW_Qpval <- het$Q_pval[2]
        res$egger_pleiotropy_pval <- pleio$pval
      }
      res$expName <- expName
      res$outName <- outName
      
      return(res)
    }
  }
}

# PRSS1 (the same for AMY2A and CELA2A)
# GWAS summary data is available at https://www.jianguoyun.com/p/DZPUMToQ2IWVDBiHhqkFIAA
# the file names are Instruments_PRSS1.txt.gz, Instruments_AMY2A.txt.gz and Instruments_CELA2A.txt.gz 
pqtl_data <- as.data.frame(data.table::fread("Instruments_PRSS1.txt.gz"))
anno_data <- as.data.frame(data.table::fread("assocvariants.annotated.txt.gz"))
pqtl_data$eaf <- anno_data$effectAlleleFreq[match(pqtl_data$Name, anno_data$Name)]

PRSS1 <- TwoSampleMR::format_data(pqtl_data, "exposure", snp_col = "rsids", beta_col = "Beta", se_col = "SE",
                                  effect_allele_col = "effectAllele", other_allele_col = "otherAllele", 
                                  eaf_col = "eaf", pval_col = "Pval", samplesize_col = 35559
)
PRSS1 <- PRSS1[which(PRSS1$pval.exposure < 5e-8), ]

# LD
exp_dat <- ld_clump(
  dplyr::tibble(rsid=PRSS1$SNP, pval=PRSS1$pval.exposure, id=PRSS1$id.exposure),
  clump_kb = 10000,
  clump_r2 = 0.01,
  plink_bin = "plink_Windows.exe",
  bfile = "EUR"
)
ld_PRSS1 <- PRSS1[match(exp_dat$rsid, PRSS1$SNP),]

ids<-as.character(unlist(read.table("outcome.id.txt",header=F)))
outcome_dat<- extract_outcome_data(
  snps = ld_PRSS1$SNP,
  outcomes = ids)

PRSS1_out <- format_data(outcome_dat, "outcome", snp_col = "SNP", beta_col = "beta.outcome", se_col = "se.outcome",
                        eaf_col = "eaf.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome"
)
myres <- mainMR_risk_factors(exp_dat = ld_PRSS1, out_dat = PRSS1_out, expName = "PRSS1", 
                              outName = "PRSS1_GWAS", samplesize.outcome = 10000)
fwrite(myres, paste0("PRSS1.myres.csv"))
