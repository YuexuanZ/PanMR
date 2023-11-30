## This is the source code for the analysis results of Figure 3 and the reverse MR analysis.
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

## For Figure 3
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


## Pancreatic benign neoplasm: 
## Pancreatic cancer: finn-b-CD2_BENIGN_PANCREAS_EXALLC
## Chronic pancreatitis: ieu-a-822
## Alcohol-induced chronic pancreatitis: finn-b-K11_CHRONPANC
## Acute pancreatitis: finn-b-ALCOPANCCHRON
## Acohol-induced acute pancreatitis: finn-b-K11_ACUTPANC
## Serum Iron Measurement: finn-b-ALCOPANCACU
## Ferritin: prot-a-1149
## Transferrin: prot-c-4162_54_2
## Lactotransferrin: prot-a-1808
## Transferrin Saturation: ieu-a-1051
## Pancreas fat: ebi-a-GCST90016675
## Pancreas volume: ebi-a-GCST90016669
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

## For reverse MR of PRSS1, AMY2A and CELA2A
pqtl_data <- as.data.frame(data.table::fread("Instruments_PRSS1.txt.gz"))
pqtl_data2 <- as.data.frame(data.table::fread("assocvariants.annotated.txt.gz"))
pqtl_data$eaf <- pqtl_data2$effectAlleleFreq[match(pqtl_data$Name, pqtl_data2$Name)]

pqtl <- TwoSampleMR::format_data(pqtl_data, "exposure", snp_col = "rsids", beta_col = "Beta", se_col = "SE",
                                    effect_allele_col = "effectAllele", other_allele_col = "otherAllele", 
                                    eaf_col = "eaf", pval_col = "Pval")
p_threshold <- 5e-8
kb <- 10000
r2 <- 0.001
pqtl_data_2 <- pqtl[which(pqtl$pval.exposure < p_threshold), ]
exp_dat <- ld_clump(
  dplyr::tibble(rsid=pqtl_data_2$SNP, pval=pqtl_data_2$pval.exposure, id=pqtl_data_2$id.exposure),
  clump_kb = kb,
  clump_r2 = r2,
  plink_bin = "plink_Windows.exe",
  bfile = "EUR"
)
ld_pqtl <- pqtl_data_2[match(exp_dat$rsid, pqtl_data_2$SNP),]

pic_data <- data.table::fread("Pancreas_iron_content_GCST90016676.tsv.gz")
pic <- TwoSampleMR::format_data(pic_data, "exposure", snp_col = "rsids", beta_col = "Beta", se_col = "SE",
                                    effect_allele_col = "effectAllele", other_allele_col = "otherAllele", 
                                    eaf_col = "eaf", pval_col = "Pval")
p_threshold <- 0.05
kb <- 10000
r2 <- 0.001
pic_data_2 <- pic[which(pic$pval.exposure < p_threshold), ]
exp_dat <- ld_clump(
  dplyr::tibble(rsid=pic_data_2$SNP, pval=pic_data_2$pval.exposure, id=pic_data_2$id.exposure),
  clump_kb = kb,
  clump_r2 = r2,
  plink_bin = "plink_Windows.exe",
  bfile = "EUR"
)
ld_pic <- pic_data_2[match(exp_dat$rsid, pic_data_2$SNP),]
position <- Reduce(intersect,list(ld_pqtl$exposure,ld_pic$exposure))
overlap <- data1[which(data1$exposure == position[1]),]

position <- Reduce(intersect,list(ld_pqtl$exposure,ld_pic$exposure))
overlap <- data1[which(data1$exposure == position[1]),]
pic_snp <- ld_pic[-which(ld_pic$exposure == position[1]),]
for (i in 2:length(position)){
  pic_snp <- pic_snp[-which(pic_snp$exposure == position[i]),]
}
pqtl_out2 <- pqtl_data[which(pqtl_data$variant_id%in%pic_snp$SNP), ]
pqtl_out2 <- format_data(pqtl_out2, "outcome", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error",
                        eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value"
)
myres1 <- mainMR_risk_factors(exp_dat = pic_snp, out_dat = pqtl_out2, expName = "pic", 
                              outName = "prss1", samplesize.outcome = 10000)
fwrite(myres1, paste0("pic.prss1.csv"))
