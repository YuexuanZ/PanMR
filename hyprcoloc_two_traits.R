## This is the source code for colocalization analysis
## All data are available at https://www.jianguoyun.com/p/DZPUMToQ2IWVDBiHhqkFIAA
library(usethis)
library(devtools)
library(data.table)
library(TwoSampleMR) 
library(tidyverse) 
library(vcfR)
library(dplyr)
library(data.table)

hg19 <- data.table::fread("glist-hg19")
hg19[which(hg19$V4=="PRSS1"),]    
#V1        V2        V3    V4
#1:  7 142457318 142460927 PRSS1
hg19[which(hg19$V4=="AMY2A"),]
#V1        V2        V3    V4
#1:  1 104159998 104168400 AMY2A
hg19[which(hg19$V4=="CELA2A"),]
#V1       V2       V3     V4
#1:  1 15783222 15798586 CELA2A

###################outcome
pic <- data.table::fread("Pancreas_iron_content_GCST90016676.tsv.gz")
pic_ch7 <- pic[which(pic$chromosome==7),]
# PRSS1
pic_reg <- pic_ch7[which(pic_ch7$base_pair_location >= 141957318 & pic_ch7$base_pair_location <= 142960927)]
pic_reg<-dplyr::filter(pic_reg,!duplicated(pic_reg$variant_id))
# AMY2A
pic_ch1 <- pic[which(pic$chromosome==1),]
pic_reg <- pic_ch1[which(pic_ch1$base_pair_location >= 103659998 & pic_ch1$base_pair_location <= 104668400)] 
pic_reg<-dplyr::filter(pic_reg,!duplicated(pic_reg$variant_id))
# CELA2A
pic_reg <- pic_ch1[which(pic_ch1$base_pair_location >= 15283222 & pic_ch1$base_pair_location <= 16298586)] 
pic_reg<-dplyr::filter(pic_reg,!duplicated(pic_reg$variant_id))

###################pqtl 
# PRSS1    chr7:142,749,468-142,753,072   
PRSS1 <- data.table::fread("Instruments_PRSS1.txt.gz")
PRSS1_ch7 <- PRSS1[which(PRSS1$Chrom=="chr7"),]
pqtl_reg<- PRSS1_ch7[which(PRSS1_ch7$Pos>=142249468 & PRSS1_ch7$Pos<=143253072)] 
pqtl_reg<-dplyr::filter(pqtl_reg,!duplicated(pqtl_reg$rsids))

# Amy    chr1:103,616,651-103,625,780
AMY2A <- data.table::fread("Instruments_AMY2A.txt.gz")
AMY2A_ch1 <- AMY2A[which(AMY2A$Chrom=="chr1"),]
pqtl_reg<- AMY2A_ch1[which(AMY2A_ch1$Pos>=103116651 & AMY2A_ch1$Pos<=104125780)]    
pqtl_reg<-dplyr::filter(pqtl_reg,!duplicated(pqtl_reg$rsids))

# CELA2A     chr1:15,456,728-15,472,091
CELA2A <- data.table::fread("Instruments_CELA2A.txt.gz")
CELA2A_ch1 <- CELA2A[which(CELA2A$Chrom=="chr1"),]
pqtl_reg<- CELA2A_ch1[which(CELA2A_ch1$Pos>=14956728 & CELA2A_ch1$Pos<=15972091)]    
pqtl_reg<-dplyr::filter(pqtl_reg,!duplicated(pqtl_reg$rsids))

################## coloc
library(coloc)
input <- merge(pqtl_reg,pic_reg,by.x="rsids",by.y="variant_id",all=FALSE,suffixes=c("_pqtl","_gwas"))
# dataset1 GWAS  dataset2 pQTL
result <- coloc.abf(dataset1=list(pvalues=input$Pval, type="quant", N=25617, snp=input$rsids),
                    dataset2=list(pvalues=input$p_value, type="quant", N=35287, snp=input$rsids), MAF=input$ImpMAF,
                    p1=1e-6,p2=0.1,
                    p12=0.1) 
