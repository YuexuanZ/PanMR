## This is the source code for colocalization analysis
## All data are available at https://www.jianguoyun.com/p/DZPUMToQ2IWVDBiHhqkFIAA
library(usethis)
library(devtools)
library(hyprcoloc)
library(data.table)
library(TwoSampleMR) 
library(tidyverse) 
library(vcfR)
library(dplyr)
library(data.table)

hg19 <- data.table::fread("D:/2023/Tasks/pQTL/2_Locus_Iron_pancreatic/glist-hg19")
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
pic <- data.table::fread("D:/2023/Tasks/pQTL/2_Locus_Iron_pancreatic/GCST90016676_buildGRCh37.tsv.gz")
pic_ch7 <- pic[which(pic$chromosome==7),]
pic_reg1 <- pic_ch7[which(pic_ch7$base_pair_location >= 142457318 & pic_ch7$base_pair_location <= 142460927)]   # 结局在PRSS1蛋白的hg19区域
pic_reg2 <- pic_ch7[which(pic_ch7$base_pair_location >= 141957318 & pic_ch7$base_pair_location <= 142960927)]   # 500kb
pic_reg2<-dplyr::filter(pic_reg2,!duplicated(pic_reg2$variant_id))
pic_reg3 <- pic_ch7[which(pic_ch7$base_pair_location >= 142157318 & pic_ch7$base_pair_location <= 142760927)]   # 300kb
pic_reg3<-dplyr::filter(pic_reg3,!duplicated(pic_reg3$variant_id))
pic_reg4 <- pic_ch7[which(pic_ch7$base_pair_location >= 142257318 & pic_ch7$base_pair_location <= 142660927)]   # 200kb
pic_reg4<-dplyr::filter(pic_reg4,!duplicated(pic_reg4$variant_id))
pic_reg5 <- pic_ch7[which(pic_ch7$base_pair_location >= 141457318 & pic_ch7$base_pair_location <= 143460927)]   # 1Mb
pic_reg5<-dplyr::filter(pic_reg5,!duplicated(pic_reg5$variant_id))
pic_reg6 <- pic_ch7[which(pic_ch7$base_pair_location >= 142457318 & pic_ch7$base_pair_location <= 142960927)]   # +500kb
pic_reg6<-dplyr::filter(pic_reg6,!duplicated(pic_reg6$variant_id))
# AMY2A
pic_ch1 <- pic[which(pic$chromosome==1),]
pic_reg1 <- pic_ch1[which(pic_ch1$base_pair_location >= 103159998 & pic_ch1$base_pair_location <= 105168400)]   # 1Mb
pic_reg1<-dplyr::filter(pic_reg1,!duplicated(pic_reg1$variant_id))
# CELA2A
pic_reg2 <- pic_ch1[which(pic_ch1$base_pair_location >= 14783222 & pic_ch1$base_pair_location <= 16798586)]   # 1Mb
pic_reg2<-dplyr::filter(pic_reg2,!duplicated(pic_reg2$variant_id))

###################pqtl 
# PRSS1    chr7:142,749,468-142,753,072   
PRSS1 <- data.table::fread("./3049_61_PRSS1_Trypsin.txt.gz")
PRSS1_ch7 <- PRSS1[which(PRSS1$Chrom=="chr7"),]
pqtl_reg1<- PRSS1_ch7[which(PRSS1_ch7$Pos>=142749468 & PRSS1_ch7$Pos<=142753072)]    # PRSS1蛋白在hg38的区域
pqtl_reg2<- PRSS1_ch7[which(PRSS1_ch7$Pos>=142249468 & PRSS1_ch7$Pos<=143253072)]    # 500kb
pqtl_reg2<-dplyr::filter(pqtl_reg2,!duplicated(pqtl_reg2$rsids))
pqtl_reg3<- PRSS1_ch7[which(PRSS1_ch7$Pos>=142449468 & PRSS1_ch7$Pos<=143053072)]    # 300kb
pqtl_reg3<-dplyr::filter(pqtl_reg3,!duplicated(pqtl_reg3$rsids))
pqtl_reg4<- PRSS1_ch7[which(PRSS1_ch7$Pos>=142549468 & PRSS1_ch7$Pos<=142953072)]    # 200kb
pqtl_reg4<-dplyr::filter(pqtl_reg4,!duplicated(pqtl_reg4$rsids))
pqtl_reg5<- PRSS1_ch7[which(PRSS1_ch7$Pos>=141749468 & PRSS1_ch7$Pos<=143753072)]    # 200kb
pqtl_reg5<-dplyr::filter(pqtl_reg5,!duplicated(pqtl_reg5$rsids))

# Amy    chr1:103,616,651-103,625,780
AMY2A <- data.table::fread("./18917_53_AMY2A_Pancreatic_alpha_amylase.txt.gz")
AMY2A_ch1 <- AMY2A[which(AMY2A$Chrom=="chr1"),]
pqtl_reg1<- AMY2A_ch1[which(AMY2A_ch1$Pos>=102616651 & AMY2A_ch1$Pos<=104625780)]    # 1Mb
pqtl_reg1<-dplyr::filter(pqtl_reg1,!duplicated(pqtl_reg1$rsids))

# CELA2A     chr1:15,456,728-15,472,091
CELA2A <- data.table::fread("./7140_1_CELA2A_ELA2A.txt.gz")
CELA2A_ch1 <- CELA2A[which(CELA2A$Chrom=="chr1"),]
pqtl_reg2<- CELA2A_ch1[which(CELA2A_ch1$Pos>=14456728 & CELA2A_ch1$Pos<=16472091)]    # 1Mb
pqtl_reg2<-dplyr::filter(pqtl_reg2,!duplicated(pqtl_reg2$rsids))

############################################################ hyprcoloc
df.merge=merge(pqtl_reg5,pic_reg5,by.x="rsids",by.y="variant_id")
#df.merge=merge(pqtl_data,pic,by.x="rsids",by.y="rsids")
betas=data.frame(SNPs=df.merge$rsids,T1=df.merge$Beta,T2=df.merge$beta)
ses=data.frame(SNPs=df.merge$rsids,T1=df.merge$SE,T2=df.merge$standard_error)

rownames(betas)<-betas$SNPs
rsid <- rownames(betas)
betas=subset(betas,select=-SNPs)
rownames(ses)<-ses$SNPs
ses=subset(ses,select=-SNPs)
traits <- paste0("T", 1:dim(betas)[2]);

# 这个时候需要注意的是要将我们创建的dataframe转换成矩阵
betas=data.matrix(betas)
ses=data.matrix(ses)
#矩阵之间才能来计算
res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid,snpscores = TRUE,
                 bb.selection = "regional",reg.steps = 1,reg.thresh = 0.5,align.thresh = 0.5,
                 prior.1 = 1e-4,prior.c = 0.02,prior.12 = 0.02,sensitivity = FALSE,sense.1 = 1,sense.2 = 2)
res[[1]]
cred.sets(res, value = 0.95);
ld.matrix["rs533406","rs7404039"]

########################################################### coloc
library(coloc)
input <- merge(pqtl_reg2,pic_reg2,by.x="rsids",by.y="variant_id",all=FALSE,suffixes=c("_pqtl","_gwas"))
# dataset1 GWAS  dataset2 pQTL
result <- coloc.abf(dataset1=list(pvalues=input$Pval, type="quant", N=25617, snp=input$rsids),     # GWAS # cc:case-control type
                    dataset2=list(pvalues=input$p_value, type="quant", N=35287, snp=input$rsids), MAF=input$ImpMAF,
                    p1=1e-6,p2=0.1,
                    p12=0.1)    # quant:连续型
# PRSS1
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#9.59e-08  4.06e-02  1.97e-08  8.32e-03  9.51e-01 
#[1] "PP abf for shared variant: 95.1%"
# AMY2A
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.48e-295  8.35e-02 1.32e-295  7.43e-02  8.42e-01 
#[1] "PP abf for shared variant: 84.2%"
# CELA2A
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.63e-24  2.87e-02  6.45e-25  7.06e-03  9.64e-01 
#[1] "PP abf for shared variant: 96.4%"

result
library(dplyr)
need_result=result$results %>% filter(SNP.PP.H4 > 0.80)
# rs62470638 PRSS1




########################################################## PLOT


