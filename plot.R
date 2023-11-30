## Forest Plot
## The source data are the outputs of main_MR.R
BiocManager::install("forestplot")
library(forestplot)
library(forestploter)
library(grid)

############ PRSS1
Data_str <- read.csv("data_1.csv")
## Odds Ratio and 95%CI
Data_str$`Odds Ratio (95%CI)` <- ifelse(is.na(Data_str$OR), "",
                                        sprintf("%.2f (%.2f - %.2f)",
                                                Data_str$OR, Data_str$lower, Data_str$upper))

Data_str$` ` <- paste(rep(" ", 25), collapse = " ")
Data_str$P <- ifelse(is.na(Data_str$P), "", Data_str$P)
Data_str[3,5]<-"0.00"
tm <- forest_theme(base_size = 9, 
                   ci_pch = 16,   
                   ci_col = "#4575b4", 
                   ci_lty = 1,   
                   ci_lwd = 1.5,    
                   ci_Theight = 0.2, 
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   arrow_length = 0.05,
                   arrow_cex = 0.7)

p1 <- forest(Data_str[, c(1,7,6,5)],
             est = Data_str$OR,
             lower = Data_str$lower,
             upper = Data_str$upper,
             ci_column = 2,
             ref_line = 1,   
             arrow_lab = c("PRSS1 downregulation raises risk", "PRSS1 upregulation raises risk"),
             xlim = c(0.5,1.5),
             ticks_at = c(0.5,0.75,1,1.25,1.5),
             theme = tm)
p1 <- edit_plot(p1,
                row = 1,
                gp = gpar(fontface = "bold"))
p1 <- add_border(p1, part = "header", row = 2, where = "top")
p1

############ PRSS1

Data_str2 <- read.csv("data_2.csv")
Data_str2$`Beta (95%CI)` <- ifelse(is.na(Data_str2$Beta), "",
                                   sprintf("%.2f (%.2f - %.2f)",
                                           Data_str2$Beta, Data_str2$lower, Data_str2$upper))
Data_str2$` ` <- paste(rep(" ", 25), collapse = " ")
Data_str2[c(1,7),5]<-""
Data_str2[4,5]<-"0.10"
Data_str2[c(9,8),5]<-"0.00"

tm <- forest_theme(base_size = 9, 
                   ci_pch = 16,
                   ci_col = "#762a83",
                   ci_lty = 1,   
                   ci_lwd = 1.5,  
                   ci_Theight = 0.2, 
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   arrow_length = 0.05,
                   arrow_cex = 0.7) 
p2 <- forest(Data_str2[, c(1,7,6,5)],
             est = Data_str2$Beta,
             lower = Data_str2$lower,
             upper = Data_str2$upper,
             ci_column = 2,
             ref_line = 0,   
             arrow_lab = c("PRSS1 downregulation raises risk", "PRSS1 upregulation raises risk"),
             xlim = c(-0.25,0.25),
             ticks_at = c(-0.25,-0.1,0,0.1,0.25),
             theme = tm)
p2 <- edit_plot(p2,
                row = c(1,7),
                gp = gpar(fontface = "bold"))

p2 <- add_border(p2, part = "header", row = 2, where = "top")
p2

############## CELA2A
Data_str3 <- read.csv("data_3.csv")
Data_str3$`Odds Ratio (95%CI)` <- ifelse(is.na(Data_str3$OR), "",
                                         sprintf("%.2f (%.2f - %.2f)",
                                                 Data_str3$OR, Data_str3$lower, Data_str3$upper))
Data_str3$` ` <- paste(rep(" ", 25), collapse = " ")
Data_str3[1,5]<-""
Data_str3[4,5]<-"0.90"
Data_str3[7,5]<-"0.60"
tm <- forest_theme(base_size = 9,
                   ci_pch = 16,  
                   ci_col = "#4575b4",
                   ci_lty = 1, 
                   ci_lwd = 1.5,   
                   ci_Theight = 0.2, 
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   arrow_length = 0.05,
                   arrow_cex = 0.7)
p3 <- forest(Data_str3[, c(1,7,6,5)],
             est = Data_str3$OR,
             lower = Data_str3$lower,
             upper = Data_str3$upper,
             ci_column = 2,
             ref_line = 1,  
             arrow_lab = c("CELA2A downregulation raises risk", "CELA2A upregulation raises risk"),
             xlim = c(0.5,1.5),
             ticks_at = c(0.5,0.75,1,1.25,1.5),
             theme = tm)
p3 <- edit_plot(p3,
                row = 1,
                gp = gpar(fontface = "bold"))
p3 <- add_border(p3, part = "header", row = 2, where = "top")
p3

############## CELA2A
Data_str4 <- read.csv("data_4.csv")
Data_str4$`Beta (95%CI)` <- ifelse(is.na(Data_str4$Beta), "",
                                   sprintf("%.2f (%.2f - %.2f)",
                                           Data_str4$Beta, Data_str4$lower, Data_str4$upper))
Data_str4$` ` <- paste(rep(" ", 25), collapse = " ")
Data_str4[c(1,7),5]<-""
Data_str4[4,5]<-"0.20"
Data_str4[c(9,8),5]<-"0.00"

tm <- forest_theme(base_size = 9, 
                   ci_pch = 16,
                   ci_col = "#762a83",
                   ci_lty = 1,  
                   ci_lwd = 1.5,  
                   ci_Theight = 0.2,
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   arrow_length = 0.05,
                   arrow_cex = 0.7) 
p4 <- forest(Data_str4[, c(1,7,6,5)],
             est = Data_str4$Beta,
             lower = Data_str4$lower,
             upper = Data_str4$upper,
             ci_column = 2,
             ref_line = 0,
             arrow_lab = c("CELA2A downregulation raises risk", "CELA2A upregulation raises risk"),
             xlim = c(-0.2,1),
             ticks_at = c(-0.2,0,0.5,1),
             theme = tm)
p4 <- edit_plot(p4,
                row = c(1,7),
                gp = gpar(fontface = "bold"))

p4 <- add_border(p4, part = "header", row = 2, where = "top")
p4


############# AMY2A 
Data_str5 <- read.csv("data_5.csv")
Data_str5$`Odds Ratio (95%CI)` <- ifelse(is.na(Data_str5$OR), "",
                                         sprintf("%.2f (%.2f - %.2f)",
                                                 Data_str5$OR, Data_str5$lower, Data_str5$upper))
Data_str5$` ` <- paste(rep(" ", 25), collapse = " ")
Data_str5[1,5]<-""
Data_str5[3,5]<-"0.00"
tm <- forest_theme(base_size = 9, 
                   ci_pch = 16, 
                   ci_col = "#4575b4", 
                   ci_lty = 1,   
                   ci_lwd = 1.5, 
                   ci_Theight = 0.2, 
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   arrow_length = 0.05,
                   arrow_cex = 0.7)

p5 <- forest(Data_str5[, c(1,7,6,5)],
             est = Data_str5$OR,
             lower = Data_str5$lower,
             upper = Data_str5$upper,
             ci_column = 2,
             ref_line = 1,  
             arrow_lab = c("AMY2A downregulation raises risk", "AMY2A upregulation raises risk"),
             xlim = c(0.25,1.75),
             ticks_at = c(0.5,1,1.5),
             theme = tm)
p5 <- edit_plot(p5,
                row = 1,
                gp = gpar(fontface = "bold"))
p5 <- add_border(p5, part = "header", row = 2, where = "top")
p5

################### AMY2A 
Data_str6 <- read.csv("data_6.csv")
Data_str6$`Beta (95%CI)` <- ifelse(is.na(Data_str6$Beta), "",
                                   sprintf("%.2f (%.2f - %.2f)",
                                           Data_str6$Beta, Data_str6$lower, Data_str6$upper))
Data_str6$` ` <- paste(rep(" ", 25), collapse = " ")
Data_str6[c(1,7),5]<-""
Data_str6[c(9,8),5]<-"0.00"

tm <- forest_theme(base_size = 9, 
                   ci_pch = 16,  
                   ci_col = "#762a83",
                   ci_lty = 1,  
                   ci_lwd = 1.5,   
                   ci_Theight = 0.2, 
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   arrow_length = 0.05,
                   arrow_cex = 0.7) 
p6 <- forest(Data_str6[, c(1,7,6,5)],
             est = Data_str6$Beta,
             lower = Data_str6$lower,
             upper = Data_str6$upper,
             ci_column = 2,
             ref_line = 0,  
             arrow_lab = c("AMY2A downregulation raises risk", "AMY2A upregulation raises risk"),
             xlim = c(-0.1,0.2),
             ticks_at = c(-0.1,0,0.1,0.2),
             theme = tm)
p6 <- edit_plot(p6,
                row = c(1,7),
                gp = gpar(fontface = "bold"))

p6 <- add_border(p6, part = "header", row = 2, where = "top")
p6

## Manhattan Plot
library(data.table)
library(CMplot)
library(openxlsx)

hg19_data <- data.table::fread("glist-hg19")
snp_data <- data.table::fread("SNP.txt")
MR_data <- read.xlsx("result_Iceland.xlsx",sheet = "mr")
MR_data1 <- merge(MR_data, hg19_data, by.x = "exposure", by.y = "V4")
MR_data2 <- MR_data1[,c(1,10,11,9)]
protein1 <- MR_data2[order(MR_data2$pval),]
protein_top <- protein1$exposure[1:25]

set.seed(666666)
CMplot(MR_data2, plot.type="m",LOG10=TRUE,
       col=c("#3E0A52", "#423D77"),  # color
       highlight=protein_top,
       highlight.cex=1,highlight.pch=c(15:17), 
       highlight.text=protein_top,      
       amplify=FALSE,file="tiff",file.name="PIC",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
