#####################################################################################################
#########testing phecode association::::covid case vs ALL case 46
#####################################################################################################

library(broom)
library(dplyr)
library(reshape2)

load("phecode_05_24_bd_pheno.RData")
covid3 <- read.table("ukb_covid_out3.txt",header = TRUE,sep="\t")
covid3$ID <- as.character(covid3$ID)

bd_pheno$IID <- as.character(bd_pheno$IID)
new_covid <- merge(bd_pheno,covid3,by.x="IID",by.y="ID",all.x=T)
new_covid <-within(new_covid,{result3 <- ifelse(is.na(result3),0,result3)})

> dim(new_covid)
[1] 468418   1888
table(new_covid$result3)
0      1
467151   1267

bm <- new_covid
colnames(bm) <- make.names(colnames(bm),unique = TRUE)

# remove those with only True or False
bm <- bm[,which(sapply(bm, function(x)
  length(unique(x[!is.na(x)])) >= 2))]

> dim(bm)
[1] 468418   1714
# bm$center <- as.factor(bm$center)
# define variable
len1 <- length(colnames(bm)) - 1
myvars <- colnames(bm)[39:len1]

fit <- lapply(myvars, function(dvar)
  tidy(glm(as.formula(paste0("result3", " ~ ", paste0(dvar," "),
                             "Sex + Assessment_centres+Age+ PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + PCA19 + PCA20")),
           family = binomial(link = "logit"),data=bm)))

myresults = data.frame()
len2 <- length(fit)
for (i in 1:len2){
  temp = fit[[i]][c(2),]
  temp[c(1),1] <- c("Main")
  temp$phenotype = fit[[i]]$term[2]
  myresults = rbind(myresults,temp)
}
#

bd <- myresults
library(stringr)
bd$phecode <- str_extract(bd$phenotype,"\\-*\\d+\\.*\\d*")

library(PheWAS)
bd$phenotype <- bd$phecode
bd <- bd[,-6]
bd <- addPhecodeInfo(bd)

or_value <- data.frame(or_main = exp(bd$estimate))
low_value <- data.frame(low_main =
                          round(exp(bd$estimate-1.96*bd$std.error),2))
high_value <- data.frame(high_main =
                           round(exp(bd$estimate+1.96*bd$std.error),2))
or_ci <- data.frame(ci_main = paste0('(',low_value$low_main,',',high_value$high_main,')'))

bd <- cbind(bd,or_value,or_ci)
# Count case/control
sumtrue <- function(x){
  x <- na.omit(x)
  length(x[which(x)])
}

sumfalse <- function(x){
  x <- na.omit(x)
  length(x[which(!x)])
}
#allP <- allP[which(allP$sex == 'Female'),]
allphenotype = data.frame(case = apply(bm[,39:len1], 2, sumtrue),
                          control = apply(bm[,39:len1], 2, sumfalse))
allphenotype$phenotype <- row.names(allphenotype)
allphenotype$phenotype <- gsub("X","",allphenotype$phenotype)

bd <- bd %>% inner_join(allphenotype,by = c("phecode" = "phenotype"))

bm1 <- bm[which(bm$result3 == 1),]
bd_phenohenotype = data.frame(positive_case = apply(bm1[,39:len1], 2, sumtrue),
                              positive_control = apply(bm1[,39:len1], 2, sumfalse))
bd_phenohenotype$phenotype <- row.names(bd_phenohenotype)
bd_phenohenotype$phenotype <- gsub("X","",bd_phenohenotype$phenotype)
bd <- bd %>% inner_join(bd_phenohenotype,by = c("phecode" = "phenotype"))

bm2 <- bm[which(bm$result3 == 0),]
bd_phenohenotype = data.frame(nonpositive_case = apply(bm2[,39:len1], 2, sumtrue),
                             nonposituve_control = apply(bm2[,39:len1], 2, sumfalse))
bd_phenohenotype$phenotype <- row.names(bd_phenohenotype)
bd_phenohenotype$phenotype <- gsub("X","",bd_phenohenotype$phenotype)

bd <- bd %>% inner_join(bd_phenohenotype,by = c("phecode" = "phenotype"))

bd$phenotype <- paste0("'",bd$phecode)
write.table(bd,file = "covid_positive_vs_all_46.txt",sep = "\t",row.names = FALSE)



####################################plot phecode inpatient 46######################################################
# Set workspace
setwd("/Users/yelab/Desktop/phewas3/")
# Load package for reading xlsx without using java
library(tidyverse)
# Read data
data <- readxl::read_xlsx("inpatient_phecode_asso_46.xlsx",sheet = "Sheet1",col_names = TRUE)

# Manually give the most significant trait a small p value
#data$p.value_Main[1] <- 5.68e-13
data$p.value_fdr<-p.adjust(data$p.value, method = "fdr", n = length(data$p.value))
p_fdr_min <- unlist(data[which(data$p.value_fdr<0.05&data$p.value_fdr>0.04),5])
max(data[which(data$p.value_fdr<0.05),5])

# Load package for drawing
library(ggplot2)
library(RColorBrewer)
# Set color
colourCount = length(unique(mtcars$hp))
nb.cols <- 18

mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
values = colorRampPalette(brewer.pal(8,"Accent"))(colourCount)

# Data preprocessing, adding POS == x, order for beautiful scatter
data <- data[order(data$description),]
data <- data[order(data$group),]
data$POS <- seq(1,nrow(data),1)
# Set where to add x text
group <- tapply(data$POS, data$group, median)
# Change data type
data$or_main <- as.numeric(data$or_main)
# Set shape
data$size <- ifelse(data$or_main > 1,24,25)
# Change data type
data$p.value <-as.numeric(data$p.value)
# Set where to divide to 2 pics
data$bin <- ifelse(-log(data$p.value,10) > 25,1,0)
# Set for hlines
difgrid <- data.frame(bin = c(0,0),values = c(-log10(0.05),-log10(0.05/122)))
# Set layer to distinguish
my_limits <- function(x) { if (max(x) < 24) c(-0.2,10) else c(140,155)}
my_breaks <- function(x) { if (max(x) < 24) seq(0,10,2) else seq(140,152,2)}

# Draw
p <- ggplot(data,aes(POS,-log10(p.value))) +
  theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "transparent")) +
  geom_point(aes(color = as.factor(group),fill = as.factor(group),size = abs(or_main)*0.5),shape = data$size,show.legend = FALSE) +
  facet_grid(cols = vars(bin),scales = "free",space = "free") +
  theme(strip.background = element_blank(), strip.text.x  = element_blank()) +
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  labs(x = NULL, y = expression(''~-log[10]~'(P)')) +
  theme(axis.text.x = element_text(color = "black", size = 12 ,angle = 45, hjust = 1, vjust = 1)) +
  scale_x_continuous(breaks = group,labels = names(group),limits = c(-3,125),expand = c(0,0)) +
  scale_y_continuous(breaks = my_breaks,limits = my_limits, expand = c(0,0)) +
  # geom_hline(yintercept = c(-log10(0.05/850)), color = c('gray'), linetype = "dotted",size = 0.5) +
  geom_hline(data = difgrid,aes(yintercept = values), color = c('gray','red'),linetype = "dashed",size = 0.5)



# Add thoese needing description
# Load librry for text
library(ggrepel)
data$annotate <- ifelse(data$p.value < 0.05,"yes","no")
p +
  geom_text_repel(data=data[data$annotate=="yes",], aes(label=as.factor(description)),  color = "black",
                  size=3, force=1.3,show.legend = FALSE,direction = 'y',point.padding = 0.5) +
  ggtitle("phecode inpatient vs non-inpatient plot") +
  theme(plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5))
