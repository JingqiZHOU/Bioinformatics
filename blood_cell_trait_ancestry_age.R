###############2020.8.16#######################################################################
###############To identify blood cell traits that are different across ancestry groups#########

###############1.summary statistics#############################################################
library(broom)
load("/scratch/jz64675/neald/Blood_08_05_unQC_bd_pheno.RData")
bd_pheno$FID<-as.character(bd_pheno$FID)
colnames(bd_pheno)[31:61] <- gsub("[(]", "", colnames(bd_pheno)[31:61])
colnames(bd_pheno)[31:61] <- gsub("[)]", "", colnames(bd_pheno)[31:61])

EURanc<-bd_pheno[bd_pheno$Race %in% c("White","British"),]
AFIanc<-bd_pheno[bd_pheno$Race %in% c("African","Black or Black British","Any other Black Background"),]
ASIanc<-bd_pheno[bd_pheno$Race %in% c("Asian or Asian British","Chinese","Pakistani","Any other Asian background","Indian","Bangladeshi"),]

> dim(EURanc)
[1] 443174     82
> dim(AFIanc)
[1] 3421   82
> dim(ASIanc)
[1] 11456    82

stat<-function(x,a){
  mean<-apply(x[,31:61],2,mean,na.rm=TRUE)
  sd<-apply(x[,31:61],2,sd,na.rm=TRUE)
  sta<-cbind(mean,sd)
  colnames(sta)<-c(paste0(a,"_mean"),paste0(a,"_sd"))
  sta
}

#EURanc_mean<-apply(EURanc[,31:61],2,mean,na.rm=TRUE)
#EURanc_sd<-apply(EURanc[,31:61],2,sd,na.rm=TRUE)
#eur_sta <- cbind(EURanc_mean,EURanc_sd)

eur_sta<-stat(EURanc,"eur")
afr_sta<-stat(AFIanc,"african")
asi_sta<-stat(ASIanc,"asian")
sta_combine<-cbind(eur_sta,afr_sta,asi_sta)
write.table(sta_combine,file="sta_combine.txt",quote = FALSE,sep="\t",row.names = TRUE)

###############2. wilcoxon test#############################################################
library(dplyr)
library(tidyverse)
library(rstatix)
library(ggpubr)

EURanc$popu<-"EUR"
AFIanc$popu<-"AFI"
afr_com<-rbind(EURanc,AFIanc)

#afr_com %>% group_by(popu) %>% get_summary_stats(Platelet_count, type = "median_iqr")
#stat.test <- afr_com %>% wilcox_test(Monocyte_count ~ popu) %>% add_significance()

myresults<-list()
myvars <- colnames(afr_com)[31:61]
for (i in 1:length(myvars)){
  myresults[[i]] <- afr_com %>% wilcox_test(as.formula(paste0(myvars[i], "~", "popu")) %>% add_significance())
}
res<-do.call(rbind, myresults)
write.table(res,file="african_wilcox.txt",sep="\t",row.names = FALSE)


ASIanc$popu<-"ASI"
asi_com<-rbind(EURanc,ASIanc)
myresults_asian<-list()
myvars <- colnames(asi_com)[31:61]
for (i in 1:length(myvars)){
  myresults_asian[[i]] <- asi_com %>% wilcox_test(as.formula(paste0(myvars[i], "~", "popu")) %>% add_significance())
}
res1<-do.call(rbind, myresults_asian)
write.table(res1,file="asian_wilcox.txt",sep="\t",row.names = FALSE)

###############3. logit regression#############################################################
lm_reg <- function(x){
afr_com <- x
colnames(afr_com)[31:61] <- gsub("[(]", "", colnames(afr_com)[31:61])
colnames(afr_com)[31:61] <- gsub("[)]", "", colnames(afr_com)[31:61])
blood_afr_com<- scale(afr_com[,31:61],center=T,scale=T)
new_afr_com <- cbind(afr_com[,-c(31:61)],blood_afr_com)

myvars <- colnames(afr_com)[31:61]
fit <- lapply(myvars, function(dvar)
tidy(lm(as.formula(paste0("popu", " ~ ", paste0(dvar), "+Sex + Age")) ,data=new_afr_com)))

myresults = data.frame()
len2 <- length(fit)
for (i in 1:len2){
  temp = fit[[i]][c(2),]
  temp[c(1),1] <- c("Trait")
  temp$phenotype = fit[[i]]$term[2]
  myresults = rbind(myresults,temp)
}
myresults
}

res<-lm_reg(asi_com)
write.table(res,file="blood_asian.txt",sep="\t",row.names = FALSE)



###############4. blood trait Age regression#############################################################
library(broom)
load("/scratch/jz64675/neald/Blood_08_05_unQC_bd_pheno.RData")
bd_pheno$FID<-as.character(bd_pheno$FID)
colnames(bd_pheno)[31:61] <- gsub("[(]", "", colnames(bd_pheno)[31:61])
colnames(bd_pheno)[31:61] <- gsub("[)]", "", colnames(bd_pheno)[31:61])

blood_com<- scale(bd_pheno[,31:61],center=T,scale=T)
new_com <- cbind(bd_pheno[,-c(31:61)],blood_com)

myvars <- colnames(bd_pheno)[31:61]

fit <- lapply(myvars, function(dvar)
tidy(lm(as.formula(paste0(paste0(dvar), " ~ ", "Sex + Age")),family = binomial(link = "logit"),data=new_com)))

fit <- lapply(myvars, function(dvar)
tidy(lm(as.formula(paste0(paste0(dvar), " ~ ", "Sex + Age+Sex*Age")),data=new_com)))


myresults = data.frame()

len2 <- length(fit)
for (i in 1:len2){
  temp = fit[[i]][2,]  ######line for covariant
  temp[c(1),1] <- myvars[i]
  temp$phenotype = fit[[i]]$term[2]  #####line for covariant
  myresults = rbind(myresults,temp)
}
myresults

write.table(myresults,file="blood_sex1_0913_2.txt",sep="\t",row.names = FALSE)

save(new_com,file="new_blood.RData")
save(bd_pheno,file="bd_pheno_blood.RData")

#########################scatter Plot#########################
new_com$Sex <- as.factor(new_com$Sex)
library(ggplot2)
ggplot(new_com, aes(x=Age, y=Monocyte_count, color=Sex)) +
  geom_point()
   +
  geom_smooth(method=lm)
# Remove confidence intervals
# Extend the regression lines
ggplot(mtcars, aes(x=wt, y=mpg, color=cyl, shape=cyl)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

dev.off()


####################5. blood_SNP_interaction###########
module load R/3.6.2-foss-2018a-X11-20180131-GACRC
R
file.remove('.RData')

library(broom)

load("/scratch/jz64675/MR/covid4all/biomarker_covatiate.RData")
load("/scratch/jz64675/neald/query_phecov.RData")
load("/scratch/jz64675/neald/Blood_08_05_unQC_bd_pheno.RData")

bd_pheno$FID<-as.character(bd_pheno$FID)
bd_db2$id <-as.character(bd_db2$id)
blood_combine <- merge(bd_db2[,c(1,1866,1847:1848)],bd_pheno,by.x="id",by.y="FID")

blood_scale<- scale(blood_combine[,34:64],center=T,scale=T)
new_blood <- cbind(blood_combine[,-c(34:64)],blood_scale)

new_blood$Assessment_centres <- as.factor(new_blood$Assessment_centres)
#new_blood$rs67959919 <- as.character(new_blood$rs67959919)

colnames(new_blood)[55:85] <- gsub("[(]", "", colnames(new_blood)[55:85])
colnames(new_blood)[55:85] <- gsub("[)]", "", colnames(new_blood)[55:85])
myvars <- colnames(new_blood)[55:85]

fit <- lapply(myvars, function(dvar)
    tidy(lm(as.formula(paste0(paste0(dvar), " ~ ", "rs67959919",
    "+rs67959919*Sex + Sex + Array+Assessment_centres+age+ PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10 "))
    ,data=new_blood)))

myresults = data.frame()
len2 <- length(fit)
for (i in 1:len2){
  temp = fit[[i]][36,]
  temp[c(1),1] <- c("Main")
  temp$phenotype = fit[[i]]$term[36]
  temp$covariant = myvars[i]
  myresults = rbind(myresults,temp)
}
#
write.table(myresults,file="blood_norm_rs67959919_sex.txt",sep="\t",row.names = FALSE)


#######################
load("/scratch/jz64675/MR/covid4all/biomarker_covatiate.RData")
load("/scratch/jz64675/neald/query_phecov.RData")
load("/scratch/jz64675/neald/Blood_08_05_unQC_bd_pheno.RData")

bd_pheno$FID<-as.character(bd_pheno$FID)
bd_db2$id <-as.character(bd_db2$id)
blood_combine <- merge(bd_db2[,c(1,1866,1847:1848)],bd_pheno,by.x="id",by.y="FID")

blood_scale<- scale(blood_combine[,34:64],center=T,scale=T)
new_blood <- cbind(blood_combine[,-c(34:64)],blood_scale)

library(dplyr)
myvar <- colnames(blood_combine)[34:64]
fit <- lapply(myvar, function(dvar)
blood_combine %>% group_by(Sex) %>% summarize(mean(dvar,na.rm=T),sd(dvar,na.rm=T)) )

myresults = data.frame()
for(i in 1:length(myvar)){
a<-blood_combine %>% group_by(Sex) %>% summarize(mean(myvar[i],na.rm=T),sd(myvar[i],na.rm=T))
a2<-t(a)
myresults<-rbind(myresults,a2)
}
write.table(myresults,file="blood_mean.txt",sep="\t")


a<-blood_combine %>% group_by(Sex) %>% summarise_at(vars(Platelet_count:Lymphocyte_percentage),list(m=mean,s=sd),na.rm=T)
a2<-t(a)
write.table(a2,file="blood_mean.txt",sep="\t")

b <-blood_combine %>% summarise_at(vars(Platelet_count:Lymphocyte_percentage),list(m=mean,s=sd),na.rm=T)
b2<-t(b)
write.table(b2,file="blood_mean2.txt",sep="\t")


library(dplyr)
myvars <- colnames(new_blood)[55:85]
###310998
myresults = data.frame()
a < -colSums(is.na(new_blood[55:85]))
a<-as.data.frame(a)
write.table(a,file="blood_sample_num.txt",sep="\t")


myvars <- colnames(new_biom2)[47:76]
> dim(new_biom2)
[1] 389620     83
a <- colSums(is.na(new_biom2[47:76]))
a<-as.data.frame(a)
write.table(a,file="blood_biomarker_sample_num.txt",sep="\t")
