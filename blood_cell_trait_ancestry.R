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
