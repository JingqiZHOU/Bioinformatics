module load PLINK/2.00-alpha2-x86_64-20191128

#########1.query snp ##############
plink2 --pfile /scratch/jz64675/ukb/ukb_imp_chr3 --extract querysnp.txt --recode vcf --out qury_chr3_genotype
sed '1,6d' qury_chr3_genotype.vcf > query_snp.txt
awk '{for (i=1; i<=NF; i++)  { a[NR,i] = $i}}
      NF>p { p = NF }
      END {
      for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        } print str }
      }' query_snp.txt > snp.txt
sed '1,2d' snp.txt > snp2.txt

library(dplyr)
library(stringr)
bd <- read.table("/scratch/jz64675/neald/snp2.txt",header = TRUE)
bd_temp <- bd[1:2,]
rownames(bd_temp) <- bd_temp[,1]
bd_temp <- bd_temp[,-1]

temp <- bd[c(-1:-6),]
temp1 <- as.data.frame(str_split_fixed(as.character(temp$ID),"_",2))
temp$ID <- temp1$V1

temp <- apply(temp,2,function(x){
  y <- str_replace_all(x,"0/0","0")
  y <- str_replace_all(y,"0/1","1")
  y <- str_replace_all(y,"1/1","2")
  y
})

temp[,2:dim(temp)[2]] <- apply(temp[,2:dim(temp)[2]],2,as.numeric)
bd_temp <- as.data.frame(temp)
save(bd_temp,file="query_snp.RData")


#########2.do phewas############
load("/scratch/jz64675/MR/tmprss2_out2/phenotype_cov2.RData")

###[1846:1860]  covariant

bd_temp$ID <- as.character(bd_temp$ID)
phenotype2$id <- as.character(phenotype2$id)

bd_db2 <- phenotype2 %>%inner_join (bd_temp, c("id"="ID"))
total_population_number <- dim(bd_db2)

#bd_db_male<-bd_db2[bd_db2$sex=="Male",]
#bd_db_female<-bd_db_ace[bd_db_ace$sex=="Female",]

genotype_name <- colnames(bd_db2)[1861:1867]
phenotype_name <- colnames(bd_db2)[2:1845]
covariate_name <- colnames(bd_db2)[c(1846:1847,1849:1860)]

save(bd_db2,covariate_name,phenotype_name,genotype_name,file="query_phecov.RData")

#####################
###########################################
library(PheWAS)
library(xlsx)
load("/scratch/jz64675/neald/query_phecov.RData")
casefilter = 20
query <- data.frame()

for (j in seq(1,7,3)) {
  if(j < 7) {k = j + 2}
   else {k = k+1}
query_PheWAS_R <- phewas_ext(bd_db2,phenotypes=phenotype_name,genotypes=genotype_name[j:k],covariates = covariate_name, cores=4)
query_PheWAS_R <- addPhecodeInfo(query_PheWAS_R)
query <- rbind(query,query_PheWAS_R)

save(query_PheWAS_R, file = paste0("query2PhewasR", j,".RData"))

snp_table <- names(table(query_PheWAS_R$snp))
query.Result <- list()

# In order to avoid memory exhausting
for (i in 1:length(snp_table)){
  query.Result[[i]] <- query_PheWAS_R[which(query_PheWAS_R$snp == snp_table[i]),]
  query.Result[[i]] <- query.Result[[i]][which(query.Result[[i]]$n_cases >=casefilter),]
  query.Result[[i]]$padjust <- p.adjust(query.Result[[i]]$p,"fdr",length(query.Result[[i]]$p))
  if(i == 1) {write.xlsx(query.Result[[i]],file = paste0("query_Phewas",j, ".xlsx"),sheetName = snp_table[i],row.names = FALSE)}
  else{write.xlsx(query.Result[[i]],file = paste0("query_Phewas",j,".xlsx"),sheetName = snp_table[i],append = TRUE,row.names = FALSE)}
}
print("finish one cycle")
}


####################
#PBS -S /bin/bash
#PBS -N query
#PBS -l nodes=1:ppn=10
#PBS -q ye_q
#PBS -l walltime=120:00:00
#PBS -l mem=20gb

cd /scratch/jz64675/neald

module load R/3.6.2-foss-2018a-X11-20180131-GACRC

R CMD BATCH 2phewas.R
