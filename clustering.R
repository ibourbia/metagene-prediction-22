library(tidyverse)
library(dplyr)
library(ggplot2)
library(dendextend)
library(psych)
library(Hmisc)
library(caret)
library(ggfortify)
library(miceadds)
library(finalfit)
library(broom)
library(OddsPlotty)
library(survival)
library(survminer)
library(enrichR)

# Clean variables #
#rm(list=ls())

# Set Seed so that same sample can be reproduced in future also #
set.seed(101) 


# Correlations between genes expression
# Might saturate memory. If so a restart is in order.
discovery_generic <- exp_mat_generic[,sample(ncol(exp_mat_generic), ncol(exp_mat_generic)/2) ]
confirmation_generic <- exp_mat_generic[, -which(names(exp_mat_generic) %in% names(discovery_generic))]
discovery_generic <- as.data.frame(t(discovery_generic))
confirmation_generic <- as.data.frame(t(confirmation_generic))


correlations <- cor(discovery_generic)
#load("data\\correlation.RData")
distances <- as.dist(1-abs(correlations))
correlation_conf <- cor(confirmation_generic) 
#load("correlation_conf.RData")
distances_conf <- as.dist(1-abs(correlation_conf))

#In case of memory issues
#save(correlations,file="data\\correlation.RData")
#save(correlation_conf,file="data\\correlation_conf.RData")
#Attempt to ease the sourcing
remove(correlations,correlation_conf)
gc()
# Hierarchical clustering discovery
hc <- hclust(distances,method="average")
#Cluster cutting at tree height = 0.6 and minimal cluster size of 25 (according to paper)
cluster_labels <- cutree(hc,h=0.6)
cut <- data.frame(cluster_labels)

#Keeping only clusters with a specific size
clusters_discovery <- data.frame(cluster_sizing(hc))
colnames(clusters_discovery)<- c("cluster_labels")
clusters_discovery$gene<-rownames(clusters_discovery)

#Hierarchical clustering in confirmation set
hc_conf <- hclust(distances_conf,method="average")
cluster_labels_conf <- cutree(hc_conf,h=0.6)
cut_conf <- data.frame(cluster_labels_conf)

#Sizing
clusters_conf <-data.frame(cluster_sizing(hc_conf))
colnames(clusters_conf)<- c("cluster_labels")
clusters_conf$gene<-rownames(clusters_conf)

#Cluster matching between discovery and confirmation sets

cluster_map <- {}
for (i in 1:nrow(clusters_discovery)) {
  cluster <- clusters_discovery$cluster_labels[i]
  if (cluster %in% clusters_conf$cluster_labels) {
    cluster_map[cluster] <- clusters_conf$cluster_labels[which(clusters_conf$cluster_labels == cluster)]
  } else {
    cluster_map[cluster] <- NA
  }
}

cluster_map<- cluster_map[!is.na(cluster_map)]
confirmed_clusters<- cut %>% filter(cluster_labels %in% cluster_map)%>% as.data.frame()
confirmed_clusters$gene <- rownames(confirmed_clusters)
####################t#############################################
# Identification of clusters enriched in proliferation and      #
# ER-related genes                                              #
# Could have been done with EnrichR,but I could not connect to  #
# the database whatsoever                                       #
# Update 2023 : Connexion worked, updating the code             #
#################################################################

#ER related metagene is cluster containing ESR1 and related genes

if(!("ESR1"%in% row.names(confirmed_clusters))){
  stop("ER-related cluster does not seem to exist") #worth stopping the sourcing to inspect
}else{
ER_cluster<- confirmed_clusters %>% filter(cluster_labels==confirmed_clusters["ESR1",]$cluster_labels)
}

# gene_list <- confirmed_clusters%>% rownames()%>%as.list()
# 
# for(i in unique(confirmed_clusters$cluster_labels)){
#   metagenes <- confirmed_clusters%>%
#     filter(cluster_labels==i)%>%
#     rownames()%>%
#     as.vector()
#   write.table(metagenes,paste("data\\cluster_",i,".txt",sep=""),row.names=FALSE,sep="\t", quote = FALSE)
#   print(metagenes)
#   enriched_gene_sets <- enrichr( metagenes,
#                                 databases = c("GO_BP", "GO_MF", "KEGG"))
# }
# EnrichR does not find anything -> manual check up. Cluster 29 is proliferation
# Associated with MYC 
proliferation_cluster <- confirmed_clusters %>% filter(cluster_labels==29)
#genes_discovery <- arrange(filter(cut, cluster_labels %in% cutted_clusters),desc(cluster_labels))
#genes_conf <- arrange(filter(cut_conf, cluster_labels_conf %in% cutted_clusters_conf),desc(cluster_labels_conf))
rows_29 <- filter(genes_conf,cluster_labels_conf == 29) #Proliferation cluster
rows_29 <- cbind(rownames(rows_29))
rows_72 <- filter(genes_conf,cluster_labels_conf == 72) #ER_related cluster
rows_72 <- cbind(rownames(rows_72))

write.csv(p_values,"data\\p_values.csv")
write.csv(proliferation,"data\\proliferation.csv")
write.csv(exp_mat_prog,"data\\exp_mat_prog.csv")
write.csv(rows_29,"data\\cluster29.csv")
write.csv(pheno_prog,"data\\pheno_prog.csv")
write.csv(exp_mat_29,"data\\exp_mat29.csv")


########################################
# Cross validation                     #
########################################

#take 90 % of samples, check how good genes predict the outcome, and take top 10
#Supplementary data no 6c
cross_validation<-setNames(data.frame(matrix(ncol = 4, nrow = 100)), c("Repeat","pvalues","Metagene_size","Genes"))
exp_mat_proliferation <- exp_mat_prog[proliferation_cluster$gene, which(names(exp_mat_prog) %in% rownames(proliferation))]

for (r in seq(1,100)){
  #Training set and test set
  training_set <- exp_mat_proliferation[,sample(ncol(exp_mat_proliferation), ncol(exp_mat_proliferation)*0.9)]
  test_set<- exp_mat_proliferation[proliferation_cluster$gene, -which(names(exp_mat_proliferation) %in% names(training_set))]
  test_set <- data.frame(t(test_set))
  
  test_set <- merge(test_set,pheno_prog[,c("t.dmfs.norm.years","event","time")],by=0)
  training_set <- data.frame(t(training_set))
  training_set1 <- merge(training_set,pheno_prog[,c("t.dmfs.norm.years","event","time")],by=0)
  
  #training_set <- exp_mat_prog[rows_29,sample(ncol(exp_mat_prog), ncol(exp_mat_prog)*0.9) ]
  
  #test_set<- exp_mat_prog[rows_29, -which(names(exp_mat_prog) %in% names(training_set))]
  #test_set <- data.frame(t(test_set))
  #test_set <- merge(test_set,pheno_prog[,c("t.dmfs.norm.years","event","time")],by=0)
  #test_set <- data.frame(t(test_set))
  
  #training_set <- data.frame(t(training_set))
#training_set1 <- merge(training_set,pheno_prog[,c("t.dmfs.norm.years","event","time")],by=0)

  #Cox model fitting
  p_values <- data.frame(row.names=proliferation_cluster$gene)
  p_values$pvalue <- NA
  for(gene in proliferation_cluster$gene){
  cox <- coxph(Surv(training_set1$time,training_set1$event)~training_set1[,gene])
  p_values[gene,]<- summary(cox)$sctest[3]
  #summary(cox)$coefficients[, 5]
}

p_values <- p_values %>% arrange(desc(pvalue))
sizes <-round(seq(10,length(proliferation_cluster$gene),length.out=17))
p_values_sizes <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(p_values_sizes) <- c("size","pvalue")

#calculate row means of genes 
for(i in seq(1,17)){
  test_set1<- test_set
  metagene <-tail(p_values,sizes[i])
  test_set1<- as.data.frame(test_set1[,c(rownames(metagene),"t.dmfs.norm.years","event","time")])
  test_set1$mean <- test_set1%>% rowMeans()
  tert <- quantile(test_set1$mean, c(0:3/3))
  test_set1$tertiles <- with(test_set1, cut(mean, tert, include.lowest = T, labels = c(0, 1, 2)))
  cox_tert <- survdiff(Surv(test_set1$time,test_set1$event)~test_set1$tertiles)
  print(cox_tert)
  p_values_sizes<-p_values_sizes%>%add_row(size=sizes[i],pvalue=1 - pchisq(cox_tert$chisq, length(cox_tert$n) - 1))
}
p_values_sizes <- p_values_sizes %>% arrange(desc(p_values_sizes$pvalue))
genes <-rownames(tail(p_values,p_values_sizes[which.min(p_values_sizes$pvalue),]$size))
row <- c(r,p_values_sizes[which.min(p_values_sizes$pvalue),]$pvalue,length(genes),toString(genes))
cross_validation[r,]<- row
#cross_validation <- cross_validation %>% add_row(Repeat=r,pvalues=p_values_sizes[which.min(p_values_sizes$pvalue),]$pvalue,Genes=toString(rownames(tail(p_values,p_values_sizes[which.min(p_values_sizes$pvalue),]$size))))
#cross_validation[,"Genes"]<-toString(rownames(tail(p_values,p_values_sizes[which.min(p_values_sizes$pvalue),]$size)))

}
#DMFS % = proportion that survived in each tertile


metagene_proliferation <- strsplit(cross_validation[which.min(cross_validation$pvalues),]$Genes,split=", ")
exp_mat_proliferation <- data.frame(t(exp_mat_proliferation))
exp_mat_proliferation <- merge(exp_mat_proliferation,proliferation[,c("t.dmfs.norm.years","event","time")],by=0)


table10 <- as.data.frame(exp_mat_proliferation[,c(metagene_proliferation[[1]],"t.dmfs.norm.years","event","time")])
table10$mean <- table10%>% rowMeans()
tert <- quantile(table10$mean, c(0:3/3))
table10$tertiles <- with(table10, cut(mean, tert, include.lowest = T, labels = c(0, 1, 2)))

low_tertile <- table10[table10$tertiles==0,]
mid_tertile <-table10[table10$tertiles==1,]
high_tertile <- table10[table10$tertiles==2,]
which(table10$tertiles==1)
table(low_tertile[,"t.dmfs.norm.years"]>=5)
table(mid_tertile[,"t.dmfs.norm.years"]>=5)
table(high_tertile[,"t.dmfs.norm.years"]>=5)

#comparison : logistic regression , independent variables = tertiles (as levels), outcome= PR, low=baselin

#For TAM : 
exp_mat_midhigh <- data.frame(t(exp_mat_tam))
exp_mat_midhigh <- exp_mat_midhigh[,which(names(exp_mat_midhigh) %in% metagene_proliferation[[1]])]
exp_mat_midhigh$mean <- exp_mat_midhigh%>%rowMeans()
tert <- quantile(exp_mat_midhigh$mean, c(0:3/3))
exp_mat_midhigh$tertiles <- with(exp_mat_midhigh, cut(mean, tert, include.lowest = T, labels = c(0, 1, 2)))
midhigh_proliferation_tam <- rownames(exp_mat_midhigh[exp_mat_midhigh$tertiles==1|exp_mat_midhigh$tertiles==2,])
exp_mat_midhigh <- exp_mat_tam[,midhigh_proliferation_tam]


cross_validation_tam<-setNames(data.frame(matrix(ncol = 4, nrow = 100)), c("Repeat","pvalues","Metagene_size","Genes"))
for (r in seq(1,100)){
  training_set <- exp_mat_midhigh[ER_cluster$gene,sample(ncol(exp_mat_midhigh), ncol(exp_mat_midhigh)*0.9)]
  test_set<- exp_mat_midhigh[ER_cluster$gene, -which(names(exp_mat_midhigh) %in% names(training_set))]
  test_set <- data.frame(t(test_set))
  
  test_set <- merge(test_set,pheno_tam[,c("t.dmfs.norm.years","event","time")],by=0)
  test_set <- na.omit(test_set)
  
  training_set <- data.frame(t(training_set))
  training_set1 <- merge(training_set,pheno_tam[,c("t.dmfs.norm.years","event","time")],by=0)
  
  
  p_values <- data.frame(row.names=ER_cluster$gene)
  p_values$pvalue <- NA
  for(gene in ER_cluster$gene){
    cox <- coxph(Surv(training_set1$time,training_set1$event)~training_set1[,gene])
    p_values[gene,]<- summary(cox)$sctest[3]#summary(cox)$coefficients[, 5]
  }
  
  p_values <- p_values %>% arrange(desc(pvalue))
  sizes <-round(seq(10,length(ER_cluster$gene),length.out=17))
  p_values_sizes <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(p_values_sizes) <- c("size","pvalue")
  
  #calculate row means of genes 
  for(i in seq(1,17)){
    test_set1<- test_set
    metagene <-tail(p_values,sizes[i])
    test_set1<- as.data.frame(test_set1[,c(rownames(metagene),"t.dmfs.norm.years","event","time")])
    test_set1$mean <- test_set1[,rownames(metagene),]%>% rowMeans()
    tert <- quantile(test_set1$mean, c(0:3/3))
    test_set1$tertiles <- with(test_set1, cut(mean, tert, include.lowest = T, labels = c(0, 1, 2)))
    cox_tert <- survdiff(Surv(test_set1$time,test_set1$event)~test_set1$tertiles)
    print(cox_tert)
    p_values_sizes<-p_values_sizes%>%add_row(size=sizes[i],pvalue=1 - pchisq(cox_tert$chisq, length(cox_tert$n) - 1))
  }
  p_values_sizes <- p_values_sizes %>% arrange(desc(p_values_sizes$pvalue))
  genes <-rownames(tail(p_values,p_values_sizes[which.min(p_values_sizes$pvalue),]$size))
  row <- c(r,p_values_sizes[which.min(p_values_sizes$pvalue),]$pvalue,length(genes),toString(genes))
  cross_validation_tam[r,]<- row
  #cross_validation <- cross_validation %>% add_row(Repeat=r,pvalues=p_values_sizes[which.min(p_values_sizes$pvalue),]$pvalue,Genes=toString(rownames(tail(p_values,p_values_sizes[which.min(p_values_sizes$pvalue),]$size))))
  #cross_validation[,"Genes"]<-toString(rownames(tail(p_values,p_values_sizes[which.min(p_values_sizes$pvalue),]$size)))
}
metagene_er_related <- strsplit(cross_validation_tam[which.min(cross_validation_tam$pvalues),]$Genes,split=", ")[[1]]
exp_mat_midhigh <- t(exp_mat_midhigh)
exp_mat_midhigh <- merge(exp_mat_midhigh,pheno_tam[,c("t.dmfs.norm.years","event","time")],by=0)
table15 <- na.omit(as.data.frame(exp_mat_midhigh[,c(metagene_er_related,"t.dmfs.norm.years","event","time")]))
table15$mean <- table15%>% rowMeans()
tert <- quantile(table15$mean, c(0:3/3))
table15$tertiles <- with(table15, cut(mean, tert, include.lowest = T, labels = c(0, 1, 2)))

low_tertile <- table15[table15$tertiles==0,]
mid_tertile <-table15[table15$tertiles==1,]
high_tertile <- table15[table15$tertiles==2,]
which(table15$tertiles==1)
table(low_tertile[,"t.dmfs.norm.years"]>=5)
table(mid_tertile[,"t.dmfs.norm.years"]>=5)
table(high_tertile[,"t.dmfs.norm.years"]>=5)

chemo_horak <- (pheno_chemo[pheno_chemo$Data_Set_Name=="Horak" & pheno_chemo$ER == 1 & pheno_chemo$HER2 == 0 & !is.na(pheno_chemo$pCR),])
exp_mat_horak <- data.frame(t(exp_mat_chemo[c(metagene_er_related,metagene_proliferation[[1]]),which(names(exp_mat_chemo) %in% rownames(chemo_horak))]))
exp_mat_horak <- tertiles(exp_mat_horak)
chemo_horak$tertiles <- exp_mat_horak$tertiles
factor(chemo_horak$tertiles)
model <- glm(pCR~tertiles,data=chemo_horak,family=binomial(link="logit"))
confint(model)
exp(coef(model))

chemo_hatzis<- pheno_chemo[ (pheno_chemo$Data_Set_Name == "Hatzis_val" | pheno_chemo$Data_Set_Name == "Hatzis_disc") & pheno_chemo$ER == 1 & pheno_chemo$HER2 == 0 & !is.na(pheno_chemo$pCR),]
chemo_hatzis <- chemo_hatzis[!is.na(chemo_hatzis$pCR),]
exp_mat_hatzis_prolif<- data.frame(t(exp_mat_chemo[c(metagene_proliferation[[1]]),which(names(exp_mat_chemo) %in% rownames(chemo_hatzis))]))
exp_mat_hatzis_prolif <- highlow(exp_mat_hatzis_prolif)
chemo_hatzis$tertilesPRO <- exp_mat_hatzis_prolif$tertiles

exp_mat_hatzis_er <-data.frame(t(exp_mat_chemo[metagene_er_related,which(names(exp_mat_chemo) %in% rownames(chemo_hatzis))]))
exp_mat_hatzis_er <- highlow(exp_mat_hatzis_er)
chemo_hatzis$tertilesER<- exp_mat_hatzis_er$tertiles

chemo_hatzis$tertiles <- NA
chemo_hatzis <- chemo_hatzis %>% 
  mutate(tertiles= case_when(
    tertilesER == 1 & tertilesPRO==0 ~'Low',
    tertilesER==0 & tertilesPRO==0 ~'Intermediate',
    tertilesER==1 & tertilesPRO==1 ~'Intermediate',
    tertilesER==0 & tertilesPRO==1 ~'High'
  )
)

chemo_petel  <- pheno_chemo[pheno_chemo$Data_Set_Name == "Petel" & pheno_chemo$ER == 1 & pheno_chemo$HER2 == 0 ,]
chemo_petel <- chemo_petel[!is.na(chemo_petel$StudyType),]
exp_mat_petel_prolif<- data.frame(t(exp_mat_chemo[c(metagene_proliferation[[1]]),which(names(exp_mat_chemo) %in% rownames(chemo_petel))]))
exp_mat_petel_prolif <- highlow(exp_mat_petel_prolif)
chemo_petel$tertilesPRO <- exp_mat_petel_prolif$tertiles

exp_mat_petel_er <-data.frame(t(exp_mat_chemo[metagene_er_related,which(names(exp_mat_chemo) %in% rownames(chemo_petel))]))
exp_mat_petel_er <- highlow(exp_mat_petel_er)
chemo_petel$tertilesER<- exp_mat_petel_er$tertiles

chemo_petel$tertiles <- NA
chemo_petel <- chemo_petel %>% 
  mutate(tertiles= case_when(
    tertilesER == 1 & tertilesPRO==0 ~'Low',
    tertilesER==0 & tertilesPRO==0 ~'Intermediate',
    tertilesER==1 & tertilesPRO==1 ~'Intermediate',
    tertilesER==0 & tertilesPRO==1 ~'High'
  )
  )


chemo_hatzis$tertiles <- factor(chemo_hatzis$tertiles)
chemo_hatzis$tertiles <-relevel(chemo_hatzis$tertiles,ref = 'Low')
model <- glm(pCR~tertiles,data=chemo_hatzis,family=binomial(link="logit"))
OR_hatzis <- cbind(coef(model),confint(model))
confint(model)
coef(model)
xtable(exp(cbind("Odds ratio" = coef(model), confint.default(model, level = 0.95))))

kaplan <- rbind(chemo_hatzis,chemo_petel)
kaplan_meier <- survfit(Surv(time,event)~tertiles,data=kaplan)
ggsurvplot(kaplan_meier,pval=TRUE,data=kaplan,risk.table=TRUE)


chemo_er <- pheno_chemo[(pheno_chemo$Data_Set_Name=="Petel"|pheno_chemo$Data_Set_Name=="Hatzis_val"|pheno_chemo$Data_Set_Name=="Hatzis_disc")&
                        pheno_chemo$ER == 1 & pheno_chemo$HER2 == 0 & !is.na(pheno_chemo$t.dmfs.norm.years) & !is.na(pheno_chemo$node) & !is.na(pheno_chemo$grade),]
exp_mat_er_prolif<- data.frame(t(exp_mat_chemo[c(metagene_proliferation[[1]]),which(names(exp_mat_chemo) %in% rownames(chemo_er))]))
exp_mat_er_prolif <- highlow(exp_mat_er_prolif)
chemo_er$tertilesPRO <- exp_mat_er_prolif$tertiles

exp_mat_er_er <-data.frame(t(exp_mat_chemo[metagene_er_related,which(names(exp_mat_chemo) %in% rownames(chemo_er))]))
exp_mat_er_er <- highlow(exp_mat_er_er)
chemo_er$tertilesER<- exp_mat_er_er$tertiles

chemo_er$tertiles <- NA
chemo_er <- chemo_er %>% 
  mutate(tertiles= case_when(
    tertilesER == 1 & tertilesPRO==0 ~'Low',
    tertilesER==0 & tertilesPRO==0 ~'Intermediate',
    tertilesER==1 & tertilesPRO==1 ~'Intermediate',
    tertilesER==0 & tertilesPRO==1 ~'High'
  )
)

chemo_er$age_threshold <- NA
chemo_er <- chemo_er %>%
  mutate(age_threshold = case_when(
    age >50 ~ 1,
    age <= 50  ~ 0
  ))
chemo_er$tertiles<- factor(chemo_er$tertiles)
chemo_er$tertiles <- relevel(chemo_er$tertiles,ref='Low')
chemo_er$age_threshold<-factor(chemo_er$age_threshold)
chemo_er$node <- factor(chemo_er$node)
chemo_er <- chemo_er %>%
  mutate(grade= case_when(
    grade ==2 ~0,
    grade==3 ~1,
    grade==1~0
  ))
chemo_er$grade <- factor(chemo_er$grade)
multi_cox <- coxph(Surv(time,event)~tertiles+node+grade+age_threshold,data=chemo_er)
summary(multi_cox)
print(cbind('HR'=exp(coef(multi_cox)),confint.default(multi_cox)))
xtable(cbind('HR'=exp(coef(multi_cox)),summary(multi_cox)$conf.int))
      