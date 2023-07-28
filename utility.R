# Utilities #

# Compute tertitles for cox survival analysis
tertiles <- function(df) {
  df$mean <- df%>% rowMeans()
  tert <- quantile(df$mean, c(0:3/3))
  df$tertiles <- with(df, cut(mean, tert, include.lowest = T, labels = c(0, 1, 2)))
  return(df)
}

highlow <- function(df){
  df$mean <- df%>% rowMeans()
  tert <- quantile(df$mean, c(0:2/2))
  print(tert)
  df$tertiles <- with(df, cut(mean, tert, include.lowest = T, labels = c(0, 1)))
  return(df)
}

# Extract pheno data from assay

get_pheno_data <- function(df){
  pheno <- phenoData(df)
  pheno <- pheno@data
  return(pheno)
}

#Compute the ER and HER Status according to metagenes (given by paper). 
#Paper metagenes and cutoff :
#C6orf97, ESR1, EVL, ABAT, SLC39A6, GATA3,SCUBE2) cutoff= 7.5
#ERBB2, PGAP3, STARD3, GRB7, PNMT, PSMD3, GSDMB, RPL19, FGFR4, CAP1, cutoff= 8.35
ER_HER_status <- function(exp_mat) {
  #Idea : use R.oo to have pass by reference behavior for expr_means arguments
  
  #ER +/- prediction
  expr_means<- exp_mat %>% dplyr::filter(row.names(exp_mat) %in% c("CCDC170", "ESR1", "EVL", "ABAT", "SLC39A6", "GATA3","SCUBE2"))%>% colMeans()%>%as.data.frame
  colnames(expr_means)=c("ER_cutoff")
  #HER+/- prediction
  expr_means$HR_cutoff<-exp_mat %>% dplyr::filter(row.names(exp_mat) %in% c("ERBB2", "PGAP3", "STARD3", "GRB7", "PNMT","PSMD3", "GSDMB", "RPL19", "FGFR4", "CAP1"))%>% colMeans()%>%as.data.frame()
  
  # Predict status based on mean gene expression level cutoff
  #1 = +, 0 = -
  expr_means<-expr_means %>%
    mutate(ER_status = case_when(
      expr_means$ER_cutoff>7.5~ 1,
      expr_means$ER_cutoff<=7.5~ 0
    ))
  
  expr_means<-expr_means%>%
    mutate(HER2_status = case_when(
      expr_means$HR_cutoff>8.35~ 1,
      expr_means$HR_cutoff<=8.35~ 0
    ))
  
  # Subset of the mean expression level of genes based on ER & HER status  
  # Only interested in ER+HER- status 
  erherstatus<-subset(expr_means,expr_means$ER_status==1 & expr_means$HER2_status==0)
  
  return(expr_means)
}

# Keep only the clusters with a minimal size
cluster_sizing <- function(hc,size=25){
  cluster_labels <- cutree(hc,h=0.6)
  cluster_sizes <- table(cluster_labels) 
  sized_clusters <- which(cluster_sizes >= size)
  clusters <- cluster_labels[cluster_labels %in% sized_clusters]
  return(clusters)
}
