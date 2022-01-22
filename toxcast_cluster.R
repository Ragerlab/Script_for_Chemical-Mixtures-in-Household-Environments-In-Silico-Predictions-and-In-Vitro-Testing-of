rm(list=ls())

#load libraries
library(cluster)
library(vegan)
library(gplots)
library(factoextra)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(janitor)

#set working directory
setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/ToxCast_Household_Chemicals/Experiment_1")

#read in the chemical/ exposure source category presence/absence dataframe and DTXSID/ name mapping
chems <- read_csv("output/presence_absence_binary_df.csv")
name_ref <- read_csv("output/ChemExpodb_name_DTXSID_ref.csv")


#merge in chemical name by DTXSID to the presence/absence dataframe
chems<- merge(chems, name_ref, by="DTXSID", all.x = TRUE)

#make ToxCast chemicals status a boolean
chems$tc <- as.logical(chems$tc)

#make a dataframe of just ToxCast chemicals
just_tc <- chems %>% filter(tc==TRUE)

#pull the list of names and DTXSIDS of just the ToxCast chemicals
names <- just_tc$true_chemname
ids <- just_tc$DTXSID


###########################  Chemical Clustering  ####################################################

#prepare dataframe of ToxCast chemicals for clustering 
chems_reduced <- just_tc %>% select(!c(DTXSID,count,tc, true_chemname))
 
#create a distance matrix of the reduced dataframe to cluster ToxCast chemicals
D_chems <- vegdist(as.matrix(chems_reduced),method = "jaccard")
 

#determine the optimal number of clusters of chemicals.
fviz_nbclust(chems_reduced, diss = D_chems, method = "wss", FUN=hcut, hc_func="diana", k.max=34)
fviz_nbclust(chems_reduced, diss = D_chems, method = "silhouette", FUN=hcut, hc_func="diana", k.max=34)
 

#cluster and look at dendrogram
cluster_chems<-diana(D_chems, diss=TRUE)
pltree(cluster_chems, labels = names_tc)
 
 
#get the cluster assignments for each chemical considering 17 clusters
ncut <-17
cluster_assignments_chems <- cutree(cluster_chems, k = ncut)
 
#number of chems in each cluster
k_chems<-table(cluster_assignments_chems)
 
#check to make sure everything is still in the same order. It is.
identical(just_tc %>% select(!c(DTXSID, count, tc,true_chemname)), chems_reduced)
 
#now we can reintroduce DTXSIDs and chemical names as well as assign the cluster number. Arrange the df by cluster number
#and add a numerical index column for identifying where cluster groupings and creating cluster separators in the heatmap later on
chems_reduced$DTXSID <- ids
chems_reduced$true_chemname <- names
chems_reduced$cluster<-cluster_assignments_chems
chems_reduced <- chems_reduced %>% arrange(cluster)
chems_reduced$index <- as.numeric(rownames(chems_reduced))
 

###########################  Exposure Source Category Clustering   ####################################################

#transpose initial presence/absence df so we can cluster the exposure source categories
cT <- just_tc %>% select(!c(count,tc, true_chemname))
cT <- t(cT)
cT <- as.data.frame(cT %>% row_to_names(row_number = 1))

#Not all exposure source categories are reflected in this dataset so remove the ones that are absent for all (ToxCast) chemicals
cT <- cT %>% mutate_all(as.character) %>% mutate_all(as.numeric) %>% rownames_to_column("DTXSID") %>% adorn_totals(where = "col")
cT <- cT %>% column_to_rownames("DTXSID")
cT <- cT %>% select(c("Total",colnames(cT)[1:ncol(cT)-1])) %>% filter(Total>0) %>% select(!Total)
escs <- rownames(cT)

#create distance matrix from the filtered transposed df to cluster the exposure source categories
D_esc <- vegdist(as.matrix(cT),method = "jaccard")


# look for optimal number of clusters
fviz_nbclust(cT, diss = D_esc, method = "wss", FUN=hcut, hc_func="diana", k.max=29)
fviz_nbclust(cT, diss = D_esc, method = "silhouette", FUN=hcut, hc_func="diana", k.max=29)

#cluster and look at dendrogram
cluster_escs<-diana(D_esc, diss=TRUE)
plot(cluster_escs)


#let's go with 13 clusters of exposure source categories
ncut <- 13
cluster_assignments_escs <- cutree(cluster_escs, k = ncut)

#number of exposure source categories in each cluster
k_escs<-table(cluster_assignments_escs)

#Arrange the df by cluster number and add a numerical index column for creating cluster separators later on
cT$esc <- escs
cT$cluster <- cluster_assignments_escs

cT <- cT %>% arrange(cluster)

cT$index <- 1:nrow(cT)

new_col_order <- cT$esc

#set the appropriate order for the exposure source categories in the heatmap, based on the clustering results, while
#also including the chemical identifier information
new_col_order <- c(c("DTXSID", "true_chemname", "cluster", "index"),new_col_order)


###################################### Make Heatmap #################################################################

#arrange the filtered ToxCast presence absence data in the correct order based on the exposure source category column order we just set
chems_reduced <- chems_reduced[new_col_order]

#write out the organized, clustered presence/absence dataframe
clust_df <- as.data.frame(chems_reduced)
clust_df <- clust_df %>% select(!index)
write.csv(clust_df,"output/ToxCast_chems_cluster_assignments.csv", row.names = FALSE) 
 
#make row and column separators so we can split the heatmap into the appropriate clusters
seprows <- chems_reduced %>% group_by(cluster) %>% slice_max(n=1, order_by=index)
seprows <- seprows$index
 
sepcols <- cT %>% group_by(cluster) %>% slice_max(n=1, order_by=index)
sepcols <- sepcols$index
 
#get chemical names to use as labels on heatmap
names <- chems_reduced$true_chemname

#remove extraneous columns 
chems_reduced <- chems_reduced %>% select(!c(DTXSID,cluster, true_chemname, index))
 
#make heatmap
pheatmap(as.matrix(chems_reduced), main="148 ToxCast Chemicals",
          cluster_rows=FALSE, cluster_cols = FALSE, show_rownames = TRUE,
          legend = FALSE, cellheight = 7, cellwidth = 18, fontsize_co1 = 5, fontsize_row = 5,  angle_col = 45, gaps_row = seprows, gaps_col = sepcols,
          labels_row=names,filename = "ToxCast_chems_clustered_heatmap.png", height = 20, width = 16)
 
 
 

 
 
 
 
 