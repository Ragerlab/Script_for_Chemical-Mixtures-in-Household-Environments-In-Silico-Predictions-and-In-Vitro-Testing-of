library(cluster)
library(vegan)
library(gplots)
library(factoextra)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

rm(list=ls())

#read in binary dataframe and name matching file
chems <- read_csv("derived_data/presence_absence_binary_df_DTXSID.csv")
name_ref <- read_csv("derived_data/ChemExpodb_name_DTXSID_ref.csv")

#merge in chemical name by DTXSID
chems<- merge(chems, name_ref, by="DTXSID", all.x = TRUE)

#make reference lists of the DTXSID and the chemical names for the chems dataframe
ids <- chems$DTXSID
names <- chems$true_chemname

#isoloate the the presence/absence data
chems_reduced <- chems %>% select(!c(DTXSID,count, true_chemname))

# create a distance matrix of the reduced dataframe
D <- vegdist(as.matrix(chems_reduced),method = "jaccard")

#make a dataframe of the matrix and write it to a csv
D_df <- data.frame(as.matrix(D))

colnames(D_df) <- ids
rownames(D_df) <- ids

write.csv(D_df, "derived_data/distance_matrix.csv")

#melt the matrix dataframe, fitler it to omit self correlated values, and write it to a csv
D_df_melt <- D_df %>% mutate(DTXSID_1=rownames(D_df)) %>% as_tibble() %>% 
  pivot_longer(!DTXSID_1, names_to = "DTXSID_2", values_to = "distance")%>% 
  left_join(name_ref, by=c("DTXSID_1"="DTXSID"))%>% 
  left_join(name_ref, by=c("DTXSID_2"="DTXSID"), suffix = c("_1","_2")) %>% 
  select(DTXSID_1, true_chemname_1, DTXSID_2, true_chemname_2, distance) %>% 
  filter(DTXSID_1!=DTXSID_2) %>% 
  arrange(distance)


write_csv(D_df_melt, "derived_data/melted_filtered_distance_matrix.csv")




dev.off()

#determine the optimal number of clusters. Looks like 20 could be reasonable.
fviz_nbclust(chems_reduced, diss = D, method = "wss", FUN=hcut, hc_func="diana", k.max=34)

fviz_nbclust(chems_reduced, diss = D, method = "silhouette", FUN=hcut, hc_func="diana", k.max=34)

#cluster and visualize the data
cluster2<-diana(D, diss=TRUE)
pltree(cluster2, labels = names)
# plot(cluster2, labels=names)


#get the assignments assuming 20 clusters
ncut <- 20
cluster_assignments <- cutree(cluster2, k = ncut)

#number of chems in each cluster
k2<-table(cluster_assignments)


#check to make sure everything is still in the same order. It is.
identical(chems %>% select(!c(DTXSID, count, true_chemname)), chems_reduced)

#Now we can reintroduce DTXSIDs and chemical names as well as assign the cluster number.
chems_reduced$DTXSID <- ids
chems_reduced$true_chemname <- names
chems_reduced$cluster<-cluster_assignments
chems_reduced <- chems_reduced %>% arrange(cluster)


#make a new reference of chemical names since we have rearranged the dataframe order by cluster
names <- chems_reduced$true_chemname

#make the row identifier the DTXSID
rownames(chems_reduced) <- chems_reduced$DTXSID


#Again, reduce the dataframe back down to the presence/absence data
chems_reduced <- chems_reduced %>% select(!c(DTXSID,cluster, true_chemname))

#make a heatmap using the distance matrix. Note clustering is off for both rows and columns as this independently
#reclusters the data.
pheatmap(as.matrix(chems_reduced), main="Chems", clustering_distance_rows = D, 
         cluster_rows=FALSE, cluster_cols = FALSE, labels_row = names, legend = TRUE, 
         cellheight = 15, cellwidth = 15, fontsize_co1 = 5,  angle_col = 45,  
         filename = "figures/toxcast_priority_cluster.png", height = 10, width = 16)





