## Script_for_Chemical-Mixtures-in-Household-Environments-In-Silico-Predictions-and-In-Vitro-Testing-of-Potential-Joint-Toxicities-in-Human-Liver-Cells
Script associated with the manuscript titled 'Chemical Mixtures in Household Environments: In Silico Predictions and In Vitro Testing of Potential Joint Toxicities in Human Liver Cells' in preparation for submission.

![Cluster5_HeatMapImage](https://user-images.githubusercontent.com/72747901/146393635-815c7716-b7f1-4052-9e00-f4a14a46e9bc.png)



## Script descriptions

**make_presence_absence.py**- Reads in the chemicals and associated keyword sets from the CPCat Chemical List Presence dataset **factotum_listPresence_092320_updated_chemnames_020621.csv**, the mapping of CPCat keywords to exposure source categories **keyword_esc.csv**, and ToxCast chemicals of interest from the **PositiveChemicals** sheet of **ToxCast_PPARg_DataPull_050321.xlsx**. Chemicals from both the CPCat and ToxCast datasets are identified as DTXSIDs. It is noted which of the chemicals contained in the CPCat dataset are also contained in the ToxCast dataset. Exposure source categories are mapped to all chemicals in the CPCat dataset based on associated keywords. Duplicate chemical/ exposure source category pairs are dropped and a filter is applied to only keep chemicals that are linked to at least two unique exposure source categories. An nxm dataframe is produced where n is the number of chemicals that pass the >=2 exposure source category filter, and m is the number of exposure source categories. A 1 in the dataframe indicates that there is an association between the chemical and the exposure source category, and 0 indicates there is not an association. The dataframe is written to the file **presence_absence_binary_df.csv**, which serves as one of the inputs to **toxcast_cluster.R**. Additionally, a reference file **name_dtxsid_ref.csv** is produced which maps the true_chemname to the DTXSID for each chemical, according to the original CPCat dataset. Both of these files are written to an *output* folder.

**toxcast_cluster.R**- reads in the binary dataframe **presence_absence_binary_df.csv** and the name/DTXSID mapping **name_dtxsid_ref.csv** produced by **make_presence_absence.py**. A dataframe of just the 148 chemicals from the CPCat dataset that were also identified as relevant ToxCast chemicals is produced. A distance matrix of the presence/absence data for these 148 chemical is generated by computing a jaccard distance for each combination of chemicals. Divisive hierarchical clustering is performed on the distance matrix to cluster the chemicals. A jaccard distance matrix is then generated for the exposure source categories and hierarchical clustering is employed to cluster exposure source categories. Clustering of the exposure source categories was performed to better organize the the final heatmap, and thereby aid in the interpretation of results. The clustered, organized presence/absence dataframe is written to **Toxcast_chems_cluster_assignments.csv**. A heatamp of this dataframe is then generated and saved as **Toxcast_chems_clustered_heatamp.png**. Additionally, a heatmap of just the chemicals in cluster 5, the cluster identified as being of greatest interest, is made and saved as **cluster_5_heatmap.png**
