#Load packages
import pandas as pd
import numpy as np
import os

#set working directory
os.chdir("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/ToxCast_Household_Chemicals/Experiment_1")


#read in CPCat DTXSID/keyword sets
kw_sets=pd.read_csv("input/factotum_listPresence_092320_updated_chemnames_020621.csv", dtype=str)

#make reference table of CPCat DTXSID/true chemical name pairs
name_dtxsid=kw_sets[["DTXSID","true_chemname"]]
name_dtxsid.drop_duplicates(inplace=True)
name_dtxsid.to_csv("output/name_dtxsid_ref.csv", index=False)


#read in keyword to exposure source category mapping
esc=pd.read_csv("input/keyword_esc.csv", usecols=[0,1,4])

#list of keywords that, when present in a set, warrant the chemical record being dropped from the CPCat DTXSID/keyword dataset
kw_1=list(esc.loc[esc.Status==1].Keyword)
kw_1_search=list(x.split()[0] for x in kw_1)

#list of keywords that, when present in a set, warrant removal from the set, but the other keywords remain in the CPCat  DTXSID/keyword dataset
kw_2=list(esc.loc[esc.Status==2].Keyword)

#dictionary of keywords to consider and the exposure source category they map to for analysis
kw2esc=esc.loc[esc.Status==0]
kw2esc=dict(zip(kw2esc.Keyword,kw2esc.Exposure_Source_Cat))

#read in relevant ToxCast data
toxcast=pd.read_excel("input/ToxCast_PPARg_DataPull_050321.xlsx", sheet_name="PositiveChemicals")
toxcast=toxcast[["ChemicalName","CASRN","DSSToxID"]]

#make list of ToxCast chemical DTXSIDs
toxcast_dtxsid=list(set(toxcast.DSSToxID))

#KW Filter 1: Remove CPCat records that don't have a kw set
kw_sets=kw_sets.loc[pd.isnull(kw_sets.listSets)==False]
kw_sets=kw_sets.loc[pd.isnull(kw_sets.true_chemname)==False]
kw_sets=kw_sets[["DTXSID","listSets"]]
kw_sets.reset_index(inplace=True, drop=True)

#KW Filter 2: Remove CPCat records that contain keywords that are found in 'kw_1'
kw_sets=kw_sets.loc[kw_sets.listSets.str.contains("|".join(kw_1_search))==False]


#KW Filter 3: Remove instances of keywords found in 'kw_2', but retain the record and other keywords in the CPCat dataset
kw_sets.loc[kw_sets.listSets.str.contains(";"),["listSets"]]=kw_sets.listSets.str.split(";")
kw_sets=kw_sets.rename(columns={"listSets":"keyword"}).explode("keyword")
kw_sets.keyword=kw_sets.keyword.str.strip()
kw_sets=kw_sets.loc[kw_sets.keyword.isin(kw_2)==False]


#Map keywords to exposure source categories and remove duplicate associations of CPCat DTXSIDS and exposure source categories
chem_esc=kw_sets.copy()
chem_esc["esc"]=chem_esc.keyword.replace(kw2esc)

chem_esc=chem_esc[["DTXSID","esc"]]
chem_esc.drop_duplicates(inplace=True)


#make a binary dataframe showing whether any DTXSID in the CPCat dataset is (1), or is not (0) associated with an exposure source category
pres_abs=pd.pivot_table(chem_esc,index="DTXSID",columns="esc",aggfunc=np.count_nonzero, fill_value=0).astype(bool)*1
pres_abs.reset_index(inplace=True)
pres_abs.columns.name=None
pres_abs.columns=[x.strip(" ") for x in pres_abs.columns]


#add count column with number of associated exposure source categories for each DTXSID and select DTXSIDs with at least 2 associated exposure source categories
pres_abs["count"]=np.sum(pres_abs,axis=1)
pres_abs= pres_abs.loc[pres_abs["count"]>=2]


#Identify which DTXSIDs in the CPCat dataset are also in the ToxCast dataset
pres_abs["tc"]=False
pres_abs.loc[pres_abs.DTXSID.isin(toxcast_dtxsid), ["tc"]]=True


#write the presence/absence df of all CPCat DTXSIDs that map to at least two exposure, along with the there ToxCast status, to a csv.
pres_abs.to_csv("output/presence_absence_binary_df.csv", index=False)
