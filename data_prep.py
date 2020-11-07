import pandas as pd
import os
import numpy as np


#read in ToxCast priority chemicals
toxcast=pd.read_excel("inputs//TopChemicals_102820.xlsx")

#read in ChemExpodb list presence chemical keyword sets
CEdb=pd.read_csv("inputs//all_factotum_DTXSID_listPresence_2020-09-23.csv", dtype=str)

#get list of overlapping DTXSIDs
dtxsid_both=list(set(toxcast.DSSToxID) & set(CEdb.DTXSID))

#Only keep ChemExpodb records for overlapping DTXSIDs
CEdb=CEdb.loc[CEdb.DTXSID.isin(dtxsid_both)]

#reference of chemical name DTXSID mapping from ChemExpodb
name_dtxsid=CEdb[["DTXSID","true_chemname"]]
name_dtxsid.drop_duplicates(inplace=True)
name_dtxsid.to_csv("derived_data//ChemExpodb_name_DTXSID_ref.csv", index=False)


#read in the keyword/bin mapping
bins=pd.read_csv("derived_data//keyword_bins.csv", usecols=[0,1,4])

#list of keywords that, when present in a set, warrant the chemical record being dropped
kw_1=list(bins.loc[bins.Status==1].Keyword)
kw_1_search=list(x.split()[0] for x in kw_1)

#list of keywords that, when present in a set, warrant removal from the set, but the other keywords remain.
kw_2=list(bins.loc[bins.Status==2].Keyword)

#dictionary of keywords to consider and the bin they map to for analysis
kw2bin=bins.loc[bins.Status==0]
kw2bin=dict(zip(kw2bin.Keyword,kw2bin.Master_Bin))

"""
KW Filter 1: Remove ChemExpodb chemicals that don't have a kw set

"""

CEdb=CEdb.loc[pd.isnull(CEdb.listSets)==False]
CEdb=CEdb[["DTXSID","listSets"]]
CEdb.reset_index(inplace=True, drop=True)

"""
KW Filter 2: Remove ChemExpodb chemical records that contain keywords in kw_1.

"""
CEdb=CEdb.loc[CEdb.listSets.str.contains("|".join(kw_1_search))==False]



"""
KW Filter 3: Remove instances of keywords in kw_2

"""

CEdb.loc[CEdb.listSets.str.contains(";"),["listSets"]]=CEdb.listSets.str.split(";")
CEdb=CEdb.rename(columns={"listSets":"keyword"}).explode("keyword")
CEdb.keyword=CEdb.keyword.str.strip()
CEdb=CEdb.loc[CEdb.keyword.isin(kw_2)==False]


"""
Map keywords to bins.

"""

CEdb["bin"]=CEdb.keyword.replace(kw2bin)

final=CEdb[["DTXSID","bin"]]
final.drop_duplicates(inplace=True)





"""
Create binary df of DTXSID and bin based on whether the DTXSID is associated with the bin or not

"""

presence_df=pd.pivot_table(final,index="DTXSID",columns="bin",aggfunc=np.count_nonzero).fillna(0).astype(bool)*1
presence_df.reset_index(inplace=True)
presence_df.columns.name=None
presence_df.columns=[x.strip(" ") for x in presence_df.columns]


#add count column with number of used keywords just as a reference
presence_df["count"]=np.sum(presence_df,axis=1)



"""
write df to csv.
"""

presence_df.to_csv("derived_data//presence_absence_binary_df_DTXSID.csv", index=False)
