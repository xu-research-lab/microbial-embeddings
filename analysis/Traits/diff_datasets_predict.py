#！/bin/python

import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

from sklearn.metrics import roc_curve, auc 
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier  
from sklearn.model_selection import cross_val_score 
from sklearn.metrics import roc_auc_score, roc_curve

rf_model = RandomForestClassifier(n_estimators=1000, random_state=0, oob_score=True, n_jobs=-1, class_weight="balanced")
fid_8 = []
for i in [1,2,3,4,5]:
    co_embedding = pd.read_csv(f"xxx/embedding/subset_table_8w_100_{i}.txt",
                          header=None, sep=" ", low_memory=False, index_col=0)
    co_embedding.drop("<unk>", inplace=True)
    fid_8 = fid_8 + co_embedding.index.to_list()

auc = []
traits_type = []
group = []
tax = []
datasize = []

datasize_dict = {"1w": "10,000", "2w":"20,000", "4w":"40,000", "8w":"80,000", "16w":"160,000", "21w": "210,000"}
for n in [1,2,3,4,5]:
    for t in ["1w", "2w", "4w", "8w", "16w", "21w"]:
        if t != "21w":
            co_embedding = pd.read_csv(f"xxx/subset_table_{t}_100_{n}.txt",
                                       header=None, sep=" ", low_memory=False, index_col=0)
            co_embedding.drop("<unk>", inplace=True)
        else:
            co_embedding = pd.read_csv("../../data/social_niche_embedding_100.txt",
                          header=None, sep=" ", low_memory=False, index_col=0)
            co_embedding.drop("<unk>", inplace=True)
        fid_8 = np.unique(fid_8)
        inter_id = np.intersect1d(fid_8, co_embedding.index.values)
        co_embedding = co_embedding.loc[inter_id, ]
        fid = co_embedding.index.values
        taxonomy = pd.read_csv("../Pretraining_data_profile/Data/taxmap_slv_ssu_ref_nr_138.2.txt", sep="\t", low_memory=False)
        acc = []
        for i in range(taxonomy.shape[0]):
            temp = taxonomy.iloc[i]
            acc.append(f"{temp[0]}.{temp[1]}.{temp[2]}")
        
        taxonomy = taxonomy.loc[:, "path"].str.split(';', expand=True)
        taxonomy.index = acc
        taxonomy = taxonomy.iloc[:, 0: 7]
        taxonomy.columns = ["k", "p", "c", "o", "f", "g", "s"]
        taxonomy = taxonomy.loc[fid]
    
        traits =  pd.read_csv("Data/trait_predcit.csv", index_col=0)
        inter_id = np.intersect1d(fid, traits.index.values)
        traits = traits.loc[inter_id]
        traits = traits.astype(int)
        traits[traits.values == 3] = 1
        
        # Oxygen Preference
        # Initialize the column with 0
        traits['Oxygen_Preference'] = 0
        # Set values based on conditions
        traits.loc[traits['Aerobe'] == 1, 'Oxygen_Preference'] = "aerobic"
        traits.loc[traits['Facultative'] == 1, 'Oxygen_Preference'] = "facultatively"
        traits.loc[traits['Anaerobe'] == 1, 'Oxygen_Preference'] = "anaerobic"
        # Set to NA (Not a Number) if the row sum is not 1
        oxygen_cols = ['Aerobe', 'Facultative', 'Anaerobe']
        traits.loc[traits[oxygen_cols].sum(axis=1) != 1, 'Oxygen_Preference'] = np.nan
        
        # Gram Status
        # Initialize the column with 0
        traits['Gram_Status'] = 0
        # Set values based on conditions
        traits.loc[traits['Gram negative'] == 1, 'Gram_Status'] = "negative"
        traits.loc[traits['Gram positive'] == 1, 'Gram_Status'] = "positive"
        # Set to NA if the row sum is not 1
        gram_cols = ['Gram negative', 'Gram positive']
        traits.loc[traits[gram_cols].sum(axis=1) != 1, 'Gram_Status'] = np.nan
        
        # Cell Shape
        # Initialize the column with 0
        traits['cell_shape'] = 0
        # Set values based on conditions
        traits.loc[traits['Coccus'] == 1, 'cell_shape'] = 1
        traits.loc[traits['Bacillus or coccobacillus'] == 1, 'cell_shape'] = 2
        # Set to NA if the row sum is not 1
        shape_cols = ['Coccus', 'Bacillus or coccobacillus']
        traits.loc[traits[shape_cols].sum(axis=1) != 1, 'cell_shape'] = np.nan
    
        traits_traitor = traits[["Oxygen_Preference", "Gram_Status", "Motile", "Spore formation"]]
        traits_traitor.columns = ["Oxygen_Preference", "Gram_Status", "Motility", "Spore_Formation"]
    
        # --- 1. Load and prepare the traits data ---
        # Read the CSV, drop duplicates based on 'X16s_ID', and set it as the index
        traits = pd.read_csv("Data/bacDive.csv")
        traits.drop_duplicates(subset='16s_ID', inplace=True)
        traits.set_index('16s_ID', inplace=True)
        
        # Extract the accession number by splitting the index string at the period '.'
        accessions_num = np.array([str(x.split('.')[0]) for x in fid])
        
        # Create a helper DataFrame to map the full embed_id to the shortened accession number
        df_map = pd.DataFrame({
            'accessions': accessions_num,
            'embed_id': fid
        })
        
        # --- 3. Align traits and embedding data ---
        # Find the common accession numbers between the two datasets
        inter_id = np.intersect1d(traits.index, df_map['accessions'])
        
        # Filter the mapping and traits DataFrames to keep only the common entries
        df_map = df_map[df_map['accessions'].isin(inter_id)]
        traits = traits.loc[df_map['accessions']]
        
        # Update the traits index to match the full embedding ID for consistency
        traits.index = df_map['embed_id']
        
        # --- 4. Clean and standardize data in the traits DataFrame ---
        # Create a dictionary for all values that need to be replaced
        replace_dict = {"NA": np.nan, "-": "no", "+": "yes",
            "": np.nan, "+;NA": np.nan, "coccus-shaped": "coccus",
            "rod-shaped": "rod", "mixed": np.nan, "negative;variable": np.nan,
            "no;yes": np.nan, "negative;positive": np.nan, "variable": np.nan
                        
        }
        traits.replace(replace_dict, inplace=True)
        
        # Standardize the 'cell_shape' column: keep only 'coccus' or 'rod', set others to NaN
        valid_shapes = ['coccus', 'rod']
        traits['cell_shape'] = traits['cell_shape'].where(traits['cell_shape'].isin(valid_shapes), np.nan)
        
        # --- 5. Create the 'Oxygen.Preference' column ---
        # Define conditions and corresponding choices for oxygen preference
        conditions = [
            traits['aerobe'] == 1,
            traits['facultative.anaerobe'] == 1,
            traits['anaerobe'] == 1
        ]
        choices = ['aerobic', 'facultatively', 'anaerobic']
        
        # Use np.select (similar to R's case_when) to create the new column
        # The default value is NaN for anything that doesn't meet a condition
        traits['Oxygen.Preference'] = np.select(conditions, choices, default=np.nan)
        
        # --- 6. Final cleanup ---
        # Define columns to remove
        remove_cols = ['X16s_ID', 'aerobe', 'facultative.anaerobe', 'anaerobe']
        # Drop the specified columns; 'errors='ignore'' prevents an error if a column is already gone
        traits.drop(columns=remove_cols, inplace=True, errors='ignore')
        
        # Drop any column where all values are missing (NaN)
        traits.dropna(axis=1, how='all', inplace=True)
    
        traits_bacdive = traits[["Oxygen.Preference", "gram_stain", "motility", "spore_formation"]]
        traits_bacdive.columns = ["Oxygen_Preference", "Gram_Status", "Motility", "Spore_Formation"]
        traits_bacdive[traits_bacdive == "yes"] = 1
        traits_bacdive[traits_bacdive == "no"] = 0

        j = "Oxygen_Preference"
        phylum_id = ["Bacillota", "Bacteroidota", "Actinomycetota", "Pseudomonadota"]
        for i in phylum_id:
            temp = traits_bacdive.loc[traits_bacdive.loc[:, j].values != 'nan']
            temp = temp.loc[temp.loc[:, j].values == temp.loc[:, j].values]
            tax_bacdive = taxonomy.loc[temp.index.values]
            phylum = tax_bacdive.p.unique()
            test_phylum = list(phylum[[i not in phylum_id for i in phylum]]) + [i]
            test_id = tax_bacdive.loc[[i in test_phylum for i in tax_bacdive.p]].index.values
            test_id = test_id[traits_bacdive.loc[test_id, j].values == traits_bacdive.loc[test_id, j].values]
            
            tax_traitor = taxonomy.loc[traits_traitor.index.values]
            phylum = tax_traitor.p.unique()
            train_id = tax_traitor.loc[[i not in test_phylum for i in tax_traitor.p]].index.values
        
            train_id = train_id[traits_traitor.loc[train_id, j].values == traits_traitor.loc[train_id, j].values]
            X_train = co_embedding.loc[train_id]
            y_train = traits_traitor.loc[train_id, j].values
            X_test = co_embedding.loc[test_id]
            y_test = traits_bacdive.loc[test_id, j].values
        
            rf_model.fit(X_train, y_train)
            y_pred = rf_model.predict(X_test)
            y_pred_proba = rf_model.predict_proba(X_test)
            
            countsunique_elements, counts = np.unique(y_test, return_counts=True)
            if np.all(counts >= 4) and len(countsunique_elements) > 2:
                auc.append(roc_auc_score(y_test, y_pred_proba, multi_class='ovr', average='macro'))
                traits_type.append(j)
                group.append(f"times_{n}")
                tax.append(i)
                datasize.append(datasize_dict[t])


        j = "Gram_Status"
        phylum_id = ["Bacillota", "Bacteroidota", "Actinomycetota"]
        for i in phylum_id:
            temp = traits_bacdive.loc[traits_bacdive.loc[:, j].values != 'nan']
            temp = temp.loc[temp.loc[:, j].values == temp.loc[:, j].values]
            tax_bacdive = taxonomy.loc[temp.index.values]
            phylum = tax_bacdive.p.unique()
            test_phylum = list(phylum[[i not in phylum_id for i in phylum]]) + [i]
            test_id = tax_bacdive.loc[[i in test_phylum for i in tax_bacdive.p]].index.values
            test_id = test_id[traits_bacdive.loc[test_id, j].values == traits_bacdive.loc[test_id, j].values]
            
            tax_traitor = taxonomy.loc[traits_traitor.index.values]
            phylum = tax_traitor.p.unique()
            train_id = tax_traitor.loc[[i not in test_phylum for i in tax_traitor.p]].index.values
        
            train_id = train_id[traits_traitor.loc[train_id, j].values == traits_traitor.loc[train_id, j].values]
            X_train = co_embedding.loc[train_id]
            y_train = traits_traitor.loc[train_id, j].values
            X_test = co_embedding.loc[test_id]
            y_test = traits_bacdive.loc[test_id, j].values
            
            countsunique_elements, counts = np.unique(y_test, return_counts=True)
            if np.all(counts > 5) and len(countsunique_elements) > 1:
                rf_model.fit(X_train, y_train)
                y_pred = rf_model.predict(X_test)
                y_pred_proba = rf_model.predict_proba(X_test)
                
                auc.append(roc_auc_score(y_test, y_pred_proba[:,1]))
                traits_type.append(j)
                group.append(f"times_{n}")
                tax.append(i)
                datasize.append(datasize_dict[t])

                    
        j = "Motility"
        phylum_id = ["Bacillota", "Bacteroidota", "Actinomycetota"]
        for i in phylum_id:
            temp = traits_bacdive.loc[traits_bacdive.loc[:, j].values != 'nan']
            temp = temp.loc[temp.loc[:, j].values == temp.loc[:, j].values]
            tax_bacdive = taxonomy.loc[temp.index.values]
            phylum = tax_bacdive.p.unique()
            test_phylum = list(phylum[[i not in phylum_id for i in phylum]]) + [i]
            test_id = tax_bacdive.loc[[i in test_phylum for i in tax_bacdive.p]].index.values
            test_id = test_id[traits_bacdive.loc[test_id, j].values == traits_bacdive.loc[test_id, j].values]
            
            tax_traitor = taxonomy.loc[traits_traitor.index.values]
            phylum = tax_traitor.p.unique()
            train_id = tax_traitor.loc[[i not in test_phylum for i in tax_traitor.p]].index.values
        
            train_id = train_id[traits_traitor.loc[train_id, j].values == traits_traitor.loc[train_id, j].values]
            X_train = co_embedding.loc[train_id]
            y_train = traits_traitor.loc[train_id, j].values
            X_test = co_embedding.loc[test_id]
            y_test = traits_bacdive.loc[test_id, j].values.astype(int)
            
            countsunique_elements, counts = np.unique(y_test, return_counts=True)
            if np.all(counts > 5) and len(countsunique_elements) > 1:
                rf_model.fit(X_train, y_train)
                y_pred = rf_model.predict(X_test)
                y_pred_proba = rf_model.predict_proba(X_test)
                
                auc.append(roc_auc_score(y_test, y_pred_proba[:, 1]))
                traits_type.append(j)
                group.append(f"times_{n}")
                tax.append(i)
                datasize.append(datasize_dict[t])

                    
        j = "Spore_Formation"
        phylum_id = ["Bacillota", "Bacteroidota", "Actinomycetota", "Pseudomonadota"]
        for i in phylum_id:
            temp = traits_bacdive.loc[traits_bacdive.loc[:, j].values != 'nan']
            temp = temp.loc[temp.loc[:, j].values == temp.loc[:, j].values]
            tax_bacdive = taxonomy.loc[temp.index.values]
            phylum = tax_bacdive.p.unique()
            test_phylum = list(phylum[[i not in phylum_id for i in phylum]]) + [i]
            test_id = tax_bacdive.loc[[i in test_phylum for i in tax_bacdive.p]].index.values
            test_id = test_id[traits_bacdive.loc[test_id, j].values == traits_bacdive.loc[test_id, j].values]
            
            tax_traitor = taxonomy.loc[traits_traitor.index.values]
            phylum = tax_traitor.p.unique()
            train_id = tax_traitor.loc[[i not in test_phylum for i in tax_traitor.p]].index.values
        
            train_id = train_id[traits_traitor.loc[train_id, j].values == traits_traitor.loc[train_id, j].values]
            X_train = co_embedding.loc[train_id]
            y_train = traits_traitor.loc[train_id, j].values
            X_test = co_embedding.loc[test_id]
            y_test = traits_bacdive.loc[test_id, j].values.astype(int)
            
            countsunique_elements, counts = np.unique(y_test, return_counts=True)
            if np.all(counts > 5) and len(countsunique_elements) > 1:
                rf_model.fit(X_train, y_train)
                y_pred = rf_model.predict(X_test)
                y_pred_proba = rf_model.predict_proba(X_test)

                auc.append(roc_auc_score(y_test, y_pred_proba[:,1]))
                traits_type.append(j)
                group.append(f"times_{n}")
                tax.append(i)
                datasize.append(datasize_dict[t])

res = pd.DataFrame({"auc":auc, "traits_type":traits_type, "group":group, "tax":tax, "datasize":datasize})
res.to_csv("Data/auc_res.csv", index=None)


