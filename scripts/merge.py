######################################################################################
#                                                                                    #
# This script merges two MongoDB collections into a single CSV file.                 #
# It connects to a MongoDB instance, retrieves data from the specified collections,  #
# and merges them based on a common field ('entry').                                 #
#                                                                                    #
######################################################################################


from pymongo import MongoClient
import pandas as pd

client = MongoClient('mongodb://172.17.0.2:27017/')
db = client['gh32']

protein_entries = db['protein_entries']
taxon_entries = db['taxon_entries']

# Converter as collections em DataFrames do Pandas
df_protein = pd.DataFrame(list(protein_entries.find()))
df_taxon = pd.DataFrame(list(taxon_entries.find()))

# Remover o campo '_id' (opcional)
df_protein = df_protein.drop(columns=['_id'], errors='ignore')
df_taxon = df_taxon.drop(columns=['_id'], errors='ignore')

# Unir as duas collections com base no campo 'entry'
df_merged = pd.merge(df_protein, df_taxon, on='entry', how='inner')


df_merged.to_csv('merge_protein_entries_taxon_entries.csv', index=False)

print('CSV criado.')