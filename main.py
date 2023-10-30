from DataBase import SQLdb
import pandas as pd


execute = SQLdb

execute.create_table('inulinases.db')
def read_csv(file):

    list_string = []


    df = pd.read_csv(file)
    df = df[['Entry', 'Entry Name', 'Protein names', 'Gene Names', 'Organism',
    'Length', 'EC number', 'PubMed ID', 'PDB', 'AlphaFoldDB', 'Pfam']]

    for i, row in df.iterrows():
        x = (f'{row["Entry"]}, {row["Entry Name"]}, '
             f'{row["Protein names"]}, {row["Gene Names"]}, '
             f'{row["Organism"]}, {row["Length"]}, {row["EC number"]}, '
             f'{row["PubMed ID"]}, {row["PDB"]}, {row["AlphaFoldDB"]}, '
             f'{row["Pfam"]}')
        print(x)
        list_string.append(x)

    for x in list_string:
        print(x)

    return list_string

read_csv('inulinases.csv')

execute.insert_dataset()
