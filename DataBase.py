import sqlite3

class SQLdb:
    def __init__(self):
        self.connection = sqlite3.connect('inulinases.db')
        self.cursor = self.connection.cursor()

    def create_table(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS inulinases (
                id INTEGER PRIMARY KEY,
                Entry TEXT NOT NULL,
                Entry_Name TEXT NOT NULL,
                Protein_Names TEXT NOT NULL,
                Gene_Names TEXT,  -- Specify the data type for 'Gene Names'
                Organism TEXT NOT NULL,
                Length TEXT NOT NULL,
                EC_Number TEXT,  -- Allow null values for 'EC_Number'
                PubMed_ID TEXT,  -- Allow null values for 'PubMed_ID'
                PDB TEXT,
                AlphaFoldDB TEXT,
                Pfam TEXT
            )
        ''')
        self.connection.commit()

    def insert_data(self, entry, entry_name, Protein_Names,
                    Gene_Names, Organism, Length, EC_Number,
                    PubMed_ID, PDB, AlphaFoldDB, Pfam):
        self.cursor.execute("INSERT INTO inulinases "
                            "(Entry, Entry_Name, Protein_Names, Gene_Names,"
                            "Organism, Length, EC_Number, PubMed_ID, PDB,"
                            "AlphaFoldDB, Pfam) VALUES ("
                            "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?"
                            ")",
                            (entry, entry_name, Protein_Names,
                             Gene_Names, Organism, Length, EC_Number,
                             PubMed_ID, PDB, AlphaFoldDB, Pfam))
        self.connection.commit()

    def insert_dataset(self, list_strings):
        for _str in list_strings:
            _data = _str.split(',')
            if len(_data) == 11:
                (entry, entry_name, Protein_Names, Gene_Names, Organism,
                 Length, EC_Number, PubMed_ID, PDB, AlphaFoldDB, Pfam) = _data
                self.insert_data(entry.strip(), entry_name.strip(),
                                 Protein_Names.strip(), Gene_Names.strip(),
                                 Organism.strip(), Length.strip(),
                                 EC_Number.strip(), PubMed_ID.strip(),
                                 PDB.strip(), AlphaFoldDB.strip(), Pfam.strip())

    def remove_entry(self, entry_id):
        self.cursor.execute("DELETE FROM inulinases WHERE id=?", (entry_id,))
        self.connection.commit()

    def close_connection(self):
        self.connection.close()

"""
if __name__ == "__main__":
    db = SQLdb()
    db.create_table()
    # You need to provide data to insert_data and insert_dataset.
    # Example:
    # db.insert_data("value1", "value2", "value3", "value4", "value5", "value6", "value7", "value8", "value9", "value10", "value11")
    # db.insert_dataset(["value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,value11"])
    db.close_connection()
"""