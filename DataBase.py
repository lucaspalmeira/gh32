import sqlite3


class SQLdb:
    def __init__(self, name_bank):
        self.connection = sqlite3.connect(name_bank)
        self.cursor = self.connection.cursor()

    def create_table(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS inulinases (
                id INTEGER PRIMARY KEY,
                Entry TEXT NOT NULL,
                Entry_Name TEXT NOT NULL,
                Protein_Names TEXT NOT NULL,
                Gene Names,
                Organism NOT NULL,
                Length NOT NULL,
                EC_Number,
                PubMed_ID,
                PDB,
                AlphaFoldDB,
                Pfam
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

    def insert_dataset(self, list_):
        for i in list_:
            self.insert_data(i['entry'], i['entry_name'], i['Protein_Names'],
                             i['Gene_Names'], i['Organism'], i['Length'],
                             i['EC_Number'], i['PubMed_ID'], i['PDB'],
                             i['AlphaFoldDB'], i['Pfam'])

    def remover_cliente(self, entry_id):
        self.cursor.execute("DELETE FROM clientes WHERE id=?", (entry_id,))
        self.connection.commit()

    def close_connection(self):
        self.connection.close()


# Exemplo de uso da classe
if __name__ == "__main__":
    db = SQLdb('inulinases.db')
    db.create_table()
    db.insert_data()
    db.insert_dataset()
    db.close_connection()
