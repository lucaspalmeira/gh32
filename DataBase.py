import sqlite3
import pymongo

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

    def consult_select(self):
        result = self.cursor.execute("SELECT * FROM inulinases")
        for row in result:
            print(row)

    def select_all_by_column_index(self, column_index):
        query = f"SELECT * FROM inulinases"
        self.cursor.execute(query)
        values = [row[column_index] for row in self.cursor.fetchall()]

        return values

    def close_connection(self):
        self.connection.close()


class MONGOdb:
    def __int__(self):
        def __init__(self, entry, entry_name, protein_names, gene_names,
                     organism, length, ec_number, pubmed_id, pdb,
                     alphafold_db, pfam):
            self.entry = entry
            self.entry_name = entry_name
            self.protein_names = protein_names
            self.gene_names = gene_names
            self.organism = organism
            self.length = length
            self.ec_number = ec_number
            self.pubmed_id = pubmed_id
            self.pdb = pdb
            self.alphafold_db = alphafold_db
            self.pfam = pfam

        def to_dict(self):
            return {
                "entry": self.entry,
                "entry_name": self.entry_name,
                "protein_names": self.protein_names,
                "gene_names": self.gene_names,
                "organism": self.organism,
                "length": self.length,
                "ec_number": self.ec_number,
                "pubmed_id": self.pubmed_id,
                "pdb": self.pdb,
                "alphafold_db": self.alphafold_db,
                "pfam": self.pfam
            }

        def __str__(self):
            return f"ProteinEntry(entry={self.entry}, entry_name={self.entry_name}, organism={self.organism})"

# Conecte-se ao servidor MongoDB (por padrão, ele se conectará ao servidor local)
client = pymongo.MongoClient("mongodb://localhost:27017/")

# Crie um banco de dados (se não existir)
db = client["proteins_database"]

# Crie uma coleção para armazenar as informações de proteínas
collection = db["protein_entries"]

# Crie uma instância da classe ProteinEntry com os dados
protein = ProteinEntry(
    entry="P12345",
    entry_name="ProteinA",
    protein_names=["Protein A"],
    gene_names=["Gene A"],
    organism="Example Organism",
    length=300,
    ec_number="2.1.1.1",
    pubmed_id="12345",
    pdb="PDB123",
    alphafold_db="AlphaFoldDB",
    pfam="PF12345"
)

# Insira os dados no banco de dados
entry_data = protein.to_dict()
insertion_result = collection.insert_one(entry_data)

print("ID do documento inserido:", insertion_result.inserted_id)

# Consultar o banco de dados
result = collection.find_one({"entry": "P12345"})
if result:
    retrieved_protein = ProteinEntry(**result)
    print("Proteína recuperada:", retrieved_protein)

# Feche a conexão com o servidor MongoDB
client.close()


if __name__ == "__main__":
    db = SQLdb()
    db.create_table()
    # Você precisa fornecer dados para insert_data e insert_dataset.
    # Exemplo:
    # db.insert_data("value1", "value2", "value3", "value4", "value5", "value6", "value7", "value8", "value9", "value10", "value11")
    # db.insert_dataset(["value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,value11"])
    db.close_connection()