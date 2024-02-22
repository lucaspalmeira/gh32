import pymongo


class ProteinEntry:
    def __init__(self, entry=None, entry_name=None,
                 protein_names=None, gene_names=None, organism=None,
                 length=None, ec_number=None, kinetics=None,
                 ph_dependence=None, temperature_dependence=None, mass=None,
                 keywords=None, alphafold_db=None, pdb=None, pfam=None,
                 signal_peptide=None, pubmed_id=None, doi_id=None):

        # Cada variável recebe um valor
        self.entry = entry
        self.entry_name = entry_name
        self.protein_names = protein_names
        self.gene_names = gene_names
        self.organism = organism
        self.length = length
        self.ec_number = ec_number
        self.kinetics = kinetics
        self.ph_dependence = ph_dependence
        self.temperature_dependence = temperature_dependence
        self.mass = mass
        self.keywords = keywords
        self.alphafold_db = alphafold_db
        self.pdb = pdb
        self.pfam = pfam
        self.signal_peptide = signal_peptide
        self.pubmed_id = pubmed_id
        self.doi_id = doi_id

    def to_dict(self):
        # Um dicionário é retornado
        return {
            "entry": self.entry,
            "entry_name": self.entry_name,
            "protein_names": self.protein_names,
            "gene_names": self.gene_names,
            "organism": self.organism,
            "length": self.length,
            "ec_number": self.ec_number,
            "kinetics": self.kinetics,
            "ph_dependence": self.ph_dependence,
            "temperature_dependence": self.temperature_dependence,
            "mass": self.mass,
            "keywords": self.keywords,
            "alphafold_db": self.alphafold_db,
            "pdb": self.pdb,
            "pfam": self.pfam,
            "signal_p": self.signal_peptide,
            "pubmed_id": self.pubmed_id,
            "doi_id": self.doi_id
        }

class MongoDB:
    def __init__(self):
        self.client = None
        self.db = None
        self.collection = None

    def connect_to_mongodb(self):
        try:
            self.client = pymongo.MongoClient("mongodb://localhost:27017/")
            self.db = self.client["inulinases_database"]
            # getting database
            self.collection = self.db["protein_entries"]
            # getting collection
            print("Connected to MongoDB")
        except pymongo.errors.ConnectionFailure:
            print("Failed to connect to MongoDB")

    def insert_data(self, protein):
        entry_data = protein.to_dict()
        insertion_result = self.collection.insert_one(entry_data)
        return insertion_result.inserted_id

    def retrieve_protein(self, entry_id):
        result = self.collection.find_one({"entry": entry_id})
        if result:
            valid_arguments = {key: result[key] for key in ProteinEntry.__init__.__code__.co_varnames if key in result}
            return ProteinEntry(**valid_arguments)
        return None

    def close_connection(self):
        if self.client:
            self.client.close()
            print("Connection to MongoDB closed")


if __name__ == "__main__":
    db = MongoDB()
    db.connect_to_mongodb()

    protein_data = {
        "entry": "P12345",
        "entry_name": "PROT123",
        "protein_names": ["Protein A", "Protein B"],
        "gene_names": ["Gene A", "Gene B"],
        "organism": "Example Organism",
        "length": 300,
        "ec_number": "1.2.3.4",
        "kinetics": 5.5,
        "ph_dependence": 7.5,
        "temperature_dependence": 37.0,
        "mass": 50.5,
        "keywords": ["inu", "gh32"],
        "alphafold_db": "ASHI086",
        "pdb": ["1XTC", "1PZM"],
        "pfam": ["PFAM1", "PFAM2"],
        "signal_peptide": "1...25",
        "pubmed_id": "PMID12345",
        "doi_id": "DOI123"
    }

    protein = ProteinEntry(**protein_data)
    db.insert_data(protein)

    retrieved_protein = db.retrieve_protein("P12345")
    if retrieved_protein:
        print("Retrieved Protein Entry:")
        print(retrieved_protein.to_dict())

    db.close_connection()
