import pymongo


class ProteinEntry:
    def __init__(self, entry=None, entry_name=None,
                 protein_names=None, gene_names=None, organism=None,
                 length=None, ec_number=None, kinetics=None,
                 ph_dependence=None, temperature_dependence=None, mass=None,
                 keywords=None, alphafold_db=None, pdb=None, pfam=None,
                 signal_peptide=None, pubmed_id=None, doi_id=None):

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
        return {
            'entry': self.entry,
            'entry_name': self.entry_name,
            'protein_names': self.protein_names,
            'gene_names': self.gene_names,
            'organism': self.organism,
            'length': self.length,
            'ec_number': self.ec_number,
            'kinetics': self.kinetics,
            'ph_dependence': self.ph_dependence,
            'temperature_dependence': self.temperature_dependence,
            'mass': self.mass,
            'keywords': self.keywords,
            'alphafold_db': self.alphafold_db,
            'pdb': self.pdb,
            'pfam': self.pfam,
            'signal_p': self.signal_peptide,
            'pubmed_id': self.pubmed_id,
            'doi_id': self.doi_id
        }


class SeqsEntry:
    def __init__(self, entry, seq):
        self.entry = entry
        self.seq = seq

    def to_dict(self):
        return {
        'entry': self.entry,
        'seq': self.seq
        }


class DescEntry:
    def __init__(self, entry, aa_percent=None, charge_ph=None,
                 isoelectric_point=None, aromaticity=None,
                 flexibility=None, gravy=None, sec_struc_frac=None):

        self.entry = entry
        self.aa_percent = aa_percent
        self.charge_ph = charge_ph
        self.isoelectric_point = isoelectric_point
        self.aromaticity = aromaticity
        self.flexibility = flexibility
        self.gravy = gravy
        self.sec_struc_frac = sec_struc_frac

    def to_dict(self):
        return {
            'entry': self.entry,
            'aa_percent': self.aa_percent,
            'charge_ph': self.charge_ph,
            'isoelectric_point': self.isoelectric_point,
            'aromaticity': self.aromaticity,
            'flexibility': self.flexibility,
            'grav': self.gravy,
            'sec_struc_frac': self.sec_struc_frac,
        }


class TaxonEntry:
    def __init__(self, entry, Superkingdom=None, Kingdom=None, Phylum=None,
                 Class=None, Order=None, Family=None, Genus=None):
        self.entry = entry
        self.Superkingdom = Superkingdom
        self.Kingdom = Kingdom
        self.Phylum = Phylum
        self.Class = Class
        self.Order = Order
        self.Family = Family
        self.Genus = Genus

    def to_dict(self):
        return {
            'entry': self.entry,
            'Superkingdom': self.Superkingdom,
            'Kingdom': self.Kingdom,
            'Phylum': self.Phylum,
            'Class': self.Class,
            'Order': self.Order,
            'Family': self.Family,
            'Genus': self.Genus
        }


class MongoDB:
    def __init__(self):
        self.client = None
        self.db = None
        self.protein_collection = None
        self.seqs_collection = None
        self.desc_collection = None
        self.taxon_collection = None
        self.collections_to_filter = None

    def connect_to_mongodb(self):
        try:
            self.client = pymongo.MongoClient('mongodb://localhost:27017/')
            self.db = self.client['gh32']
            self.protein_collection = self.db['protein_entries']
            self.seqs_collection = self.db['seqs_entries']
            self.desc_collection = self.db['desc_entries']
            self.taxon_collection = self.db['taxon_entries']
            
            self.collections_to_filter = ['protein_entries',
                                          'desc_entries',
                                          'seqs_entries',
                                          'taxon_entries']
            print('Connected to MongoDB')
        except pymongo.errors.ConnectionFailure:
            print('Failed to connect to MongoDB')

    def db_exists(self):
        list_db = self.client.list_database_names()
        return list_db

    def get_entries(self):
        entries = []
        collection = self.protein_collection
        values = collection.find({}, {'entry': 1, '_id': 0})
        for dic in values:
            entries.append(dic['entry'])
        return entries

    def insert_data(self, protein):
        entry_data = protein.to_dict()
        insertion_result = self.protein_collection.insert_one(entry_data)
        return insertion_result.inserted_id

    def retrieve_protein(self, entry_id):
        result = self.protein_collection.find_one({'entry': entry_id})
        if result:
            valid_arguments = {key: result[key] for key in
                               ProteinEntry.__init__.__code__.co_varnames if
                               key in result}
            return ProteinEntry(**valid_arguments)
        return None

    def insert_seqs_data(self, seqs_entry):
        entry_data = seqs_entry.to_dict()
        insertion_result = self.seqs_collection.insert_one(entry_data)
        return insertion_result.inserted_id

    def retrieve_seq(self, entry_id):
        result = self.seqs_collection.find_one({'entry': entry_id})
        if result:
            valid_arguments = {key: result[key] for key in
                               SeqsEntry.__init__.__code__.co_varnames if
                               key in result}
            return SeqsEntry(**valid_arguments)
        return None

    def insert_desc_data(self, desc_entry):
        desc_data = desc_entry.to_dict()
        insertion_result = self.desc_collection.insert_one(desc_data)
        return insertion_result.inserted_id

    def retrieve_desc(self, entry_id):
        result = self.desc_collection.find_one({'entry': entry_id})
        if result:
            valid_arguments = {key: result[key] for key in
                               DescEntry.__init__.__code__.co_varnames if
                               key in result}
            return DescEntry(**valid_arguments)
        return None

    def insert_taxon_data(self, taxon_entry):
        taxon_data = taxon_entry.to_dict()
        insertion_result = self.taxon_collection.insert_one(taxon_data)
        return insertion_result.inserted_id

    def retrieve_taxon(self, entry_id):
        result = self.taxon_collection.find_one({'entry': entry_id})
        if result:
            valid_arguments = {key: result[key] for key in
                               TaxonEntry.__init__.__code__.co_varnames if
                               key in result}
            return TaxonEntry(**valid_arguments)
        return None

    def update_ec_number(self, dataframe):

        for index, row in dataframe.iterrows():
            entry = row['entry']
            new_ec_number = row['ec_number']

            document = self.protein_collection.find_one({"entry": entry})

            if document:
                self.protein_collection.update_one(
                    {"_id": document["_id"]},
                    {"$set": {"ec_number": new_ec_number}}
                )

    def remove_enzymes_not_in_gh32(self, entries):
        for collection_name in self.collections_to_filter:
            collection = self.db[collection_name]
            result = collection.delete_many({'entry': {'$in': entries}})
            print(f'Removed {result.deleted_count} documents '
                  f'from {collection_name}.')

    def close_connection(self):
        if self.client:
            self.client.close()
            print('Connection to MongoDB closed')
