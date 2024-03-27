import re
import pandas as pd
import requests
from io import StringIO
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqUtils import ProtParam
from DataBase import MongoDB, ProteinEntry, SeqsEntry, DescEntry, TaxonEntry


def get_uniprot_data(query):
    tsv_base_url = 'https://rest.uniprot.org/uniprotkb/stream'
    params = {
        'fields': 'accession,reviewed,id,protein_name,gene_names,'
                  'organism_name,length,ec,kinetics,ph_dependence,'
                  'mass,gene_orf,gene_primary,gene_oln,gene_synonym,'
                  'keyword,keywordid,annotation_score,lit_pubmed_id,'
                  'lit_doi_id,protein_families,xref_alphafolddb,xref_pdb,'
                  'xref_geneid,xref_brenda,xref_pfam,go_f,go_p,ft_signal,'
                  'xref_cazy,xref_signalink,temp_dependence,organism_id',
        'format': 'tsv',
        'query': f'({query})'
    }

    try:
        response = requests.get(tsv_base_url, params=params)
        response.raise_for_status()

        data = pd.read_csv(StringIO(response.text), delimiter='\t')
        data.to_csv('inulinases.csv', index=False)

        return data

    except requests.exceptions.RequestException as error:
        print(f'Falha na solicitação: {error}')
        return None


def get_uniprot_fasta(query):
    fasta_base_url = 'https://rest.uniprot.org/uniprotkb/stream'
    params = {
        'format': 'fasta',
        'query': f'({query})'
    }

    try:
        response = requests.get(fasta_base_url, params=params)
        response.raise_for_status()

        with open('inulinases_uniprot.fasta', 'w') as file:
            file.write(response.text)

        return 'inulinases_uniprot.fasta'

    except requests.exceptions.RequestException as error:
        print(f'Falha na solicitação: {error}')
        return None


def filter_inulinases():
    df = get_uniprot_data('inulinase')
    inulinase_selection = df[
        df['Protein names'].str.contains('inulinase|Inulinase', case=False)]
    inulinase_selection = inulinase_selection.reset_index(drop=True)
    inulinase_selection.to_csv('inulinases_clean.csv')
    return inulinase_selection


def only_inulinases():
    fasta = get_uniprot_fasta('inulinase')
    df = filter_inulinases()
    seqs_id = df['Entry']
    seqs_id = set(seqs_id)
    seqs_found = []

    with open(fasta, 'r') as arq:
        for record in SeqIO.parse(arq, 'fasta'):
            entry = record.id.split('|')[1]
            if entry in seqs_id:
                record.id = f'{entry}'
                record.description = ''
                seqs_found.append(record)

    with open('only_inulinases.fasta', 'w') as out:
        SeqIO.write(seqs_found, out, 'fasta')

    return seqs_found


def csv_to_dict():
    df = filter_inulinases()
    entries_list = []

    def km(string):
        if 'KM=' in string:
            regex = r'KM=(\d+\.\d+)\s+mM'
            match = re.search(regex, string)
            if match:
                return match.group(1)
        return None

    def pH(string):
        if 'pH' in string:
            regex = r'Optimum pH is (\d+\.\d+)'
            match = re.search(regex, string)
            if match:
                return match.group(1)
        return None

    def temp(string):
        if 'temperature is' in string:
            regex = r'temperature is (\d+)'
            match = re.search(regex, string)
            if match:
                return match.group(1)
        return None

    def keywords_split(string):
        if pd.isna(string):
            return None
        return string.split(';')

    def pfam_split(string):
        if pd.isna(string):
            return None
        pfam_list = string.split(';')
        return [pfam.strip() for pfam in pfam_list if pfam.strip()]

    def pdb_split(string):
        if pd.isna(string) or string.lower() == 'nan':
            return None
        list_pdb = string.split(';')
        return [pdb.strip() for pdb in list_pdb if pdb.strip()]

    def signalp_split(string):
        if pd.isna(string) or str(string).lower() == 'nan':
            return None
        else:
            regex = r'SIGNAL (\d+\..\d+)'
            match = re.search(regex, str(string))
            if match:
                return match.group(1)

    def pubmed_split(string):
        if pd.isna(string) or string.lower() == 'nan':
            return None
        list_pubmed = string.split(';')
        return [pubmed.strip() for pubmed in list_pubmed if pubmed.strip()]

    def doi_split(string):
        if pd.isna(string) or string.lower() == 'nan':
            return None
        list_doi = string.split(';')
        return [doi.strip() for doi in list_doi if doi.strip()]

    for index, row in df.iterrows():
        entry_dict = {
            "entry": row['Entry'],
            "entry_name": row['Entry Name'],
            "protein_names": row['Protein names'],
            "gene_names": row['Gene Names'],
            "organism": row['Organism'],
            "length": row['Length'],
            "ec_number": row['EC number'],
            "kinetics": km(str(row['Kinetics'])),
            "ph_dependence": pH(str(row['pH dependence'])),
            "temperature_dependence": temp(
                str(row['Temperature dependence'])),
            "mass": row['Mass'],
            "keywords": keywords_split(row['Keywords']),
            "alphafold_db": row['AlphaFoldDB'],
            "pdb": pdb_split(row['PDB']),
            "pfam": pfam_split(row['Pfam']),
            "signal_peptide": signalp_split(row['Signal peptide']),
            "pubmed_id": pubmed_split(row['PubMed ID']),
            "doi_id": doi_split(row['DOI ID'])
        }

        entries_list.append(entry_dict)

    return entries_list


def fasta_to_dict():
    list_seqs = only_inulinases()
    seqs_list = []

    for record in list_seqs:
        entry = f'>{record.id}'
        fasta_seq = str(record.seq)
        seqs_data = {"header": entry, "seq": fasta_seq}
        seqs_list.append(seqs_data)
    return seqs_list


def calc_descriptor(header, seq):
    def verify_x(s):
        for i in s:
            if i == 'X':
                return "X found"
            else:
                continue

    verify = verify_x(seq)

    if verify != "X found":

        protein_param = ProtParam.ProteinAnalysis(seq)

        desc_dict = {
            'header': header[1:],
            'aa_percent': protein_param.get_amino_acids_percent(),
            'charge_ph': protein_param.charge_at_pH(7),
            'isoelectric_point': protein_param.isoelectric_point(),
            'aromaticity': protein_param.aromaticity(),
            'flexibility': protein_param.flexibility(),
            'gravy': protein_param.gravy(),
            'sec_struc_frac': protein_param.secondary_structure_fraction(),
        }

        return desc_dict

    else:
        desc_dict = {
            'header': header[1:],
            'aa_percent': None,
            'charge_ph': None,
            'isoelectric_point': None,
            'aromaticity': None,
            'flexibility': None,
            'gravy': None,
            'sec_struc_frac': None,
        }

        return desc_dict


def get_taxon():

    Entrez.email = "lspalmeira.bio@gmail.com"

    df = filter_inulinases()

    list_dict = []

    for index, row in df.iterrows():
        organism_id = row['Organism (ID)']

        rank_mapping = {
            "superkingdom": "Superkingdom",
            "kingdom": "Kingdom",
            "phylum": "Phylum",
            "class": "Class",
            "order": "Order",
            "family": "Family",
            "genus": "Genus"
        }

        try:
            handle = Entrez.efetch(db="taxonomy", id=organism_id,
                                   retmode="xml")
            taxon_record = Entrez.read(handle)

            taxon_details = {'entry': row['Entry']}
            for item in taxon_record[0]["LineageEx"]:
                rank = item["Rank"]
                if rank in rank_mapping:
                    taxon_details[rank_mapping[rank]] = item["ScientificName"]

            list_dict.append(taxon_details)

        except Exception as e:
            taxon_none = {
                'entry': row['Entry'],
                'Superkingdom': None,
                'Kingdom': None,
                'Phylum': None,
                'Class': None,
                'Order': None,
                'Family': None,
                'Genus': None
            }

            list_dict.append(taxon_none)

    return list_dict


if __name__ == "__main__":

    list_dict = csv_to_dict()
    seqs_fasta = fasta_to_dict()

    db = MongoDB()
    db.connect_to_mongodb()

    for data_dict in list_dict:
        protein = ProteinEntry(**data_dict)
        db.insert_data(protein)

    for seqs_data in seqs_fasta:
        seqs_entry = SeqsEntry(**seqs_data)
        db.insert_seqs_data(seqs_entry)

    for seqs_data in seqs_fasta:
        desc_calc = calc_descriptor(seqs_data['header'], seqs_data['seq'])
        desc_entry = DescEntry(**desc_calc)
        db.insert_desc_data(desc_entry)

    list_taxon = get_taxon()

    for dict_taxon in list_taxon:
        taxon = TaxonEntry(**dict_taxon)
        db.insert_taxon_data(taxon)
