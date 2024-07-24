import os
import pandas as pd
import requests
from io import StringIO
from Bio import SeqIO
from Bio import Entrez
import sys, errno, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep
from Bio.SeqUtils import ProtParam
from DataBase import MongoDB, ProteinEntry, SeqsEntry, DescEntry, TaxonEntry


def parse_items(items):
    if type(items) == list:
        return ",".join(items)
    return ""


def parse_member_databases(dbs):
    if type(dbs) == dict:
        return ";".join([f"{db}:{','.join(dbs[db])}" for db in dbs.keys()])
    return ""


def parse_go_terms(gos):
    if type(gos) == list:
        return ",".join([go["identifier"] for go in gos])
    return ""


def parse_locations(locations):
    if type(locations) == list:
        return ",".join(
            [",".join([f"{fragment['start']}..{fragment['end']}"
                       for fragment in location["fragments"]]) for location in
             locations])
    return ""


def parse_group_column(values, selector):
    return ",".join([parse_column(value, selector) for value in values])


def parse_column(value, selector):
    if value is None:
        return ""
    elif "member_databases" in selector:
        return parse_member_databases(value)
    elif "go_terms" in selector:
        return parse_go_terms(value)
    elif "children" in selector:
        return parse_items(value)
    elif "locations" in selector:
        return parse_locations(value)
    return str(value)


def gh32_interpro():
    # disable SSL verification to avoid config issues
    context = ssl._create_unverified_context()

    BASE_URL = ('https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/'
                'entry/InterPro/IPR001362/taxonomy/uniprot/4751/?page_size=200')

    next_url = BASE_URL
    last_page = False

    data = []

    while next_url:
        try:
            req = request.Request(next_url,
                                  headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)
            if res.status == 408:
                sleep(61)
                continue
            elif res.status == 204:
                break
            payload = json.loads(res.read().decode())
            next_url = payload["next"]
            if not next_url:
                last_page = True
        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
                raise e

        for item in payload["results"]:
            row = [
                parse_column(item["metadata"]["accession"],
                             'metadata.accession'),
                parse_column(item["metadata"]["source_database"],
                             'metadata.source_database'),
                parse_column(item["metadata"]["name"],
                             'metadata.name'),
                parse_column(item["metadata"]["source_organism"]["taxId"],
                             'metadata.source_organism.taxId'),
                parse_column(
                    item["metadata"]["source_organism"]["scientificName"],
                    'metadata.source_organism.scientificName'),
                parse_column(item["metadata"]["length"],
                             'metadata.length'),
                parse_column(item["entries"][0]["accession"],
                             'entries[0].accession'),
                parse_column(item["entries"][0]["entry_protein_locations"],
                             'entries[0].entry_protein_locations')
            ]
            data.append(row)

        if next_url:
            sleep(1)

    df = pd.DataFrame(data, columns=['Accession', 'Source Database', 'Name',
                                     'Taxonomy ID', 'Scientific Name',
                                     'Length', 'Entry Accession',
                                     'Entry Protein Locations'])

    accession = df['Accession']

    df.to_csv('data_from_interpro.csv', index=False)

    return accession


def update_db(ids_interpro, ids_in_db):
    new_ids = set(ids_interpro) - set(ids_in_db)
    new_ids = list(new_ids)
    return new_ids


def submit_id_mapping_job(ids, from_db, to_db):
    url = 'https://rest.uniprot.org/idmapping/run'
    data = {
        'ids': ','.join(ids),
        'from': from_db,
        'to': to_db
    }
    response = requests.post(url, data=data)
    if response.status_code == 200:
        response_data = response.json()
        job_id = response_data.get('jobId')
        if job_id:
            return job_id
        else:
            print("Erro: jobId não encontrado na resposta.")
            return None
    else:
        print(f"Erro ao enviar solicitação: {response.status_code}")
        return None


def extract_id_mapping_results(job_id):
    tsv_base_url = 'https://rest.uniprot.org/idmapping/uniprotkb/results/stream'
    params = {
        'fields': 'accession,reviewed,id,protein_name,gene_names,'
                  'organism_name,length,ec,kinetics,ph_dependence,'
                  'mass,gene_orf,gene_primary,gene_oln,gene_synonym,'
                  'keyword,keywordid,annotation_score,lit_pubmed_id,'
                  'lit_doi_id,protein_families,xref_alphafolddb,xref_pdb,'
                  'xref_geneid,xref_brenda,xref_pfam,go_f,go_p,ft_signal,'
                  'xref_cazy,xref_signalink,temp_dependence,organism_id',
        'format': 'tsv'
    }

    url = f"{tsv_base_url}/{job_id}"
    response = requests.get(url, params=params)
    if response.status_code == 200:
        data = response.text

        df = pd.read_csv(StringIO(data), delimiter='\t')
        df.to_csv('gh32.csv', index=False)

        print("Dados extraídos com sucesso.")

        return df
    else:
        print(f"Falha ao extrair dados. Status code: {response.status_code}")


def get_uniprot_fasta(job_id):
    fasta_base_url = (f'https://rest.uniprot.org/idmapping/'
                      f'uniprotkb/results/stream/{job_id}')

    params = {
        'format': 'fasta'
    }

    try:
        response = requests.get(fasta_base_url, params=params)
        response.raise_for_status()

        with open('gh32.fasta', 'w') as file:
            file.write(response.text)

        print('Requisição de sequências realizada.')

        path_fasta_in_CLEAN_app = os.path.join(os.getcwd(),
                                               'CLEAN', 'app', 'data',
                                               'inputs', 'gh32.fasta')

        with open(path_fasta_in_CLEAN_app, 'w') as file2:
            file2.write(response.text)

        print(f'gh32.fasta criado com sucesso em: {path_fasta_in_CLEAN_app}.')

        return 'gh32.fasta'

    except requests.exceptions.RequestException as error:
        print(f'Download de sequências: falha na solicitação: {error}')
        return None


def csv_to_dict(df):
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
            'entry': row['Entry'],
            'entry_name': row['Entry Name'],
            'protein_names': row['Protein names'],
            'gene_names': row['Gene Names'],
            'organism': row['Organism'],
            'length': row['Length'],
            'ec_number': row['EC number'],
            'kinetics': km(str(row['Kinetics'])),
            'ph_dependence': pH(str(row['pH dependence'])),
            'temperature_dependence': temp(
                str(row['Temperature dependence'])),
            'mass': row['Mass'],
            'keywords': keywords_split(row['Keywords']),
            'alphafold_db': row['AlphaFoldDB'],
            'pdb': pdb_split(row['PDB']),
            'pfam': pfam_split(row['Pfam']),
            'signal_peptide': signalp_split(row['Signal peptide']),
            'pubmed_id': pubmed_split(row['PubMed ID']),
            'doi_id': doi_split(row['DOI ID'])
        }

        entries_list.append(entry_dict)

    return entries_list


def fasta_to_dict(list_seqs):
    seqs_list = []

    for record in list_seqs:
        entry = f'{record.id}'
        entry = entry.split('|')
        entry = entry[1]
        fasta_seq = str(record.seq)
        seqs_data = {'entry': entry, 'seq': fasta_seq}
        seqs_list.append(seqs_data)
    return seqs_list


def calc_descriptor(entry, seq):
    def verify_x(s):
        for i in s:
            if i == 'X':
                return 'X found'
            else:
                continue

    verify = verify_x(seq)

    if verify != 'X found':

        protein_param = ProtParam.ProteinAnalysis(seq)

        desc_dict = {
            'entry': entry,
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
            'entry': entry[1:],
            'aa_percent': None,
            'charge_ph': None,
            'isoelectric_point': None,
            'aromaticity': None,
            'flexibility': None,
            'gravy': None,
            'sec_struc_frac': None,
        }

        return desc_dict


def get_taxon(df, email):
    # recebe o dataframe de acessos do uniprot
    Entrez.email = email

    list_dict = []

    for index, row in df.iterrows():
        organism_id = row['Organism (ID)']

        rank_mapping = {
            'superkingdom': 'Superkingdom',
            'kingdom': 'Kingdom',
            'phylum': 'Phylum',
            'class': 'Class',
            'order': 'Order',
            'family': 'Family',
            'genus': 'Genus'
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


def main(job_id, db):
    if job_id:
        print(f"Solicitação enviada com sucesso. JobID: {job_id}")
        data = extract_id_mapping_results(job_id)
        list_dict = csv_to_dict(data)

        for data_dict in list_dict:
            protein = ProteinEntry(**data_dict)
            db.insert_data(protein)

        fasta_file = get_uniprot_fasta(job_id)

        if fasta_file:
            list_seqs = list(SeqIO.parse(fasta_file, "fasta"))

            seqs_fasta = fasta_to_dict(list_seqs)

            for seqs_data in seqs_fasta:
                seqs_entry = SeqsEntry(**seqs_data)
                db.insert_seqs_data(seqs_entry)

            for seqs_data in seqs_fasta:
                desc_calc = calc_descriptor(seqs_data['entry'],
                                            seqs_data['seq'])
                desc_entry = DescEntry(**desc_calc)
                db.insert_desc_data(desc_entry)

            print('Descritores físico-químicos calculados.')

        else:
            print("Falha ao obter dados em formato FASTA.")
            sys.exit(1)

        list_taxon = get_taxon(data, 'lspalmeira.bio@gmail.com')

        for dict_taxon in list_taxon:
            taxon = TaxonEntry(**dict_taxon)
            db.insert_taxon_data(taxon)

        print('Dados taxonômicos extraídos.')

    else:
        print("Falha ao enviar solicitação.")


if __name__ == "__main__":
    db = MongoDB()
    db.connect_to_mongodb()

    from_db = "UniProtKB_AC-ID"
    to_db = "UniProtKB"

    if 'gh32' in db.db_exists():

        entries_interpro = gh32_interpro()

        entries_interpro = list(entries_interpro)

        entries_db = db.get_entries()

        new_entries = update_db(entries_interpro, entries_db)

        if not new_entries:

            print('O banco GH32 já existe, mas não há nenhuma '
                  'atualização a ser feita.')

            print('Lista vazia: new_entries ==', new_entries)

            sys.exit(1)

        else:

            print(f'O banco de dados GH32 existe e {len(new_entries)} novas '
                  f'entradas Uniprot serão adicionadas ao banco.')

            job_id = submit_id_mapping_job(new_entries, from_db, to_db)

            main(job_id, db)

            sys.exit(1)

    else:
        print('O banco de dados GH32 será criado...')

        entries_interpro = gh32_interpro()
        entries_interpro = list(entries_interpro)
        job_id = submit_id_mapping_job(entries_interpro, from_db, to_db)

        main(job_id, db)

    db.close_connection()

    sys.exit(1)
