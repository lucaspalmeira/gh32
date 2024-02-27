import requests
import pandas as pd
from io import StringIO
from Bio import SeqIO
import os
import time


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

        return response.text

    except requests.exceptions.RequestException as error:
        print(f'Falha na solicitação: {error}')
        return None


def get_taxonomy_csv(file):
    data = pd.read_csv(file)

    base_url = ('http://bioinfo.icb.ufmg.br/cgi-bin/taxallnomy/'
                'taxallnomy_multi.pl')
    params = {
        'txid': ','.join(map(str, data['Organism (ID)'])),
        'rank': 'custom',
        'srank': 'superkingdom'
    }

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()

        data1 = pd.read_csv(StringIO(response.text), comment='#',
                            header=None, sep='\t')
        data1.columns = ['taxid', 'superkingdom']
        data1.to_csv('taxallnomy.csv', index=False)

        return data1

    except requests.exceptions.RequestException as error:
        print(f'Falha na solicitação para Taxallnomy: {error}')
        return None


def read_csv(file_csv):
    return pd.read_csv(file_csv)


def filter_inulinases(df):
    inulinase_selection = df[df['Protein names'].str.contains('inulinase|Inulinase', case=False)]
    inulinase_selection = inulinase_selection.reset_index(drop=True)
    inulinase_selection.to_csv('inulinases_clean.csv')
    return inulinase_selection


def only_inulinases(df, fasta):
    seqs_id = df['Entry']
    seqs_id = set(seqs_id)
    seqs_found = []

    with open(fasta, 'r') as arq:
        for record in SeqIO.parse(arq, 'fasta'):
            entry = record.id.split('|')[1]
            if entry in seqs_id:
                seqs_found.append(record)

    with open('only_inulinases.fasta', 'w') as out:
        SeqIO.write(seqs_found, out, 'fasta')


def find_taxid(taxid):
    df = read_csv('taxallnomy.csv')

    line = df[df['taxid'] == taxid]
    superkingdom_result = line.iloc[0]['superkingdom']

    if superkingdom_result == 'Eukaryota':
        return 'euk'
    elif superkingdom_result == 'Bacteria':
        return 'other'


def read_fasta(fasta):
    dict_seqs = {}

    with open(fasta, 'r') as arq:
        for record in SeqIO.parse(arq, 'fasta'):
            entry = record.id.split('|')[1]
            seq = record.seq
            dict_seqs[entry] = str(seq)

    return dict_seqs


def cs_pos(file):
    # CS pos = Clivage site position

    with open(file, 'r') as file:

        lines = file.readlines()
        line = lines[2]
        line_elements = line.strip().split()

        if 'pos:' in line_elements:

            pos_index = line_elements.index('pos:')
            cs_position = line_elements[pos_index + 1]

            return cs_position
        else:
            return 'Pep signal not found'


def write_fasta(filename, headers_and_sequences):
    with open(filename, 'w') as file:
        for header, sequence in headers_and_sequences.items():
            file.write(f'{header}\n{sequence}\n')


def process_signalp_results(entry, sequence, taxid, signal):

    header = ''
    processed_sequence = ''

    if type(signal) == float:
        cs_file = f'seq_{entry}.fasta'
        with open(cs_file, 'w') as file:
            file.write(f'> ID {entry}\n{sequence}')

        os.system(f'mv {cs_file} signalp/')
        os.system(f'mkdir signalp/seq_{entry}')

        os.system(
            f'signalp6 --fastafile signalp/seq_{entry}.fasta '
            f'--organism {find_taxid(taxid)} --output_dir '
            f'signalp/seq_{entry}/ --format txt')

        while not os.path.exists(f'signalp/seq_{entry}/'
                                 f'prediction_results.txt'):
            print('O arquivo prediction_results.txt ainda não existe.')
            time.sleep(1)

        cs_file = f'signalp/seq_{entry}/prediction_results.txt'
        cs_position = cs_pos(cs_file)

        if cs_position != 'Pep signal not found':
            x = 0
            y = int(cs_position[:-4])
            pep_signal = sequence[x:y]
            header = (f'> ID | {entry} | position: {x, y} | '
                      f'{find_taxid(taxid)} | peptideo sinal '
                      f'removido de acordo com predição do SignalP')
            print(f'\n> ID | {entry} | position: {x, y} | '
                  f'{find_taxid(taxid)} | peptideo sinal '
                  f'removido de acordo com predição do SignalP \n')
            processed_sequence = sequence.replace(pep_signal, '')
        else:
            with open('not_found_entries.txt', 'a') as file:
                file.write(f'> ID | {entry} | Pep signal not found '
                    f'| {find_taxid(taxid)} \n')

            print(f'\n> ID | {entry} | Pep signal not found | '
                  f'{find_taxid(taxid)}\n')

    else:
        x = int(signal[7:8])
        y = int(signal[10:12])
        pep_signal = sequence[x - 1:y - 1]
        header = (f'> ID | {entry} | position: {x, y} | '
                  f'{find_taxid(taxid)} | peptideo sinal '
                  f'removido de acordo com informações do uniprot')
        print(f'\n> ID | {entry} | position: {x, y} | '
              f'{find_taxid(taxid)} | peptideo sinal '
              f'removido de acordo com informações do uniprot\n')
        processed_sequence = sequence.replace(pep_signal, '')

    return header, processed_sequence


def remove_signalp(csv_data, sequences_dict):
    os.environ['MKL_THREADING_LAYER'] = 'INTEL'
    os.environ['MKL_SERVICE_FORCE_INTEL'] = '1'
    os.system('mkdir signalp')

    results_dict = {}

    if os.path.exists('not_found_entries.txt'):
        os.remove('not_found_entries.txt')

    for entry, signal, organism_id in zip(csv_data['Entry'], csv_data['Signal peptide'], csv_data['Organism (ID)']):
        sequence = sequences_dict[entry]
        organism_id = int(organism_id)
        header, processed_sequence = process_signalp_results(entry, sequence, organism_id, signal)
        results_dict[header] = processed_sequence

    return results_dict


def main():

    result = get_uniprot_data('inulinase')
    fasta_sequence = get_uniprot_fasta('inulinase')
    get_taxonomy = get_taxonomy_csv('inulinases.csv')

    if result is not None and fasta_sequence is not None and get_taxonomy is not None:
        inulinases_filtered = filter_inulinases(result)
        only_inulinases(inulinases_filtered, 'inulinases_uniprot.fasta')
        sequences_dict = read_fasta('only_inulinases.fasta')
        processed_results = remove_signalp(inulinases_filtered, sequences_dict)
        write_fasta('inulinases_no_signalp.fasta', processed_results)

    else:
        print('Error')


if __name__ == '__main__':
    main()
    os.system('muscle -align inulinases_no_signalp.fasta -output '
              'inu_no_sigp_aln.fasta')
    os.system('trimal -in inu_no_sigp_aln.fasta -out inu_no_sigp_trim.fasta '
              '-automated1')