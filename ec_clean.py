import pandas as pd
import plotly.express as px
from DataBase import MongoDB
import sys
import os
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PATH_CLEAN = BASE_DIR + '/CLEAN/app/'
sys.path.append(os.path.join(PATH_CLEAN, 'src'))
print(BASE_DIR)
print(PATH_CLEAN)
from CLEAN.utils import *
from CLEAN.infer import infer_maxsep


# extrair o EC number
def find_ec(title):
    start = title.find('EC:')
    if start == -1:
        return ''
    end = title.find('/', start)
    if end == -1:
        end = len(title)
    ec = title[start + 3:end].strip()
    return ec


# extrair o ID
def find_id(title):
    parts = title.split('|')
    if len(parts) > 1:
        id_seq = parts[1]
        return id_seq.strip()
    return ''


def find_name_enzyme(ec):
    enzyme_names = {
        '3.2.1.7': 'Endo-inulinase (EC 3.2.1.7)',
        '3.2.1.80': 'fructan beta-fructosidase (EC 3.2.1.80)',
        '3.2.1.26': 'β-fructofuranosidase / invertase (EC 3.2.1.26)',
        '3.2.1.153': 'fructan 1-exohydrolase (1-FEH) / 1-exo-inulinase / '
                     'fructan β-2,1-fructosidase (EC 3.2.1.153)',
        '3.2.1.154': 'fructan 6-exohydrolase (6-FEH) / 6-exo-levanase / '
                     'fructan β-2,6-fructosidase (EC 3.2.1.154)',
        '3.2.1.64': 'β-2,6-fructan 6-levanbiohydrolase (EC 3.2.1.64)',
        '3.2.1.65': 'Endo-levanase / β-2,6-fructanase (EC 3.2.1.65)',
        '4.2.2.16': 'Levan fructotransferase (DFA-IV-forming) (EC 4.2.2.16)',
        '3.2.1.-': '6-kestose β-2,6-fructosidase (EC 3.2.1.-)',
        '2.4.1.99': 'sucrose:sucrose 1F-fructosyltransferase (EC 2.4.1.99)',
        '2.4.1.243': 'fructan:fructan 6-fructosyltransferase (6-FFT) '
                     '(EC 2.4.1.243)',
        '2.4.1.100': 'fructan:fructan 1-fructosyltransferase (1-FFT) '
                     '(EC 2.4.1.100)',
        '2.4.1.10': 'levansucrase / sucrose:fructan 6-fructosyltransferase '
                    '(6-SFT) (EC 2.4.1.10)',
        '2.4.1.-': 'sucrose:sucrose 6-fructosyltransferase / '
                   'cycloinulo-oligosaccharide fructanotransferase / '
                   'Sucrose 6 acetate β-2,1-fructosyltransferase (EC 2.4.1.-)'
    }

    return enzyme_names.get(ec, 'Not in Glycoside Hydrolase Family 32')


def process(file):
    data_ec = []

    enzyme_count = {}

    with open(file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            EC = find_ec(line)
            ID = find_id(line)
            NAME = find_name_enzyme(EC)
            data_ec.append({'entry': ID, 'ec_number': EC,
                            'enzyme_name': NAME})

            if NAME in enzyme_count:
                enzyme_count[NAME] += 1
            else:
                enzyme_count[NAME] = 1

    return data_ec, enzyme_count


def main():

    data_ec, enzyme_count = process('gh32_maxsep.csv')

    df = pd.DataFrame(data_ec)
    df.to_csv('gh32_maxsep_clean.csv', index=False)

    db = MongoDB()
    db.connect_to_mongodb()

    # adição dos dados preditivos de EC
    db.update_ec_number(df)

    # remoção de entradas de sequêncais que não pertecem a GH32
    removed_entries = df[df['enzyme_name'] == 'Not in Glycoside Hydrolase Family 32']['entry'].tolist()
    db.remove_enzymes_not_in_gh32(removed_entries)
    db.close_connection()

    enzyme_df = pd.DataFrame(list(enzyme_count.items()),
                             columns=['NAME', 'COUNT'])

    title = 'Distribuição de Enzimas na Família 32 das Glicosídeo Hidrolases'
    fig = px.pie(enzyme_df, values='COUNT', names='NAME', width=1200,
                 height=800, title=title)

    fig.update_traces(textinfo='percent+value')

    fig.write_image('distribuicao_enzimas_gh32.png')


if __name__ == '__main__':
    train_data = 'split100'
    test_data =  'gh32'
    # Converting fasta to dummy csv file, will delete after inference
    prepare_infer_fasta(test_data)

    # Inferred results is in
    # results/[args.fasta_data].csv
    infer_maxsep(train_data, test_data, report_metrics=False, pretrained=True,
                 gmm=os.path.join(BASE_DIR, 'data', 'pretrained',
                                  'gmm_ensumble.pkl'))

    # Removing dummy csv file
    os.remove("data/"+ test_data +'.csv')
    main()
