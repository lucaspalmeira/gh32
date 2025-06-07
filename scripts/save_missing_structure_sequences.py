#####################################################
# Script para exportar em formato FASTA             #
# as sequências de proteínas cuja entrada           #
# na coleção 'protein_entries' possui os campos     #
# 'alphafold_db' E 'pdb' vazios.                    #
#                                                   #
# Para cada entrada que satisfaz esta condição,     #
# a sequência correspondente é buscada na coleção   #
# 'seqs_entries' (campo 'seq') e salva no arquivo   #
# 'proteins_missing_alphafold_pdb.fasta'.           #
#                                                   #
# O campo comum entre as collections é 'entry',     #
# que é utilizado como cabeçalho no arquivo FASTA.  #
#                                                   #
# Condições:                                        #
# alphafold_db vazio E pdb vazio → salva no FASTA   #
# Em qualquer outro caso → não salva.               #
#####################################################


import pymongo

client = pymongo.MongoClient('mongodb://172.17.0.2:27017/')
db = client['gh32']

protein_collection = db['protein_entries']
seqs_collection = db['seqs_entries']

# alphafold_db vazio E pdb vazio
query = {
    'alphafold_db': {
        '$in': [None, '', [], float('nan')]
    },
    'pdb': {
        '$in': [None, '', [], float('nan')]
    }
}

with open('proteins_missing_alphafold_pdb.fasta', 'w') as fasta_file:
    proteins_cursor = protein_collection.find(query, {'entry': 1, '_id': 0})
    
    for protein in proteins_cursor:
        entry_id = protein['entry']
        
        # Busca a sequência correspondente na seqs_entries
        seq_doc = seqs_collection.find_one({'entry': entry_id}, {'seq': 1, '_id': 0})
        
        if seq_doc:
            sequence = seq_doc['seq']
            
            # Escreve em formato FASTA
            fasta_file.write(f'>{entry_id}\n')
            
            # Quebra a sequência em linhas de no máximo 60 caracteres (padrão fasta)
            for i in range(0, len(sequence), 60):
                fasta_file.write(sequence[i:i+60] + '\n')
            
            print(f'Saved entry {entry_id} to FASTA')
        else:
            print(f'Seq not found for entry {entry_id}')