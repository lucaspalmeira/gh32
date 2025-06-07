##################################################################################
# Script to extract sequences from MongoDB based on EC numbers.                  #
# This script retrieves sequences from the 'seqs_entries' collection             #
# and saves them in separate FASTA files based on their EC numbers.              #
# It uses the 'entry' field as the header and the 'seq' field as the sequence.   #
# The EC numbers are specified in the script, and sequences are saved            #
# in files named according to their EC numbers.                                  #
##################################################################################


import pymongo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

client = pymongo.MongoClient('mongodb://172.17.0.2:27017/')
db = client['gh32']
protein_collection = db['protein_entries']
seqs_collection = db['seqs_entries']

def getseq(entry):
    seq = seqs_collection.find({'entry': entry}, {'_id': 0, 'seq': 1})
    for i in seq:
        i = (list(i.values()))
        i = i[0]
        return i
    
for data in protein_collection.find({}, {'entry': 1, '_id': 0, 'ec_number':1, 'organism':1}):
    ec_number = str(data.get('ec_number', ''))
    entry = data.get('entry', '')
    organism = data.get('organism', '')
    
    if ec_number == '3.2.1.7':
        with open('3.2.1.7.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '3.2.1.80':
        with open('3.2.1.80.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '3.2.1.26':
        with open('3.2.1.26.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '3.2.1.153':
        with open('3.2.1.153.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '3.2.1.154':
        with open('3.2.1.154.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '3.2.1.64':
        with open('3.2.1.64.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '3.2.1.65':
        with open('3.2.1.65.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '4.2.2.16':
        with open('4.2.2.16.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '3.2.1.-':
        with open('3.2.1.-.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '2.4.1.99':
        with open('2.4.1.99.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '2.4.1.243':
        with open('2.4.1.243.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '2.4.1.100':
        with open('2.4.1.100.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '2.4.1.10':
        with open('2.4.1.10.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')
            
    elif ec_number == '2.4.1.-':
        with open('2.4.1.-.fasta', 'a') as file:
            file.write(f'>{entry} {organism} {ec_number}\n{getseq(entry)}\n')

print('Concluído.')