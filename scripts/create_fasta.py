################################################################################
# Script to create a FASTA file from MongoDB entries                           #
# This script retrieves sequences from the 'seqs_entries' collection           #
# and saves them in a FASTA format file named 'gh32_all_fasta.fasta'.          #
# It uses the 'entry' field as the header and the 'seq' field as the sequence. #
################################################################################


import pymongo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

client = pymongo.MongoClient('mongodb://172.17.0.2:27017/')

db = client['gh32']

seqs_collection = db['seqs_entries']

with open('gh32_all_fasta.fasta', 'w') as file:
    for data in seqs_collection.find({}, {'entry': 1, '_id': 0, 'seq':1}):
        data = list(data.values())
        data = f'>{data[0]}\n{data[1]}'
        file.write(data+'\n')