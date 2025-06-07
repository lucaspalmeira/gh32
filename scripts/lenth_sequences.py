########################################################################
#                                                                      #
# This script reads a FASTA file containing protein sequences,         #
# and separates the sequences into two files based on their length:    #
# one file for sequences shorter than 1000 amino acids,                #
# and another for sequences longer than or equal to 1000 amino acids.  #
#                                                                      #
########################################################################

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

for seq_record in SeqIO.parse("proteins_missing_alphafold_pdb.fasta", "fasta"):
    
    my_record = SeqRecord(Seq(seq_record.seq),id=seq_record.id, description=f'lenth {len(seq_record)}')
    
    if len(seq_record) < 1000:
        with open('select_sequences_under_1k_aa.fasta', 'a') as file:
            SeqIO.write(my_record, file, 'fasta')
    else:
        with open('select_sequences_over_1k_aa.fasta', 'a') as file:
            SeqIO.write(my_record, file, 'fasta')