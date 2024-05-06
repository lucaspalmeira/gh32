import matplotlib.pyplot as plt
from Bio import SeqIO

"""
Histograma para verificar a presença de sequências que possívelmente são 
outliers devido ao seu comprimento.
"""

# lista que contém cada uma das sequências
list_seqs = [record.seq for record in SeqIO.parse('gh32.fasta',
                                                  'fasta')]

# lista que contém o comprimento de aminoácidos de cada uma das sequências
length_seqs = [len(seq) for seq in list_seqs]

plt.figure(figsize=(20, 10))

# 'bins' é o número de intervalos (3000 intervalos)
hist_vals = plt.hist(length_seqs, bins=3000, color='skyblue',
                     edgecolor='black', linewidth=0.5)

plt.title('Sequence Lengths for 3148 Protein Sequences')
plt.xlabel('Length of sequences')
plt.ylabel('Number of sequences')

# marcadores a cada 100 unidades ao longo do eixo x
# criando assim intervalos de 100 em 100
plt.xticks(range(0, max(length_seqs) + 100, 100), fontsize=7)

# hist_vals[0] contém as alturas das barras do histograma
# ou seja, o número de sequências em cada intervalo
max_y = int(max(hist_vals[0]))
plt.yticks(range(0, max_y, 5), fontsize=7)

# linhas de grade apenas para o eixo y
plt.grid(axis='y')
plt.show()
