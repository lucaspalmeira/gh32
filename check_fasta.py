from Bio import SeqIO
import argparse


def are_sequences_different(sequence1, sequence2):
    if len(sequence1) != len(sequence2):
        return True

    for aa1, aa2 in zip(sequence1, sequence2):
        if aa1 != aa2:
            return True

    return False


def get_identical_sequences(file_path):
    sequences = list(SeqIO.parse(file_path, 'fasta'))
    num_sequences = len(sequences)
    identical_sequences = []

    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            if not are_sequences_different(sequences[i].seq,
                                           sequences[j].seq):
                identical_sequences.append(sequences[i])
                identical_sequences.append(sequences[j])

    return identical_sequences


def save_to_multifasta(sequences, output_file):
    with open(output_file, 'w') as output_handle:
        SeqIO.write(sequences, output_handle, 'fasta')

parser = argparse.ArgumentParser(description='Conferir sequências repetidas')


parser.add_argument('arquivo_entrada',
                    type=str, help='Caminho para o arquivo de entrada')

parser.add_argument('arquivo_saida',
                    type=str, help='Caminho para o arquivo de entrada')

args = parser.parse_args()

identical_sequences = get_identical_sequences(args.arquivo_entrada)

if identical_sequences:
    print('As seguintes sequências não são diferentes '
          'em pelo menos um aminoácido:')

    save_to_multifasta(identical_sequences, args.arquivo_saida)

    print(f'As sequências idênticas foram salvas em {args.arquivo_saida}.')

else:
    print('Todas as sequências são diferentes em pelo menos um aminoácido.')
