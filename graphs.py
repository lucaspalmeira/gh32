import pymongo
import plotly.express as px
import pandas as pd

def verify_outliers():
    """
    Histograma para verificar a presença de sequências que possivelmente são
    outliers devido ao seu comprimento.
    """

    # Conectar ao MongoDB
    client = pymongo.MongoClient("mongodb://localhost:27017/")
    db = client["gh32"]  # Substitua pelo nome do seu banco de dados
    collection = db["seqs_entries"]

    # Buscar todas as sequências da coleção
    sequences = collection.find({}, {"seq": 1, "entry": 1, "_id": 0})

    # Lista que contém o comprimento de aminoácidos de cada uma das sequências
    length_seqs = [(record["entry"], len(record["seq"])) for record in sequences]

    # Criar um DataFrame para facilitar a plotagem com Plotly
    df = pd.DataFrame(length_seqs, columns=['Entry', 'Sequence Length'])

    # Criar o histograma com Plotly
    fig = px.histogram(df, x='Sequence Length', nbins=3000,
                       title='Sequence Lengths for Protein Sequences',
                       labels={'Sequence Length': 'Length of sequences'},
                       template='plotly_white')

    # Atualizar o layout do gráfico
    fig.update_layout(
        xaxis=dict(
            tickmode='linear',
            tick0=0,
            dtick=100,
            title='Length of sequences'
        ),
        yaxis=dict(
            title='Number of sequences'
        ),
        bargap=0.2,
        bargroupgap=0.1
    )

    # Mostrar o gráfico
    fig.show()
    fig.write_image("sequence_lengths_histogram.png")

# Chamar a função
verify_outliers()
