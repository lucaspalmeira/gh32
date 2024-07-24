#!/bin/bash

# Diretório base

BASE_DIR=$(dirname $(realpath $0))

# Executar a criação/atualização do banco

python $BASE_DIR/create_db.py

# Executar o cálculo de classificação enzimática

cd $BASE_DIR/CLEAN/app/

python CLEAN_infer_fasta.py --fasta_data gh32

# Executar filtragem dos dados de EC e adicioná-los ao banco

cd $BASE_DIR

python $BASE_DIR/ec_clean.py