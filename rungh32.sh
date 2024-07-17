#!/bin/bash

# Diretório base

BASE_DIR=$(dirname $(realpath $0))

# Executar a criação/atualização do banco

python $BASE_DIR/create_db.py

# Executar o cálculo de classificação enzimática

python $BASE_DIR/CLEAN/app/CLEAN_infer_fasta.py gh32

# Executar filtragem dos dados de EC e adicioná-los ao banco

python $BASE_DIR/ec_clean.py