#!/bin/bash

# Executar a criação/atualização do banco

python create_db.py

# Executar o cálculo de classificação enzimática

python /CLEAN/app/CLEAN_infer_fasta.py gh32

# Executar filtragem dos dados de EC e adicioná-los ao banco

python ec_clean.py