#!/bin/bash

# uso do operador ; para executar o código mesmo que o código de saída anterior seja 1
python create_db.py; python CLEAN/app/CLEAN_infer_fasta.py --fasta_data gh32; python ec_clean.py
