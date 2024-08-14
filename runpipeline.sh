#!/bin/bash

python create_db.py && python CLEAN/app/CLEAN_infer_fasta.py --fasta_data gh32 && python ec_clean.py
