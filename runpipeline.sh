#!/bin/bash

python create_db.py

docker run -it \
  -v $(pwd):/root/.cache/torch/hub/checkpoints \
  -v $(pwd)/CLEAN/app/data/inputs:/app/data/inputs \
  -v $(pwd)/results:/app/results/inputs \
  clean-image \
  /bin/bash -c 'echo Starting Execution && python /app/CLEAN_infer_fasta.py --fasta_data gh32'

python ec_clean.py
