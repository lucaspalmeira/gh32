#!/bin/bash

docker run -it --network gh32-network -v $(pwd):/gh32 gh32-pipeline /bin/bash -c 'python /gh32/app/create_db.py'; \
cd app; docker run -it -v $(pwd):/app/data/inputs -v $(pwd)/results/inputs:/app/results/inputs/ moleculemaker/clean-image-amd64 /bin/bash -c 'echo Starting Execution && python $(pwd)/CLEAN_infer_fasta.py --fasta_data gh32'; \
cd ..; docker run -it --network gh32-network -v $(pwd)/app/results/inputs:/results/inputs -v $(pwd):/gh32 gh32-pipeline /bin/bash -c 'python /gh32/app/ec_clean.py

echo "Concluído."