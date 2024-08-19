#!/bin/bash

docker run -it --network gh32-network -v $(pwd):/gh32_data_results gh32-pipeline /bin/bash -c 'python /gh32_data_results/app/create_db.py'; \
docker run -it --network gh32-network -v $(pwd):/gh32_data_results moleculemaker/clean-image-amd64 /bin/bash -c 'echo Starting Execution && python /gh32_data_results/CLEAN_infer_fasta.py --fasta_data /gh32_data_results/gh32'; \
docker run -it --network gh32-network -v $(pwd):/gh32_data_results gh32-pipeline /bin/bash -c 'python /gh32_data_results/app/ec_clean.py'

echo "Concluído."