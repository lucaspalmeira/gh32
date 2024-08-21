#!/bin/bash

docker run -it --network gh32-network -v $(pwd):/gh32 gh32-pipeline /bin/bash -c 'python /gh32/app/create_db.py'; \
docker run -it --network gh32-network -v $(pwd):/gh32 moleculemaker/clean-image-amd64 /bin/bash -c 'echo Starting Execution && python $(pwd)/CLEAN_infer_fasta.py --fasta_data gh32/app/gh32'; \
docker run -it --network gh32-network -v $(pwd):/gh32 gh32-pipeline /bin/bash -c 'python /gh32/app/ec_clean.py'

echo "Concluído."