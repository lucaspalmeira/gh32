#!/bin/bash

docker run -it --network gh32-network -v $(pwd):/gh32_data_results gh32-pipeline /bin/bash -c 'python /gh32_data_results/app/create_db.py'