#!/bin/bash

docker run -it --name gh32-pipeline -v $(pwd):/gh32_data_results gh32 /bin/bash -c 'python /gh32_data_results/app/create_db.py'