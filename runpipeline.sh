#!/bin/bash

docker run -d -it --name mongodb mongo

docker run -it --name gh32-pipeline -v $(pwd):/gh32_data_results GH32 /bin/bash -c 'python /app/create_db.py'