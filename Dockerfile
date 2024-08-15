FROM python:3.9-slim

WORKDIR /gh32

COPY . .

RUN apt-get update && apt-get install openssl
RUN rm -rf /var/lib/apt/lists/*
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN pip install kaleido

RUN python CLEAN/app/build.py install

# CMD ["tail", "-f", "/dev/null"]
CMD ["./runpipeline.sh"]