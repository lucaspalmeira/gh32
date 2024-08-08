FROM python:3.9-slim

WORKDIR /gh32

COPY . .

RUN apt-get update
RUN apt-get install -y wget git unzip
RUN rm -rf /var/lib/apt/lists/*
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p /opt/conda
RUN rm ~/miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"
RUN /opt/conda/bin/conda init bash
RUN /opt/conda/bin/conda env create -f environment.yml
RUN /opt/conda/bin/conda clean -afy
RUN pip install kaleido

SHELL ["conda", "run", "-n", "GH32", "/bin/bash", "-c"]

# Install CLEAN
RUN git clone https://github.com/tttianhao/CLEAN.git

WORKDIR /gh32/CLEAN/app/

RUN git clone https://github.com/facebookresearch/esm.git
RUN mkdir data/esm_data
RUN pip install gdown
RUN gdown --id 1gsxjSf2CtXzgW1XsennTr-TcvSoTSDtk
RUN unzip pretrained.zip -d data/pretrained
RUN python build.py install

WORKDIR /gh32

COPY . .

CMD ["conda", "run", "-n", "GH32", "bash", "-c", "python create_db.py && cd CLEAN/app/ && python CLEAN_infer_fasta.py --fasta_data gh32 && cd /gh32 && python ec_clean.py"]
