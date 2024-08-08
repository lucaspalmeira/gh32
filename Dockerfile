FROM python:3.9-slim

WORKDIR /gh32

COPY environment.yml .

RUN apt-get update && apt-get install -y wget git unzip \
    && rm -rf /var/lib/apt/lists/* \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh \
    && bash ~/miniconda.sh -b -p /opt/conda \
    && rm ~/miniconda.sh \
    && /opt/conda/bin/conda init bash \
    && /opt/conda/bin/conda install -y -c conda-forge mamba \
    && /opt/conda/bin/mamba env create -f environment.yml \
    && /opt/conda/bin/conda clean -afy


SHELL ["conda", "run", "-n", "GH32", "/bin/bash", "-c"]

# Install CLEAN
RUN git clone https://github.com/tttianhao/CLEAN.git \
    && cd CLEAN/app/ \
    && git clone https://github.com/facebookresearch/esm.git \
    && mkdir data/esm_data \
    && pip install gdown \
    && gdown --id 1gsxjSf2CtXzgW1XsennTr-TcvSoTSDtk \
    && unzip pretrained.zip -d data/pretrained \
    && python build.py install

COPY . .

CMD ["bash", "-c", "python create_db.py && cd CLEAN/app/ && python CLEAN_infer_fasta.py --fasta_data gh32 && cd /gh32 && python ec_clean.py"]
