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

CMD ["tail", "-f", "/dev/null"]