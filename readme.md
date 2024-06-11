### Using GH32

Clone o repositório
```bash
git clone https://github.com/lucaspalmeira/gh32.git
```

Crie o ambiente

```bash
conda create -n GH32 python==3.11.9 -y
conda activate GH32
```

Instale as bibliotecas
```bash
pip install -r requirements.txt
```

Start MongoDB
```bash
sudo systemctl start mongod
```

Criar banco de dados de inulinases fúngicas

Execute
```bash
python create_db.py
```


### Calcular número de comissão enzimática (EC number) com o CLEAN

Clone o repositório.
```bash
git clone https://github.com/tttianhao/CLEAN.git
```

Navegue até o diretório CLEAN/
```bash
cd CLEAN/
```

Criar variável de ambiente com a versão 3.11 do Python
```bash
conda create --name clean_env python=3.11
```

Ativar a veriável de ambiente
```bash
conda activate clan_env
```

Instale as bibliotecas
```bash
pip install -r requirements.txt
```

Instale o Pytorch
```bash
conda install pytorch torchvision torchaudio cpuonly -c pytorch
```

Compile o CLEAN
```bash
python build.py install
```

Clone o repositório ESM e crie a pasta data/esm_data
```bash
git clone https://github.com/facebookresearch/esm.git
mkdir data/esm_data
```

Preparando seus dados
Mova seu arquivo fasta (ex.: query.fasta) para a pasta 'data/inputs'.

Exemplo:

```bash
>Sequence1
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVGA
>Sequence2
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
```

Execute a predição na pasta CLEAN/app
```bash
python CLEAN_infer_fasta.py --fasta_data query
```

Os resultados estarão na pasta 'results' em formato CSV
```bash
results/query_maxsep.csv
```

Desativar a vaariável de ambiente
```bash
conda deactivate
```

### Alteração e inclusão do EC number na collection 'protein_entries'

Os resultados gerados pelo CLEAN serão lidos ('gh32_maxsep.csv') 
e cada EC será alterado ou adiocionado de acordo com a sua respectiva 
entrada ('entry').

Execute

```bash
python ec_clean.py
```

Um arquivo de saída ('gh32_maxsep_clean.csv') será criado contendo 
os dados limpos.