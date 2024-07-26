### Using GH32


<p align="center">
    <img src="https://github.com/lucaspalmeira/gh32/blob/master/gh32_project.png" alt="drawing"  width="600"/>
</p>


#### Preparando o ambiente

Clone o repositório:
```bash
git clone https://github.com/lucaspalmeira/gh32.git
```

Crie o ambiente:
```bash
conda env create -f environment.yml
conda activate GH32
```

Clone o repositório:
```bash
git clone https://github.com/tttianhao/CLEAN.git
```

Navegue até o diretório '**CLEAN/app/**':
```bash
cd CLEAN/app/
```

Clone o repositório ESM e crie o diretório '**data/esm_data**':
```bash
git clone https://github.com/facebookresearch/esm.git
mkdir data/esm_data
```

Crie o repositório '**data/pretrained**':
```bash
mkdir data/pretrained
```

<p>Faça o Download dos arquivos pré-treinados conforme as instruções de instalação do <a href="https://github.com/tttianhao/CLEAN?tab=readme-ov-file#1-install">CLEAN</a>.</p>

Os arquivos pré-treinados deverão estar em '**data/pretrained**'

Compile o CLEAN (certifique-se de estar em '**CLEAN/app**'):
```bash
python build.py install
```

### Instruções para instalar o MongoDB
<p>Para instalar o MongoDB, siga o tutorial de instalação presente na documentação <a href="https://www.mongodb.com/pt-br/docs/manual/installation/">aqui</a>.</p>

Inicie o MongoDB:
```bash
sudo systemctl start mongod
```

### Run
```bash
bash rungh32.sh
```
ou
```bash
chmod +x rungh32.sh
./rungh32.sh
```

### Backup do banco

Para realizar o backup do banco GH32, execute:
```bash
mongodump --db gh32 --out /gh32/backup/
```

Para restaurar o banco de dados, execute:
```bash
mongorestore --db gh32 /gh32/backup/
```

Desative a variável de ambiente:
```bash
conda deactivate
```