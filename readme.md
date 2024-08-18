# GH32 Database


<p align="center">
    <img src="https://github.com/lucaspalmeira/gh32/blob/master/gh32_project.png" alt="drawing"  width="800"/>
</p>

#### Clone o repositório:

```bash
git clone https://github.com/lucaspalmeira/gh32.git
```

#### Construa a imagem gh32-pipeline

```bash
cd /gh32
docker build -t gh32-pipeline .
```

#### Faça o pull da imagem do CLEAN

```bash
docker pull moleculemaker/clean-image-amd64
```

#### Faça o pull da imagem do MongoDB

```bash
docker pull mongo:latest
```

#### Crie uma rede Docker

```bash
docker network create gh32-network
```

#### Execute o conteiner mongodb

```bash
docker run -d --name mongodb --network gh32-network -p 27017:27017 mongo
```

#### Execute o pipeline

```bash
chmod +x runpipeline.sh
./runpipeline
```