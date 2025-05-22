FROM python:3.11.0-slim-buster

WORKDIR /app

RUN apt-get update && \ 
    apt-get install -y build-essential \
    wget \
    perl \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

# MAFFT
# download
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.526-1_amd64.deb -O mafft_7.526-1_amd64.deb
# install
RUN dpkg -i mafft_7.526-1_amd64.deb
# cleanup
RUN rm mafft_7.526-1_amd64.deb

COPY ./requirements.txt /app/requirements.txt

RUN pip install -r requirements.txt

COPY ./ /app/

RUN chmod +x /app/h1n1.sh /app/h5n1.sh
