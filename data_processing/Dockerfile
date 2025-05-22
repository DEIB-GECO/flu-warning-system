FROM python:3.11.0-slim-buster

WORKDIR /app

RUN apt-get update && apt-get install -y build-essential
RUN pip install --upgrade pip

COPY ./requirements.txt /app/requirements.txt

RUN pip install -r requirements.txt

COPY ./ /app/

ENTRYPOINT ["python", "-u", "warnings_analysis.py"]