FROM python:3.13-slim

WORKDIR /app

RUN apt-get update && apt-get install -y libmariadb-dev gcc && rm -rf /var/lib/apt/lists/*
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8501

CMD ["streamlit", "run", "main.py"] 