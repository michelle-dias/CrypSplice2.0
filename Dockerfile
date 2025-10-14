FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential git libssl-dev libffi-dev python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy all files including the 'crypsplice' package
COPY . .

# Upgrade pip and install Python dependencies
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN pip install .

# Set entrypoint to the CLI (or module)
ENTRYPOINT ["crypsplice"]
