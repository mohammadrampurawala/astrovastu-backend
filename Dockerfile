# ===========================================
# AstroVastu Backend Dockerfile
# ===========================================
FROM python:3.11-slim

# Prevent Python from buffering stdout/stderr
ENV PYTHONUNBUFFERED=1

# Set work directory
WORKDIR /app

# Install system dependencies (for timezonefinder, numpy, etc.)
RUN apt-get update && apt-get install -y \
    build-essential \
    libgeos-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy Python dependency list
COPY requirements.txt .

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the backend code
COPY . .

# Expose FastAPI port
EXPOSE 8000

# Default run command
CMD ["uvicorn", "astro_service_with_dasha:app", "--host", "0.0.0.0", "--port", "8000"]
