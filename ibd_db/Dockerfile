### Builder stage: install dependencies using pdm and create a project .venv
FROM python:3.13 AS builder

# Put the project venv inside the project directory so we can copy it later
ENV PDM_VENV_IN_PROJECT=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       gcc \
       libbz2-dev \
       liblzma-dev \
       libzstd-dev \
       cargo \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy lockfile and project metadata first to leverage layer caching
COPY pyproject.toml pdm.lock README.md /app/

# Install pdm (will be used to install packages according to pdm.lock)
RUN python -m pip install --upgrade pip \
    && pip install pdm

# Copy source code
COPY src/ /app/src/


# Install production dependencies from the lock into .venv in the project
RUN pdm install --prod


### Final stage: slim runtime image with the created virtualenv
FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PATH="/app/venv/bin:$PATH"

WORKDIR /app

# Copy the virtualenv created by pdm in the builder stage
COPY --from=builder /app/.venv /app/.venv

# Copy application source
COPY src/ /app/src/
COPY scripts/ /app/scripts/

# # Ensure the package source can be imported
# ENV PYTHONPATH=/app/src

ENV PATH="/app/.venv/bin:$PATH"

EXPOSE 8080

# # Run the FastAPI app via uvicorn. Module: `ibd_db.main`, app var: `app`.
CMD ["fastapi", "run", "./src/ibd_db/main.py", "--port", "8080"]
