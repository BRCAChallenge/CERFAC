# Single-stage: SeqRepo + PostgreSQL + UTA (from database dump)
FROM python:3.11-slim

ARG ASSEMBLY="GRCh38"

# Install dependencies including build tools for psycopg2
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    build-essential curl git libpq-dev libpq5 postgresql postgresql-contrib unzip \
    && rm -rf /var/lib/apt/lists/*

# Setup Python venv and VRS libraries
RUN python3 -m venv /vrs-env && \
    /vrs-env/bin/pip install --upgrade setuptools wheel pip && \
    /vrs-env/bin/pip install 'ga4gh-vrs==2.3.1' 'biocommons.seqrepo>=0.6.0' 'hgvs>=1.5.0' 'requests' 'pandas'

# Download SeqRepo
RUN curl -L -o /tmp/seqrepo.zip \
    https://github.com/ehclark/vrs-in-a-box/releases/download/seqrepofiles/seqrepo-${ASSEMBLY}.zip && \
    unzip /tmp/seqrepo.zip -d / && \
    rm /tmp/seqrepo.zip

# Create postgres user and initialize database cluster
RUN useradd -m -s /bin/bash postgres 2>/dev/null || true && \
    mkdir -p /pgdata && \
    chown postgres:postgres /pgdata && \
    chmod 0700 /pgdata && \
    su - postgres -c "/usr/lib/postgresql/*/bin/initdb -D /pgdata" && \
    sed -i 's/^local   all             all                                     peer/local   all             all                                     md5/' /pgdata/pg_hba.conf && \
    sed -i 's/^host    all             all             127.0.0.1\/32            scram-sha-256/host    all             all             127.0.0.1\/32            md5/' /pgdata/pg_hba.conf && \
    echo "host    all             all             ::1/128                 md5" >> /pgdata/pg_hba.conf && \
    /etc/init.d/postgresql start && sleep 2 && \
    su - postgres -c "psql -c \"ALTER USER postgres WITH PASSWORD 'postgres';\"" && \
    /etc/init.d/postgresql stop

# Copy pre-built UTA database dump
COPY uta_20241220.pgd /tmp/uta.pgd

# Restore UTA database from dump (text format SQL dump)
RUN /etc/init.d/postgresql start || su - postgres -c "/usr/lib/postgresql/*/bin/pg_ctl -D /pgdata start" && \
    sleep 2 && \
    su - postgres -c "/usr/lib/postgresql/*/bin/psql -c 'CREATE DATABASE uta;'" 2>&1 || true && \
    cat /tmp/uta.pgd | su - postgres -c "/usr/lib/postgresql/*/bin/psql -d uta" 2>&1 && \
    su - postgres -c "/usr/lib/postgresql/*/bin/psql -d uta -c 'ALTER DATABASE uta SET search_path TO uta_20241220, public;'" 2>&1 && \
    su - postgres -c "/usr/lib/postgresql/*/bin/pg_ctl -D /pgdata stop" || true && \
    rm /tmp/uta.pgd && \
    echo "✓ UTA database loaded from dump"

# Set runtime environment
ENV VIRTUAL_ENV=/vrs-env PATH=/vrs-env/bin:$PATH \
    GA4GH_VRS_DATAPROXY_URI=seqrepo+file:///seqrepo-${ASSEMBLY}/master \
    UTA_DB_URL=postgresql://postgres:postgres@localhost/uta/uta_20241220 \
    PYTHONUNBUFFERED=1 PGDATA=/pgdata

# Copy VRS script
COPY compute_vrs_digests.py /usr/local/bin/compute_vrs_digests.py
RUN chmod +x /usr/local/bin/compute_vrs_digests.py

# Entrypoint: start PostgreSQL, then run script with proper argument handling
RUN printf '#!/bin/bash\necho "Starting PostgreSQL + UTA..."\n/etc/init.d/postgresql start 2>&1 || true\nsleep 2\nargs=()\nfor arg in "$@"; do\n  [ "$arg" != "--" ] && args+=("$arg")\ndone\nexec python3 /usr/local/bin/compute_vrs_digests.py "${args[@]}"\n' > /entrypoint.sh && chmod +x /entrypoint.sh

WORKDIR /data
ENTRYPOINT ["/entrypoint.sh"]
CMD ["python3", "/usr/local/bin/compute_vrs_digests.py"]
