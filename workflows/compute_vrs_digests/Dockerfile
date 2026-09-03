FROM eugene75/vrs-vcf-annotator-grch38:2.3.1

# Copy our VRS digest computation script
COPY compute_vrs_digests.py /vrs-python/compute_vrs_digests.py
RUN chmod +x /vrs-python/compute_vrs_digests.py

WORKDIR /data

# Set Python to not buffer output
ENV PYTHONUNBUFFERED=1

# Remove the vcf-annotator entrypoint so we can run our script
ENTRYPOINT []
