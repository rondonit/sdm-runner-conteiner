FROM rocker/geospatial:4.3.3

RUN apt-get update && apt-get install -y --no-install-recommends \
    awscli \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('rgbif','dplyr','tidyr','usdm','maps','sdm','argparse'), repos='https://cloud.r-project.org')"

WORKDIR /app

COPY utils.r /app/utils.r
COPY r/ /app/r/
COPY run.sh /app/run.sh

RUN chmod +x /app/run.sh

ENTRYPOINT ["/app/run.sh"]
