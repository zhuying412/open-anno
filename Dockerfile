from docker.io/staphb/bcftools:1.16
RUN apt update -y && apt install -y tabix
COPY bin/openanno /usr/local/bin/openanno
RUN chmod a+x /usr/local/bin/openanno
